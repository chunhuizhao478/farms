"""
mpirun -np 8 python correlated_vonkarman_2D.py \
  --nx 1024 --ny 256 \
  --x-min -20000 --x-max 20000 \
  --y-min -1000  --y-max 1000 \
  --taper-width 1000 \
  --nu 0.5 --corr-len 200 \
  --mean 0.1 --std 0.025 \
  --seed 1234 \
  --exo static_solve_out.e \
  --out mapped_vonkarman_field.csv
"""

import argparse
import numpy as np
from mpi4py import MPI
import netCDF4
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import rfftn, irfftn
from scipy.special import gamma as gamma_function
import warnings


def matern_spectral_density_2d(kx, ky, nu, corr_len):
    """
    Spectral density for isotropic Matérn (von Kármán) in 2D up to a constant factor.
    S(k) ∝ (k0^2 + |k|^2)^(-(nu + d/2)), with d=2 and k0 = sqrt(2*nu)/corr_len.
    We only need relative shape for sampling; scale later to desired variance.
    """
    k0 = np.sqrt(2.0 * nu) / corr_len
    k2 = kx**2 + ky**2
    alpha = nu + 1.0  # d/2 = 1 for 2D
    return (k0**2 + k2) ** (-alpha)


def generate_vonkarman_field_fft(nx, ny, x_min, x_max, y_min, y_max, nu, corr_len, rng):
    """
    Generate a zero-mean unit-variance Gaussian random field on a regular grid
    with Matérn/von Kármán covariance using spectral synthesis.

    Returns (x, y, field) with shapes (nx,), (ny,), (ny, nx)
    """
    # Grids
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)

    # Frequencies for rFFT along x (last axis) and full FFT along y
    Lx = x_max - x_min
    Ly = y_max - y_min
    dkx = 2 * np.pi / Lx
    dky = 2 * np.pi / Ly

    ky = np.fft.fftfreq(ny, d=1.0 / ny) * dky  # size ny, symmetric
    kx = np.fft.rfftfreq(nx, d=1.0 / nx) * dkx  # size nx//2 + 1

    KY, KX = np.meshgrid(ky, kx, indexing="ij")  # shapes (ny, nx//2+1)
    S = matern_spectral_density_2d(KX, KY, nu, corr_len)

    # Random complex spectrum with Hermitian symmetry satisfied by rFFT
    # Real and imag parts ~ N(0, 0.5) so that variance per complex sample is 1
    real_part = rng.normal(scale=1.0, size=S.shape)
    imag_part = rng.normal(scale=1.0, size=S.shape)
    Z = (real_part + 1j * imag_part) * np.sqrt(S)

    # Enforce purely real on kx=0 and Nyquist lines per rFFT conventions
    Z[:, 0] = rng.normal(size=Z[:, 0].shape) * np.sqrt(S[:, 0])
    if nx % 2 == 0:
        Z[:, -1] = rng.normal(size=Z[:, -1].shape) * np.sqrt(S[:, -1])

    field = irfftn(Z, s=(ny, nx))

    # Normalize to unit variance (numerical scaling can vary)
    field = field - np.mean(field)
    std = np.std(field)
    if std > 0:
        field = field / std
    return x, y, field


def parallel_interpolate_and_write(
    comm,
    x,
    y,
    field,
    exo_path,
    out_csv,
    plot=False,
    mean_val=0.1,
    std_val=0.025,
    plot_prefix="mapped_vonkarman",
    oob_mode: str = "nearest",
    x_min_bound: float = -16000,
    x_max_bound: float = 16000,
    y_min_bound: float = -2000,
    y_max_bound: float = 2000,
    taper_width: float = 0.0,
    exclude_rects=None,
):
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Open Exodus/NetCDF on rank 0 and scatter coordinates
    if rank == 0:
        nc = netCDF4.Dataset(exo_path)
        x_coord = np.array(nc.variables["coordx"][:])
        y_coord = np.array(nc.variables["coordy"][:])
        nc.close()
        N = x_coord.size
    else:
        x_coord = None
        y_coord = None
        N = None

    N = comm.bcast(N, root=0)

    # Compute scatter counts and displacements
    counts = np.full(size, N // size, dtype=int)
    counts[: N % size] += 1
    displs = np.concatenate(([0], np.cumsum(counts[:-1])))

    # Create recv buffers
    recv_x = np.empty(counts[rank], dtype=float)
    recv_y = np.empty(counts[rank], dtype=float)

    # Scatterv coordinates
    if rank == 0:
        comm.Scatterv([x_coord, counts, displs, MPI.DOUBLE], recv_x, root=0)
        comm.Scatterv([y_coord, counts, displs, MPI.DOUBLE], recv_y, root=0)
    else:
        comm.Scatterv([None, counts, displs, MPI.DOUBLE], recv_x, root=0)
        comm.Scatterv([None, counts, displs, MPI.DOUBLE], recv_y, root=0)

    # Build interpolators on all ranks (read-only, small cost vs comms)
    interp_linear = RegularGridInterpolator(
        (y, x), field, bounds_error=False, fill_value=np.nan
    )
    interp_nearest = None
    if oob_mode == "nearest":
        # fill_value=None ensures nearest neighbor is used for OOB instead of NaN
        interp_nearest = RegularGridInterpolator(
            (y, x), field, method="nearest", bounds_error=False, fill_value=None
        )

    # Interpolate local chunk with OOB handling
    pts = np.column_stack((recv_y, recv_x))
    # in-bounds mask
    ib = (recv_x >= x[0]) & (recv_x <= x[-1]) & (recv_y >= y[0]) & (recv_y <= y[-1])

    if oob_mode == "nan":
        local_vals = interp_linear(pts)
    elif oob_mode == "nearest":
        local_vals = np.empty_like(recv_x)
        # in-bounds via linear
        if np.any(ib):
            local_vals[ib] = interp_linear(pts[ib])
        # out-of-bounds via nearest
        if np.any(~ib):
            local_vals[~ib] = interp_nearest(pts[~ib])
    elif oob_mode == "clip":
        pts_clipped = np.column_stack(
            (np.clip(recv_y, y[0], y[-1]), np.clip(recv_x, x[0], x[-1]))
        )
        local_vals = interp_linear(pts_clipped)
    else:
        raise ValueError("oob_mode must be one of: 'nan', 'nearest', 'clip'")

    # Map to damage variable
    mapped_local = mean_val + std_val * local_vals

    # Apply exclusion rectangles (set to zero) BEFORE taper so zeros stay zero.
    excl_zeroed = 0
    if exclude_rects:
        for ex_xmin, ex_xmax, ex_ymin, ex_ymax in exclude_rects:
            m_ex = (
                (recv_x >= ex_xmin)
                & (recv_x <= ex_xmax)
                & (recv_y >= ex_ymin)
                & (recv_y <= ex_ymax)
            )
            if np.any(m_ex):
                mapped_local[m_ex] = 0.0
                excl_zeroed += int(np.count_nonzero(m_ex))

    # Apply optional cosine taper outside the user box; width in meters.
    # If taper_width <= 0: hard cutoff to zero outside the box.
    if taper_width and taper_width > 0:
        # 1D cosine taper factor along x
        fx = np.ones_like(recv_x, dtype=float)
        # Right side (x_max .. x_max + w)
        mask_r = (recv_x > x_max_bound) & (recv_x <= x_max_bound + taper_width)
        fx[mask_r] = 0.5 * (
            1.0 + np.cos(np.pi * (recv_x[mask_r] - x_max_bound) / taper_width)
        )
        # Left side (x_min - w .. x_min)
        mask_l = (recv_x < x_min_bound) & (recv_x >= x_min_bound - taper_width)
        fx[mask_l] = 0.5 * (
            1.0 + np.cos(np.pi * (x_min_bound - recv_x[mask_l]) / taper_width)
        )
        # Beyond taper -> 0
        fx[recv_x > x_max_bound + taper_width] = 0.0
        fx[recv_x < x_min_bound - taper_width] = 0.0

        # 1D cosine taper factor along y
        fy = np.ones_like(recv_y, dtype=float)
        mask_t = (recv_y > y_max_bound) & (recv_y <= y_max_bound + taper_width)
        fy[mask_t] = 0.5 * (
            1.0 + np.cos(np.pi * (recv_y[mask_t] - y_max_bound) / taper_width)
        )
        mask_b = (recv_y < y_min_bound) & (recv_y >= y_min_bound - taper_width)
        fy[mask_b] = 0.5 * (
            1.0 + np.cos(np.pi * (y_min_bound - recv_y[mask_b]) / taper_width)
        )
        fy[recv_y > y_max_bound + taper_width] = 0.0
        fy[recv_y < y_min_bound - taper_width] = 0.0

        taper_factor = fx * fy
        mapped_local *= taper_factor
    else:
        # Force zero outside user-specified bounds (hard cutoff)
        in_box = (
            (recv_x >= x_min_bound)
            & (recv_x <= x_max_bound)
            & (recv_y >= y_min_bound)
            & (recv_y <= y_max_bound)
        )
        if np.any(~in_box):
            mapped_local[~in_box] = 0.0

    # Gatherv results
    if rank == 0:
        mapped = np.empty(N, dtype=float)
    else:
        mapped = None

    comm.Gatherv(mapped_local, [mapped, counts, displs, MPI.DOUBLE], root=0)

    # Report diagnostics on OOB ratio
    local_oob = int(np.count_nonzero(~ib))
    total_oob = comm.reduce(local_oob, op=MPI.SUM, root=0)
    # Exclusion diagnostics
    total_excl = comm.reduce(excl_zeroed, op=MPI.SUM, root=0)
    if rank == 0 and N > 0:
        frac = total_oob / float(N)
        print(f"Out-of-bounds nodes: {total_oob}/{N} ({frac:.2%}), oob_mode={oob_mode}")
        if exclude_rects:
            print(
                f"Excluded (zeroed) nodes inside rectangles: {total_excl}/{N} ({total_excl / float(N):.2%}) from {len(exclude_rects)} rectangle(s)"
            )

    # Rank 0 writes CSV and optional plot
    if rank == 0:
        output_data = np.column_stack((x_coord, y_coord, mapped))
        np.savetxt(out_csv, output_data, delimiter=",", comments="")
        print(f"Mapped data saved to {out_csv}")

        if plot:
            try:
                import matplotlib

                matplotlib.use("Agg")  # headless-safe
                import matplotlib.pyplot as plt

                # Plot 1: Grid contour of generated field
                fig1 = plt.figure(figsize=(10, 6))
                ax1 = fig1.add_subplot(111)
                c1 = ax1.contourf(
                    x, y, mean_val + std_val * field, levels=50, cmap="viridis"
                )
                fig1.colorbar(c1, ax=ax1, label="Damage Variable")
                ax1.set_title("Von Kármán Damage Variable (grid)")
                ax1.set_xlabel("X (m)")
                ax1.set_ylabel("Y (m)")
                fig1.tight_layout()
                grid_path = f"{plot_prefix}_grid.png"
                fig1.savefig(grid_path, dpi=200)
                plt.close(fig1)
                print(f"Saved grid plot to {grid_path}")

                # Plot 2: Scatter of mapped Exodus node values
                fig2 = plt.figure(figsize=(10, 6))
                ax2 = fig2.add_subplot(111)
                mask = ~np.isnan(mapped)
                sc = ax2.scatter(
                    x_coord[mask], y_coord[mask], c=mapped[mask], s=5, cmap="viridis"
                )
                fig2.colorbar(sc, ax=ax2, label="Damage Variable (mapped)")
                ax2.set_title("Mapped Damage Variable at Exodus Nodes")
                ax2.set_xlabel("X (m)")
                ax2.set_ylabel("Y (m)")
                fig2.tight_layout()
                pts_path = f"{plot_prefix}_points.png"
                fig2.savefig(pts_path, dpi=200)
                plt.close(fig2)
                print(f"Saved node scatter plot to {pts_path}")
            except Exception as e:
                warnings.warn(f"Plotting failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="MPI von Kármán field generator and mapper"
    )
    parser.add_argument("--nx", type=int, default=640)
    parser.add_argument("--ny", type=int, default=80)
    parser.add_argument("--x-min", type=float, default=-16000)
    parser.add_argument("--x-max", type=float, default=16000)
    parser.add_argument("--y-min", type=float, default=-2000)
    parser.add_argument("--y-max", type=float, default=2000)
    parser.add_argument(
        "--corr-len", type=float, default=200, help="Correlation length (m)"
    )
    parser.add_argument("--nu", type=float, default=0.5, help="Matérn smoothness nu")
    parser.add_argument("--mean", type=float, default=0.1, help="Mean damage variable")
    parser.add_argument(
        "--std", type=float, default=0.025, help="Std dev of damage variable"
    )
    parser.add_argument(
        "--exo", type=str, default="./static_solve_out.e", help="Exodus/NetCDF file"
    )
    parser.add_argument(
        "--out",
        type=str,
        default="./mapped_vonkarman_field.csv",
        help="Output CSV path",
    )
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument(
        "--plot", action="store_true", help="Enable plotting and save PNGs on rank 0"
    )
    parser.add_argument(
        "--plot-prefix",
        type=str,
        default="mapped_vonkarman",
        help="Prefix for saved plot files",
    )
    parser.add_argument(
        "--oob",
        type=str,
        default="nearest",
        choices=["nearest", "nan", "clip"],
        help="Out-of-bounds handling for interpolation",
    )
    parser.add_argument(
        "--auto-extent",
        action="store_true",
        help="Use Exodus coord bounds for grid extents (overrides x/y min/max)",
    )
    parser.add_argument(
        "--taper-width",
        type=float,
        default=0.0,
        help="Cosine taper width outside the box (m)",
    )
    parser.add_argument(
        "--exclude-rect",
        action="append",
        nargs=4,
        type=float,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Add an exclusion rectangle where damage is forced to zero; can repeat.",
    )
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Seed per-rank for independence while reproducible overall
    base_seed = args.seed if args.seed is not None else np.random.SeedSequence().entropy
    ss = np.random.SeedSequence(int(base_seed), spawn_key=(rank,))
    rng = np.random.default_rng(ss)

    # Optionally compute grid extents from Exodus coords on rank 0
    if args.auto_extent and rank == 0:
        try:
            nc = netCDF4.Dataset(args.exo)
            ex = np.array(nc.variables["coordx"][:])
            ey = np.array(nc.variables["coordy"][:])
            nc.close()
            x_min, x_max = float(np.nanmin(ex)), float(np.nanmax(ex))
            y_min, y_max = float(np.nanmin(ey)), float(np.nanmax(ey))
        except Exception as e:
            warnings.warn(f"auto-extent failed: {e}; falling back to CLI extents")
            x_min, x_max = args.x_min, args.x_max
            y_min, y_max = args.y_min, args.y_max
    elif args.auto_extent:
        x_min = x_max = y_min = y_max = None
    else:
        x_min, x_max = args.x_min, args.x_max
        y_min, y_max = args.y_min, args.y_max

    if args.auto_extent:
        x_min = comm.bcast(x_min, root=0)
        x_max = comm.bcast(x_max, root=0)
        y_min = comm.bcast(y_min, root=0)
        y_max = comm.bcast(y_max, root=0)

    # Generate field on rank 0 then broadcast to all
    if rank == 0:
        x, y, field = generate_vonkarman_field_fft(
            args.nx,
            args.ny,
            x_min if args.auto_extent else args.x_min,
            x_max if args.auto_extent else args.x_max,
            y_min if args.auto_extent else args.y_min,
            y_max if args.auto_extent else args.y_max,
            args.nu,
            args.corr_len,
            rng,
        )
    else:
        x = None
        y = None
        field = None

    # Broadcast grids and field
    x = comm.bcast(x, root=0)
    y = comm.bcast(y, root=0)

    # Broadcast field via NumPy buffer friendly Bcast
    if rank != 0:
        field = np.empty((args.ny, args.nx), dtype=float)
    comm.Bcast(field, root=0)

    # Prepare exclusion rectangles list
    exclude_rects = []
    if args.exclude_rect:
        for r in args.exclude_rect:
            # each r is list of 4 floats
            exclude_rects.append(tuple(r))

    # Interpolate and write in parallel
    parallel_interpolate_and_write(
        comm,
        x,
        y,
        field,
        args.exo,
        args.out,
        plot=args.plot,
        mean_val=args.mean,
        std_val=args.std,
        plot_prefix=args.plot_prefix,
        oob_mode=args.oob,
        x_min_bound=args.x_min,
        x_max_bound=args.x_max,
        y_min_bound=args.y_min,
        y_max_bound=args.y_max,
        taper_width=args.taper_width,
        exclude_rects=exclude_rects,
    )

    if rank == 0:
        print(f"Grid extents used: x=[{x[0]}, {x[-1]}], y=[{y[0]}, {y[-1]}]")


if __name__ == "__main__":
    main()
