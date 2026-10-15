#!/usr/bin/env python3
"""
render_figures.py -- crack snapshots and summary tables of the 16 Frontera restart runs (3PB new samples).

WHAT
    Reads, per case, the tetrahedra that cross the mid-thickness plane z = Z_MM at the last complete output
    (<case>.npz, written on Frontera by scripts/frontera_extract_slices.py), cuts them exactly at z = Z_MM and draws
    the damage d (0..1, Inferno) on that plane, axes in mm, one labelled panel per case, five figures by class.
    Writes the case-catalogue and results tables (Markdown) from the same data plus data/hard_steps.txt (crack
    onsets), data/job_case_map.txt (job id per case, from the Frontera log names) and the restart schedule
    (../../logs/restart_schedule.tsv: probe or fitted crack time, start of the fine dt window).
STEPS
    1 load_cases()     npz, onset, job id, schedule row and geometry (make_new_sample_geos.CASES) per case
    2 plane_cut()      exact z-plane cut of the tetrahedra -> triangles with linearly interpolated d
    3 draw_figure()    one figure per class group, one labelled panel per case
    4 write_tables()   Markdown tables for REPORT.md (data/tables.md)
CHECKS
    every case has its npz, onset, job id and schedule row; every panel's cut has triangles
USAGE
    python render_figures.py NPZ_DIR
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, TOOLS)
import make_new_sample_geos as geos  # noqa: E402

Z_MM = 4.0 + 1e-5               # mid-thickness plane (nudged so no mesh node lies exactly on it)
X_WINDOW_MM = (7.0, 21.0)       # panel window around the notch at x = 14 mm
D_BROKEN = 0.95                 # fully broken
HOLE_Y_MM = geos.NOTCH_TIP_MM[1] + 3.6
INK, INK2 = "#0b0b0b", "#52514e"
GROUPS = [                      # figure file, title, case ids
    ("fig_c_asym.png", "(c) two holes, not connected", ["C1", "C2", "C3"]),
    ("fig_d_conn_DL15.png", "(d) two connected holes: left D = 1.5 mm, right D = 2.0 mm",
     ["D1", "D2", "D3", "D4", "D5"]),
    ("fig_d_conn_DL20.png", "(d) two connected holes: both D = 2.0 mm", ["D6", "D7", "D8", "D9"]),
    ("fig_e_sine.png", "(e) sine-wave void", ["E1", "E2"]),
    ("fig_f_pillars.png", "(f) pillar structures (welded mismatch)", ["F1", "F2"]),
]
ORDER = [                       # report id, case tag (make_new_sample_geos.case_tag)
    ("C1", "xL3.00mm_xR2.00mm_y3.6mm_radi1.0mm"), ("C2", "xL3.00mm_xR2.50mm_y3.6mm_radi1.0mm"),
    ("C3", "xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm"),
    ("D1", "DL1.5mm_DR2.0mm_x1.00mm_y3.6mm_overlap"), ("D2", "DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_downline"),
    ("D3", "DL1.5mm_DR2.0mm_x1.50mm_y3.6mm_upline"), ("D4", "DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_downline"),
    ("D5", "DL1.5mm_DR2.0mm_x2.00mm_y3.6mm_midline"), ("D6", "DL2.0mm_DR2.0mm_x1.00mm_y3.6mm_overlap"),
    ("D7", "DL2.0mm_DR2.0mm_x1.50mm_y3.6mm_downline"), ("D8", "DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_downline"),
    ("D9", "DL2.0mm_DR2.0mm_x2.00mm_y3.6mm_upline"),
    ("E1", "D1.0mm_y1.6mm_amp2.0mm_centered"), ("E2", "D1.0mm_y1.6mm_amp2.0mm_1ampoffset"),
    ("F1", "welded_mismatch_2onN"), ("F2", "welded_mismatch_3onN"),
]
DISCOVERY_ONSET = {"C1": 21.614, "D2": 27.749, "E1": 22.720, "F1": 24.242}   # first hard step, Discovery runs


def short_name(case):
    """One-line label of a case (panel titles, tables)."""
    f = case["family"]
    if f == "2holes_asym":
        return "(c) holes xL %.2f / xR %.2f mm" % (case["dx_left"], case["dx_right"])
    if f == "2holes_connected":
        return "(d) DL%.1f/DR%.1f, x %.2f mm, %s" % (2 * case["r_left"], 2 * case["r_right"], case["dx"], case["line"])
    if f == "sinewave":
        return "(e) sine void, %s" % case["variant"]
    return "(f) pillars, %s" % case["variant"]


def geometry_words(case):
    """At most 15 words describing the case geometry (catalogue table)."""
    f = case["family"]
    if f == "2holes_asym":
        return "holes r 1.0 at x 14-%.2f and 14+%.2f, y %.1f" % (case["dx_left"], case["dx_right"], HOLE_Y_MM)
    if f == "2holes_connected":
        slot = geos.connected_slot_mm(case)
        where = "slot y %.2f-%.2f" % slot if slot else "holes touch, no slot"
        return "holes r %.2f/%.2f at x 14 -/+ %.2f, y %.1f; %s" % (case["r_left"], case["r_right"], case["dx"],
                                                                     HOLE_Y_MM, where)
    if f == "sinewave":
        return "sine tube r 0.5 at y 3.2, amplitude 2, axis x %.1f" % case["x_axis"]
    return "blocks gap y 4.0-4.5, pillars r 0.5, %s rows" % case["variant"]


def tip_note(case, tip):
    """Where the crack tip sits relative to the case's own features (results table)."""
    f = case["family"]
    if f == "2holes_asym":
        top = HOLE_Y_MM + case["radius"]
        return "between the holes" if tip < top else "above the hole tops (y %.1f)" % top
    if f == "2holes_connected":
        slot = geos.connected_slot_mm(case)
        if slot is None:
            return "above the touching holes" if tip > HOLE_Y_MM else "below the touching holes"
        if abs(tip - slot[0]) <= 0.05:
            return "held at the slot bottom"
        return "below the slot" if tip < slot[0] else "past the slot"
    if f == "sinewave":
        return "past the void band (y 2.7-3.7)" if tip > 3.7 else "in or below the void band"
    return "past the pillar gap (y 4.0-4.5)" if tip > 4.5 else "below the pillar gap"


def tag_of_log_name(name):
    """Case tag from a Frontera log name '<family>_<tag>' (data/job_case_map.txt)."""
    for family in ("2holes_asym", "2holes_connected", "sinewave", "pillars"):
        if name.startswith(family + "_"):
            return name[len(family) + 1:]
    raise ValueError("unknown family in log name %r" % name)


def load_cases(npz_dir):
    """Step 1: everything the figures and tables need, per case id."""
    onsets = {}
    for line in open(os.path.join(HERE, "data", "hard_steps.txt")):
        m = re.match(r"H (\d+) .*first=([0-9.]+)", line)
        if m:
            onsets[int(m.group(1))] = float(m.group(2))
    jobs = {}
    for line in open(os.path.join(HERE, "data", "job_case_map.txt")):
        job, name = line.split()
        jobs[tag_of_log_name(name)] = int(job)
    schedule = {}
    for line in open(os.path.join(TOOLS, "logs", "restart_schedule.tsv")):
        p = line.rstrip("\n").split("\t")
        if p[0] != "case":
            schedule[p[0]] = dict(crack=float(p[2]), source=p[3], fine=float(p[4]))
    by_tag = {geos.case_tag(c): c for c in geos.CASES}
    cases = {}
    for cid, tag in ORDER:
        case = by_tag[tag]
        folder = "case_%s_3d_%s" % (case["family"], tag)
        npz = np.load(os.path.join(npz_dir, folder + ".npz"))
        job = jobs[tag]
        cases[cid] = dict(case=case, tag=tag, job=job, onset=onsets[job], sched=schedule[tag],
                          xyz=npz["xyz"], d=npz["d"], tets=npz["tets"], t=float(npz["t"]), tip=float(npz["tip"]),
                          top=float(npz["top"]))
    return cases


def plane_cut(xyz, d, tets, z0):
    """Step 2: exact cut of tetrahedra by the plane z = z0 -> (k, 3) triangle x, y and d (triangles and split quads)."""
    s = xyz[tets, 2] - z0
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    ia, ib = [a for a, _ in edges], [b for _, b in edges]
    sa, sb = s[:, ia], s[:, ib]
    cross = (sa < 0) != (sb < 0)
    lam = np.where(cross, sa / np.where(cross, sa - sb, 1.0), np.nan)
    pa, pb = xyz[tets[:, ia]], xyz[tets[:, ib]]
    px = pa[..., 0] + lam * (pb[..., 0] - pa[..., 0])
    py = pa[..., 1] + lam * (pb[..., 1] - pa[..., 1])
    pv = d[tets[:, ia]] + lam * (d[tets[:, ib]] - d[tets[:, ia]])
    count = cross.sum(axis=1)
    order = np.argsort(~cross, axis=1, kind="stable")       # crossing edges first
    take = lambda arr, rows, k: np.take_along_axis(arr[rows], order[rows][:, :k], axis=1)
    tri3 = count == 3
    xs, ys, vs = [take(px, tri3, 3)], [take(py, tri3, 3)], [take(pv, tri3, 3)]
    quad = count == 4
    qx, qy, qv = take(px, quad, 4), take(py, quad, 4), take(pv, quad, 4)
    ang = np.arctan2(qy - qy.mean(axis=1, keepdims=True), qx - qx.mean(axis=1, keepdims=True))
    o = np.argsort(ang, axis=1)
    qx, qy, qv = (np.take_along_axis(a, o, axis=1) for a in (qx, qy, qv))
    for i, j, k in ((0, 1, 2), (0, 2, 3)):
        xs.append(qx[:, [i, j, k]])
        ys.append(qy[:, [i, j, k]])
        vs.append(qv[:, [i, j, k]])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(vs)


def draw_panel(ax, item, cid):
    """One case: the cut coloured by d, the mid-thickness tip marked, a two-line label."""
    tx, ty, tv = plane_cut(item["xyz"], item["d"], item["tets"], Z_MM)
    assert len(tx), "%s: empty cut" % cid
    tri = np.arange(3 * len(tx)).reshape(-1, 3)
    im = ax.tripcolor(tx.ravel(), ty.ravel(), tri, tv.ravel(), shading="gouraud", cmap="inferno",
                      vmin=0.0, vmax=1.0, rasterized=True)
    broken = tv.ravel() > D_BROKEN
    if broken.any():
        k = np.argmax(np.where(broken, ty.ravel(), -1))
        ax.plot(tx.ravel()[k], ty.ravel()[k], marker="v", ms=5, mfc="#00e5ff", mec="black", mew=0.5)
    ax.set_xlim(*X_WINDOW_MM)
    ax.set_ylim(0, item["top"])
    ax.set_aspect("equal")
    ax.set_title("%s   %s\njob %d   |   crack onset t = %.3f\nsnapshot t = %.1f   |   crack tip y = %.2f mm"
                 % (cid, short_name(item["case"]), item["job"], item["onset"], item["t"], item["tip"]),
                 fontsize=7, color=INK, loc="left", linespacing=1.25)
    ax.tick_params(labelsize=6, colors=INK2)
    ax.set_xlabel("x (mm)", fontsize=6.5, color=INK2, labelpad=1)
    ax.set_ylabel("y (mm)", fontsize=6.5, color=INK2, labelpad=1)
    return im


def draw_figure(cases, fname, title, ids, out_dir):
    """Step 3: one figure for a class group, two panels per row, one shared colour bar."""
    rows = (len(ids) + 1) // 2
    fig = plt.figure(figsize=(6.5, 0.4 + 2.5 * rows), dpi=200)
    gs = fig.add_gridspec(rows, 3, width_ratios=[1, 1, 0.04], hspace=0.75, wspace=0.18)
    im = None
    for k, cid in enumerate(ids):
        ax = fig.add_subplot(gs[k // 2, k % 2])
        im = draw_panel(ax, cases[cid], cid)
    cax = fig.add_subplot(gs[:, 2])
    cb = fig.colorbar(im, cax=cax, ticks=[0, 0.25, 0.5, 0.75, 1.0])
    cb.set_label("damage d (slice z = 4 mm)", fontsize=7, color=INK2)
    cb.ax.tick_params(labelsize=6, colors=INK2)
    fig.suptitle(title, fontsize=8.5, color=INK, x=0.06, ha="left")
    path = os.path.join(out_dir, fname)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def write_tables(cases, path):
    """Step 4: case catalogue and results tables (Markdown, cells of at most 15 words)."""
    lines = ["| ID | Class | Case | Geometry (mm) | STL file | Frontera job |", "|---|---|---|---|---|---|"]
    for cid, _ in ORDER:
        it = cases[cid]
        lines.append("| %s | (%s) | `%s` | %s | `%s` | %d |" % (cid, cid[0].lower(), it["tag"], geometry_words(it["case"]),
                                                           os.path.basename(it["case"]["stl"]), it["job"]))
    lines += ["", "| ID | Onset t | Estimate t (source) | Fine dt from | Frame t | Tip y (mm) | To top (mm) | Tip position |",
              "|---|---|---|---|---|---|---|---|"]
    for cid, _ in ORDER:
        it = cases[cid]
        est = "%.1f (%s)" % (it["sched"]["crack"], "probe" if it["sched"]["source"] == "measured" else "fit")
        if cid in DISCOVERY_ONSET:
            est += "; Discovery %.3f" % DISCOVERY_ONSET[cid]
        lines.append("| %s | %.3f | %s | %.1f | %.1f | %.2f | %.2f | %s |"
                     % (cid, it["onset"], est, it["sched"]["fine"], it["t"], it["tip"], it["top"] - it["tip"],
                        tip_note(it["case"], it["tip"])))
    open(path, "w").write("\n".join(lines) + "\n")
    print("wrote", path)


def main(npz_dir):
    cases = load_cases(npz_dir)
    out_dir = os.path.join(HERE, "figures")
    os.makedirs(out_dir, exist_ok=True)
    for fname, title, ids in GROUPS:
        draw_figure(cases, fname, title, ids, out_dir)
    write_tables(cases, os.path.join(HERE, "data", "tables.md"))


if __name__ == "__main__":
    main(sys.argv[1])
