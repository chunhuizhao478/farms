# Meshes for the SPE-HFTC acrylic-cube cases

Geometry from the manuscript (Table 2) and the decisions of 2026-10-21:
cube 51.11 x 51.92 x 49.14 mm, borehole diameter 8.3 mm along z through the
centre, box spanning x in [-Lx/2, Lx/2], y in [-Ly/2, Ly/2], z in [0, Lz].
`*.msh` files are git-ignored (regenerate with the commands below; gmsh 4.15
python API, `python3 -c "import gmsh"` must work). Each run also writes a
`<name>.json` with counts, areas, volume and quality.

## 3D cube (`make_spe_cube_mesh.py`)

Named physical groups (MOOSE uses the names): volume `acrylic`; surfaces
`borehole`, `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, `zmax`. Faces are
classified by their centre of mass, never by OCC tag. The fine zone is the
cylinder r <= r_fine over the full height, with a tanh transition to h_coarse.

| name | command | tets | nodes | min gamma |
| --- | --- | --- | --- | --- |
| `spe_cube_smoke` | `python3 make_spe_cube_mesh.py --h-fine 1.5e-3 --h-coarse 1.5e-3 --out spe_cube_smoke` | 0.15 M | 28.6 k | 0.22 |
| `spe_cube_h0p3` (production) | `python3 make_spe_cube_mesh.py --h-fine 0.3e-3 --h-coarse 1.0e-3 --r-fine 28e-3 --out spe_cube_h0p3 --threads 12` | 14.6 M | 2.4 M | 0.23 |

Production: 0.3 mm elements wherever r <= ~26.5 mm (covers the face centres
and the 20.45 mm extent of the paper), 1 mm in the four corners; built with
HXT in about one minute on 12 threads; 712 MB on disk (14 586 319 tets, 2 381 878 nodes). Checks in the json:
borehole area 1.2813e-3 m^2 and volume 1.2774e-4 m^3 match the analytic values.

## 2D plane-strain section (`make_spe_section_2d.py`)

The section of the cube at the electrode-gap height: rectangle minus the
borehole disk, uniform element size. Physical groups: surface `acrylic`;
curves `borehole`, `xmin`, `xmax`, `ymin`, `ymax`.

| name | command | triangles | nodes | phase-field l (2 elements) |
| --- | --- | --- | --- | --- |
| `spe_section_h0p3` | `python3 make_spe_section_2d.py --h 0.3e-3 --out spe_section_h0p3` | 67 k | 34 k | 0.6 mm |
| `spe_section_h0p15` | `python3 make_spe_section_2d.py --h 0.15e-3 --out spe_section_h0p15` | 0.27 M | 0.13 M | 0.3 mm |
| `spe_section_h0p075` | `python3 make_spe_section_2d.py --h 0.075e-3 --out spe_section_h0p075` | 1.07 M | 0.54 M | 0.15 mm |
| `spe_section_h0p0385` (not built) | `python3 make_spe_section_2d.py --h 0.0385e-3 --out spe_section_h0p0385` | ~4 M | ~2 M | 77 um (true acrylic) |

With the strength-derived toughness Gc = 8 l ft^2/(3E), the input files
recompute Gc from `l`; pass `l=` and `mesh_file=` together on the command line,
e.g. `l=1.5e-4 mesh_file=../mesh/spe_section_h0p075.msh`.

## Physical ids written by gmsh (for reference; the inputs use the names)

3D: xmin 1, ymin 2, zmax 3, ymax 4, zmin 5, xmax 6, borehole 7, acrylic 8.
2D: borehole 1, ymin 2, xmin 3, xmax 4, ymax 5, acrylic 6.
