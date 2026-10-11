#!/usr/bin/env python3
"""
make_new_sample_geos.py -- write one gmsh .geo per new three-point-bending sample.

WHAT
    Writes meshfile_<family>/mesh_<family>_3d_<tag>.geo for the 16 specimens in
    Paper_ThreePointBending/new_samples/New sample stls.  Every file follows
    meshfile_2holes/: OpenCASCADE kernel, metres, the standard 28 x 8 x 8 mm
    notched beam outline with the support/loading points 25/26/27 (their extruded
    edges carry the MOOSE line BCs), Box size fields, one Physical Volume "volume"
    (MOOSE block 1), MSH 2.2 output.  Family-specific geometry:

      2holes_asym       (c) two r = 1 mm holes 3.6 mm above the notch tip,
                            left/right offsets as in the STL file name
      2holes_connected  (d) a r = 0.75/1.0 or 1.0/1.0 mm hole pair at -+dx joined
                            by a 0.2 mm slot (up/mid/down line) or not (overlap);
                            the touching 2.0/2.0 pair is joined by a throat that
                            trims the two zero-width material cusps
      sinewave          (e) r = 0.5 mm void 1.6 mm above the notch tip whose
                            centre follows x0 + 2 sin(2 pi z / 8 mm)
      pillars           (f) two 4 mm blocks joined across a 0.5 mm gap by
                            r = 0.5 mm pillars on a staggered 1.5 mm grid

    All numbers in CASES were measured from the STLs (circle fits to the STL
    vertices, residual < 1e-6 mm; sine centreline within 0.03 mm);
    verify_new_sample_meshes.py checks every mesh against its STL.

    Notch: the STLs draw a flat-tipped 1.0 x 1.6 mm slot.  Every earlier
    meshfile_* set, the paper's setup figure and the X-ray frames use a
    semicircular tip of radius 0.5 mm, so NOTCH_ROUND = 1 keeps that tip;
    NOTCH_ROUND = 0 reproduces the STL notch exactly (one line in each .geo).

STEPS
    1. notch_outline_geo()  -- outline points/lines with notch and support points
    2. <family>_body_geo()  -- beam volume(s), family voids/pillars, Physical Volume;
                               family (d) cuts its void as an explicit inner outline
                               of the section instead (connected_void_outline_geo)
    3. size_fields_geo()    -- Box fields: notch band, feature boxes, supports
    4. write_case_geo()     -- assemble and write one file per case

CHECKS
    check_geo_volume() opens each written .geo with the gmsh API (no meshing),
    requires exactly one volume and compares its OCC volume with the analytic
    volume expected_volume_mm3() of the case.

USAGE
    python make_new_sample_geos.py [--check]
"""

import argparse
import math
import os
import sys

# ----------------------------------------------------------------------------- settings
TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
TPB_DIR = os.path.dirname(TOOLS_DIR)                  # .../threepointbending
STL_ROOT = "/Users/chunhuizhao/Documents/Paper_ThreePointBending/new_samples/New sample stls"

LC_GLOBAL = 1e-3          # m, global size (all earlier sets)
LC_REFINED = 1e-4         # m, notch/feature boxes (as meshfile_2holes; 2 l with l = 50 um)
LC_SUPPORT = 2.5e-4       # m, support boxes (all earlier sets)
LC_PILLAR_LAYER = 1.25e-4  # m, whole pillar layer of family (f)
NOTCH_ROUND = 1           # 1 = semicircular notch tip (earlier sets), 0 = flat tip (STL)

BEAM_LENGTH_MM = 28.0
BEAM_HEIGHT_MM = 8.0
BEAM_THICKNESS_MM = 8.0
NOTCH_TIP_MM = (14.0, 1.6)
NOTCH_WIDTH_MM = 1.0
CONNECTING_LINE_WIDTH_MM = 0.2   # (d) slot width measured from the STLs
CUSP_MIN_WIDTH_MM = 0.1          # (d) touching holes: trim cusps narrower than this (= LC_REFINED)

FAMILY_TITLES = {
    "2holes_asym": "(c) 2 holes not connected (same size, asymmetric offset)",
    "2holes_connected": "(d) 2 holes connected by thin line (diverging crack path)",
    "sinewave": "(e) sine wave hole",
    "pillars": "(f) pillar structures (welded mismatch)",
}

# (c) offsets: "A-B mm apart" = right hole A mm, left hole B mm from the notch tip
# (d) holes: left r_left at x = 14 - dx, right r_right at x = 14 + dx; line = slot position
# (e) x_axis = sine centreline at the faces z = 0 and z = 8 mm
# (f) row_x0/row_n: first pillar x and pillar count of the rows at z = 1, 2.5, 4, 5.5, 7 mm
CASES = [
    dict(family="2holes_asym", stl="(c) 2 holes not connected (Same size, asymmetric offset)/2mm_circle_3.6up_2.0-3.0mm apart v2.stl",
         experiment="1-Alum_2mmD_3.6up_2.0-3.0mm_apt", dx_left=3.0, dx_right=2.0, dy=3.6, radius=1.0),
    dict(family="2holes_asym", stl="(c) 2 holes not connected (Same size, asymmetric offset)/2mm_circle_3.6up_2.5-3.0mm apart v2.stl",
         experiment="2-Alum_2mmD_3.6up_2.5-3.0mm_apt", dx_left=3.0, dx_right=2.5, dy=3.6, radius=1.0),
    dict(family="2holes_asym", stl="(c) 2 holes not connected (Same size, asymmetric offset)/2mm_circle_3.6up_2.5-3.5mm apart v2.stl",
         experiment="3-Alum_2mmD_3.6up_2.5-3.5mm_apt", dx_left=3.5, dx_right=2.5, dy=3.6, radius=1.0),

    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 1/Alum_2.0_1.5mmD_y3.6mm_x1.5mm_upline.stl",
         experiment="1.1 Alum_2.0_1.5mmD_y3.6mm_x1.5mm_upline", r_left=0.75, r_right=1.0, dx=1.5, dy=3.6, line="upline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 1/Alum_2.0_1.5mmD_y3.6mm_x1.5mm_downline.stl",
         experiment="1.2 Alum_2.0_1.5mmD_y3.6mm_x1.5mm_downline", r_left=0.75, r_right=1.0, dx=1.5, dy=3.6, line="downline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 1/Alum_2.0_1.5mmD_y3.6mm_x1.0mm_overlap.stl",
         experiment="1.3 Alum_2.0_1.5mmD_y3.6mm_x1.0mm_overlap", r_left=0.75, r_right=1.0, dx=1.0, dy=3.6, line="overlap"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 2/Alum_2.0_1.5mmD_y3.6mm_x2.0mm_midline.stl",
         experiment="2.1 Alum_2.0_1.5mmD_y3.6mm_x2.0mm_midline", r_left=0.75, r_right=1.0, dx=2.0, dy=3.6, line="midline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 2/Alum_2.0_1.5mmD_y3.6mm_x2.0mm_downline.stl",
         experiment="2.2 Alum_2.0_1.5mmD_y3.6mm_x2.0mm_downline", r_left=0.75, r_right=1.0, dx=2.0, dy=3.6, line="downline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 3/Alum_2.0_2.0mmD_y3.6mm_x2.0mm_upline.stl",
         experiment="3.1 Alum_2.0_2.0mmD_y3.6mm_x2.0mm_upline", r_left=1.0, r_right=1.0, dx=2.0, dy=3.6, line="upline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 3/Alum_2.0_2.0mmD_y3.6mm_x2.0mm_downline.stl",
         experiment="3.2 Alum_2.0_2.0mmD_y3.6mm_x2.0mm_downline", r_left=1.0, r_right=1.0, dx=2.0, dy=3.6, line="downline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 3/Alum_2.0_2.0mmD_y3.6mm_x1.5mm_downline.stl",
         experiment="3.3 Alum_2.0_2.0mmD_y3.6mm_x1.5mm_downline", r_left=1.0, r_right=1.0, dx=1.5, dy=3.6, line="downline"),
    dict(family="2holes_connected", stl="(d) 2 holes connected by thin line (diverging crack path)/connected 3/Alum_2.0_2.0mmD_y3.6mm_x1.0mm_overlap.stl",
         experiment="3.4 Alum_2.0_2.0mmD_y3.6mm_x1.0mm_overlap", r_left=1.0, r_right=1.0, dx=1.0, dy=3.6, line="overlap"),

    dict(family="sinewave", stl="(e) Sine wave hole/Alum_1mmD_1.6up_sinewavevoid_Amp2_centered.stl",
         experiment="Centered sine wave", x_axis=14.0, amp=2.0, wavelength=8.0, radius=0.5, dy=1.6, variant="centered"),
    dict(family="sinewave", stl="(e) Sine wave hole/Alum_1mmD_1.6up_sinewavevoid_Amp2_1ampoffset.stl",
         experiment="Offset sine wave", x_axis=16.0, amp=2.0, wavelength=8.0, radius=0.5, dy=1.6, variant="1ampoffset"),

    dict(family="pillars", stl="(f)  Pillar structures/Welded_mismatch_2onN v2.stl",
         experiment="1-Alum_welded_mismatch_2onN_1", variant="2onN",
         row_z=[1.0, 2.5, 4.0, 5.5, 7.0], row_x0=[1.25, 2.0, 1.25, 2.0, 1.25], row_n=[18, 17, 18, 17, 18]),
    dict(family="pillars", stl="(f)  Pillar structures/Welded_mismatch_3onN v2.stl",
         experiment="2-Alum_welded_mismatch_3onN_1", variant="3onN",
         row_z=[1.0, 2.5, 4.0, 5.5, 7.0], row_x0=[2.0, 1.25, 2.0, 1.25, 2.0], row_n=[17, 18, 17, 18, 17]),
]
PILLAR_GAP_MM = (4.0, 4.5)     # (f) gap between the blocks
PILLAR_TOP_MM = 8.5            # (f) top of the upper block = loading line height
PILLAR_RADIUS_MM = 0.5
PILLAR_PITCH_MM = 1.5


# ----------------------------------------------------------------------------- naming
def case_tag(case):
    """Short case name used for mesh, deck and job file names."""
    f = case["family"]
    if f == "2holes_asym":
        return (f"xL{case['dx_left']:.2f}mm_xR{case['dx_right']:.2f}mm_"
                f"y{case['dy']:.1f}mm_radi{case['radius']:.1f}mm")
    if f == "2holes_connected":
        return (f"DL{2 * case['r_left']:.1f}mm_DR{2 * case['r_right']:.1f}mm_"
                f"x{case['dx']:.2f}mm_y{case['dy']:.1f}mm_{case['line']}")
    if f == "sinewave":
        return (f"D{2 * case['radius']:.1f}mm_y{case['dy']:.1f}mm_"
                f"amp{case['amp']:.1f}mm_{case['variant']}")
    if f == "pillars":
        return f"welded_mismatch_{case['variant']}"
    raise ValueError(f)


def mesh_basename(case):
    return f"mesh_{case['family']}_3d_{case_tag(case)}"


def geo_path(case):
    return os.path.join(TPB_DIR, f"meshfile_{case['family']}", mesh_basename(case) + ".geo")


def beam_top_mm(case):
    """Height of the loading line (top face) in mm."""
    return PILLAR_TOP_MM if case["family"] == "pillars" else BEAM_HEIGHT_MM


def mm(value):
    """Format a length given in mm as metres for the .geo file."""
    return f"{value * 1e-3:.10g}"


# ----------------------------------------------------------------------------- geometry helpers
def connected_slot_mm(case):
    """(d) slot (y0, y1) in mm, or None.  The up/down line is tangent to the smaller hole."""
    hole_y = NOTCH_TIP_MM[1] + case["dy"]
    r_min = min(case["r_left"], case["r_right"])
    w = CONNECTING_LINE_WIDTH_MM
    if case["line"] == "upline":
        return hole_y + r_min - w, hole_y + r_min
    if case["line"] == "downline":
        return hole_y - r_min, hole_y - r_min + w
    if case["line"] == "midline":
        return hole_y - w / 2, hole_y + w / 2
    return None


def holes_touch(case):
    """(d) True when the two holes touch (centre distance = sum of radii)."""
    return abs(2 * case["dx"] - case["r_left"] - case["r_right"]) < 1e-9


def throat_half_height_mm(case):
    """(d) touching equal holes: half height where the material cusp is CUSP_MIN_WIDTH_MM wide."""
    r = case["r_left"]
    assert abs(case["r_right"] - r) < 1e-12, "throat formula assumes equal radii"
    return math.sqrt(r ** 2 - (r - CUSP_MIN_WIDTH_MM / 2) ** 2)


# ----------------------------------------------------------------------------- .geo pieces
def header_geo(case):
    rel_stl = case["stl"]
    return f"""SetFactory("OpenCASCADE");

//==============================================================================
// {FAMILY_TITLES[case['family']]}
// Case       : {case_tag(case)}
// Experiment : {case['experiment']}
// Source STL : New sample stls/{rel_stl}
// Written by threepointbending/new_samples_tools/make_new_sample_geos.py;
// parameters measured from the STL (mm in the STL, metres here).
//==============================================================================

//==============================================================================
// MESH PARAMETERS
//==============================================================================
// Global mesh size parameters
lc = {LC_GLOBAL:g};          // Global mesh size (1 mm)
lc_refined = {LC_REFINED:g};   // Refined mesh size in notch/feature region (100 microns, as meshfile_2holes)
lc_support = {LC_SUPPORT:g}; // Refined mesh size at support regions (250 microns)
extrude_z = {mm(BEAM_THICKNESS_MM)};   // Total thickness in z-direction (8 mm)
Mesh.MshFileVersion = 2.2;   // MSH 2.2 ASCII, as meshfile_height / meshfile_2holes

//==============================================================================
// NOTCH (tip at x = 14 mm, y = 1.6 mm; width 1 mm)
//==============================================================================
notch_tip_x = {mm(NOTCH_TIP_MM[0])};
notch_tip_y = {mm(NOTCH_TIP_MM[1])};
// 1 = semicircular tip of radius 0.5 mm (all earlier meshfile_* sets, paper figure)
// 0 = flat tip, i.e. the rectangular 1.0 x 1.6 mm notch exactly as drawn in the STL
notch_round = {NOTCH_ROUND};
"""


def notch_outline_geo(top_mm, with_load_point, inner_geo="", inner_loops=()):
    """Outline of the notched (lower) beam at z = 0, extruded through the thickness.

    Keeps the point/line numbering of meshfile_2holes: points 25/26 (and 27) sit on the
    bottom (top) edge so the extrusion creates the support (loading) lines along z.
    inner_geo / inner_loops: curve loops cut out of the section (family (d) voids).
    """
    surface_loops = ", ".join(str(t) for t in (1, *inner_loops))
    top = mm(top_mm)
    load_pt = (f"Point(27) = {{ 0.014, {top}, 0, lc_support}};  // Loading line position\n"
               if with_load_point else "")
    top_lines = ("Line(10) = {8, 27};      // Top right to loading line\n"
                 "Line(11) = {27, 9};      // Loading line to top left corner\n"
                 if with_load_point else
                 "Line(10) = {8, 9};       // Top edge\n")
    loop = "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12" if with_load_point else "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12"
    return f"""
//==============================================================================
// BEAM OUTLINE (XY plane at z = 0), extruded in z below
//==============================================================================
Point(1) = {{      0,      0,   0, lc}};
Point(2) = {{ 0.0135,      0,   0, lc}};
If (notch_round == 1)
  Point(3) = {{ 0.0135, 0.0011,   0, lc}};
  Point(4) = {{ 0.0140, 0.0016,   0, lc}};
  Point(5) = {{ 0.0145, 0.0011,   0, lc}};
  Point(14) = {{ 0.0140, 0.0011, 0, lc}};  // Circle center for notch
Else
  Point(3) = {{ 0.0135, 0.0016,   0, lc}};
  Point(4) = {{ 0.0140, 0.0016,   0, lc}};
  Point(5) = {{ 0.0145, 0.0016,   0, lc}};
EndIf
Point(6) = {{ 0.0145,      0,   0, lc}};
Point(7) = {{ 0.0280,      0,   0, lc}};
Point(8) = {{ 0.0280, {top},   0, lc}};
Point(9) = {{      0, {top},   0, lc}};

// POINTS FOR LINE SUPPORTS (at center of original patches)
Point(25) = {{ 0.004,      0, 0, lc_support}};  // Left support line position
Point(26) = {{ 0.024,      0, 0, lc_support}};  // Right support line position
{load_pt}
// Bottom edge: left corner -> left support -> notch -> right support -> right corner
Line(1) = {{1, 25}};       // Bottom left corner to left support line
Line(2) = {{25, 2}};       // Left support line to notch start
Line(3) = {{2, 3}};        // Notch vertical left
If (notch_round == 1)
  Circle(4) = {{3, 14, 4}};  // Notch circle left
  Circle(5) = {{4, 14, 5}};  // Notch circle right
Else
  Line(4) = {{3, 4}};        // Flat notch tip, left half
  Line(5) = {{4, 5}};        // Flat notch tip, right half
EndIf
Line(6) = {{5, 6}};        // Notch vertical right
Line(7) = {{6, 26}};       // Notch end to right support line
Line(8) = {{26, 7}};       // Right support line to right corner
Line(9) = {{7, 8}};        // Right edge
{top_lines}Line(12) = {{9, 1}};       // Left edge
{inner_geo}
Curve Loop(1) = {{{loop}}};
Plane Surface(1) = {{{surface_loops}}};
"""


def extrude_lower_geo():
    """Extrude the notched outline (surface 1) through the thickness."""
    return """
// Extrude in z (tetrahedral elements; size controlled by the background field)
out[] = Extrude {0, 0, extrude_z} {
  Surface{1};
};
main_vol = out[1];
"""


def physical_volume_geo():
    return """
// Physical Volume = MOOSE block 1 (the decks' damage block); every volume of the result
all_vols[] = Volume "*";
Printf("Total volumes found: %g", #all_vols[]);
Physical Volume("volume") = {all_vols[]};
"""


def two_holes_asym_body_geo(case):
    p = case
    return f"""
//==============================================================================
// TWO CYLINDRICAL HOLES (same size, asymmetric offsets from the notch tip)
//==============================================================================
dx_left = {mm(p['dx_left'])};     // Left hole: horizontal distance from notch tip ({p['dx_left']:.2f} mm)
dx_right = {mm(p['dx_right'])};    // Right hole: horizontal distance from notch tip ({p['dx_right']:.2f} mm)
dy_left = {mm(p['dy'])};    // Left hole: vertical distance from notch tip ({p['dy']:.1f} mm)
dy_right = {mm(p['dy'])};   // Right hole: vertical distance from notch tip ({p['dy']:.1f} mm)
hole_radius_left = {mm(p['radius'])};    // Left hole radius ({p['radius']:.1f} mm)
hole_radius_right = {mm(p['radius'])};   // Right hole radius ({p['radius']:.1f} mm)

hole_left_x = notch_tip_x - dx_left;
hole_left_y = notch_tip_y + dy_left;
hole_right_x = notch_tip_x + dx_right;
hole_right_y = notch_tip_y + dy_right;

Cylinder(100) = {{hole_left_x, hole_left_y, 0, 0, 0, extrude_z, hole_radius_left}};
Cylinder(101) = {{hole_right_x, hole_right_y, 0, 0, 0, extrude_z, hole_radius_right}};
BooleanDifference{{ Volume{{main_vol}}; Delete; }}{{ Volume{{100}}; Volume{{101}}; Delete; }}
Printf("Left hole at (%g, %g), r = %g; right hole at (%g, %g), r = %g",
       hole_left_x, hole_left_y, hole_radius_left, hole_right_x, hole_right_y, hole_radius_right);
"""


def two_holes_connected_params_geo(case):
    """(d) parameter block; the void itself is cut by connected_void_outline_geo()."""
    p = case
    line_code = {"upline": 1, "midline": 2, "downline": 3, "overlap": 0}[p["line"]]
    s = f"""
//==============================================================================
// TWO HOLES (left smaller or equal) CONNECTED BY A THIN LINE (0.2 mm slot)
//==============================================================================
dx_left = {mm(p['dx'])};     // Left hole: horizontal distance from notch tip ({p['dx']:.2f} mm)
dx_right = {mm(p['dx'])};    // Right hole: horizontal distance from notch tip ({p['dx']:.2f} mm)
dy = {mm(p['dy'])};         // Both holes: vertical distance from notch tip ({p['dy']:.1f} mm)
hole_radius_left = {mm(p['r_left'])};   // Left hole radius (D = {2 * p['r_left']:.1f} mm)
hole_radius_right = {mm(p['r_right'])};  // Right hole radius (D = {2 * p['r_right']:.1f} mm)
// Connecting line: 0 = none ("overlap" samples), 1 = up, 2 = mid, 3 = down.
// Up/down lines are tangent to the top/bottom of the smaller hole; mid is centred.
line_type = {line_code};
line_width = {mm(CONNECTING_LINE_WIDTH_MM)};   // 0.2 mm

hole_left_x = notch_tip_x - dx_left;
hole_right_x = notch_tip_x + dx_right;
hole_y = notch_tip_y + dy;
"""
    if holes_touch(p):
        s += f"""
// The two holes TOUCH (centre distance = r_left + r_right): the material between them ends
// in two zero-angle cusps at the contact line, which cannot be meshed.  The X-ray shows one
// merged figure-8 void, so the holes are joined by a throat that trims each cusp where it is
// narrower than cusp_min_width (= lc_refined): for equal radii r the cusp width at height s
// from the contact line is 2 (r - Sqrt(r^2 - s^2)), so the throat half height is throat_h.
cusp_min_width = {mm(CUSP_MIN_WIDTH_MM)};
throat_h = Sqrt(hole_radius_left^2 - (hole_radius_left - cusp_min_width / 2)^2);
"""
    return s


def connected_band_mm(case):
    """(d) y-band (y0, y1) in mm of the slot or throat joining the holes; None if separate."""
    if holes_touch(case):
        yc = NOTCH_TIP_MM[1] + case["dy"]
        h = throat_half_height_mm(case)
        return yc - h, yc + h
    return connected_slot_mm(case)


def arc_break_angles(theta_start, theta_end):
    """Split a counter-clockwise arc (rad, end > start) at the multiples of pi/2 inside it,
    so every piece spans at most pi/2 (OCC three-point arcs must be shorter than pi)."""
    quarter = math.pi / 2
    angles = [theta_start]
    k = math.floor(theta_start / quarter) + 1
    while k * quarter < theta_end - 1e-9:
        if k * quarter > theta_start + 1e-9:
            angles.append(k * quarter)
        k += 1
    angles.append(theta_end)
    return angles


def connected_void_outline_geo(case):
    """(d) inner outline(s) of the section at z = 0 as explicit points, lines and arcs.

    No boolean is used: an up/down line is TANGENT to the hole it touches, and an OCC
    boolean on tangent faces left a sliver at the tangent line (first meshing pass: an STL
    vertex 29 um off the mesh, tetrahedra of quality 0.016).  Returns (geo_text, loops).
    """
    xn, yn = NOTCH_TIP_MM
    yc = yn + case["dy"]
    xl, xr = xn - case["dx"], xn + case["dx"]
    rl, rr = case["r_left"], case["r_right"]
    band = connected_band_mm(case)
    if band is None:
        return (f"""
// Hole outlines: two separate circles (no connecting line)
Circle(60) = {{{mm(xl)}, {mm(yc)}, 0, {mm(rl)}}};
Curve Loop(2) = {{60}};
Circle(61) = {{{mm(xr)}, {mm(yc)}, 0, {mm(rr)}}};
Curve Loop(3) = {{61}};
""", (2, 3))

    y0, y1 = band

    def half_chord(r, y):
        """Half chord of a hole at height y; 0 where the edge is tangent to the hole."""
        return math.sqrt(max(r * r - (y - yc) ** 2, 0.0))

    # corners where the holes meet the lower (y0) and upper (y1) edges of the line/throat
    corners = {62: (xl + half_chord(rl, y0), y0, "left hole meets the lower edge"),
               63: (xr - half_chord(rr, y0), y0, "right hole meets the lower edge"),
               64: (xr - half_chord(rr, y1), y1, "right hole meets the upper edge"),
               65: (xl + half_chord(rl, y1), y1, "left hole meets the upper edge")}

    # both arcs run counter-clockwise around the far side of their hole
    two_pi = 2 * math.pi
    start_r = math.atan2(corners[63][1] - yc, corners[63][0] - xr) % two_pi
    end_r = math.atan2(corners[64][1] - yc, corners[64][0] - xr) % two_pi
    if end_r <= start_r:
        end_r += two_pi
    start_l = math.atan2(corners[65][1] - yc, corners[65][0] - xl)
    end_l = math.atan2(corners[62][1] - yc, corners[62][0] - xl)
    if end_l <= start_l:
        end_l += two_pi

    points = [(60, xl, yc, "left hole centre"), (61, xr, yc, "right hole centre")]
    points += [(tag, x, y, note) for tag, (x, y, note) in corners.items()]
    curves = [(60, "Line", (62, 63), "lower edge of the line/throat")]
    state = {"point": 66, "curve": 61}

    def add_arc(centre_tag, xc, r, start_tag, end_tag, angles, label):
        """Arc pieces between consecutive break angles (each at most 90 degrees)."""
        tags = [start_tag]
        for a in angles[1:-1]:
            points.append((state["point"], xc + r * math.cos(a), yc + r * math.sin(a),
                           f"{label}, {math.degrees(a) % 360:.0f} deg"))
            tags.append(state["point"])
            state["point"] += 1
        tags.append(end_tag)
        for a, b in zip(tags[:-1], tags[1:]):
            curves.append((state["curve"], "Circle", (a, centre_tag, b), label))
            state["curve"] += 1

    add_arc(61, xr, rr, 63, 64, arc_break_angles(start_r, end_r), "right hole, far side")
    curves.append((state["curve"], "Line", (64, 65), "upper edge of the line/throat"))
    state["curve"] += 1
    add_arc(60, xl, rl, 65, 62, arc_break_angles(start_l, end_l), "left hole, far side")

    lines = ["", "// Hole pair + connecting line/throat as ONE closed outline cut from the section",
             "// (no boolean; points computed by connected_void_outline_geo in make_new_sample_geos.py)"]
    for tag, x, y, note in points:
        lines.append(f"Point({tag}) = {{{x * 1e-3:.15g}, {y * 1e-3:.15g}, 0, lc_refined}};  // {note}")
    for tag, kind, ends, note in curves:
        lines.append(f"{kind}({tag}) = {{{', '.join(str(e) for e in ends)}}};  // {note}")
    lines.append(f"Curve Loop(2) = {{{', '.join(str(c[0]) for c in curves)}}};")
    return "\n".join(lines) + "\n", (2,)


def sinewave_body_geo(case):
    p = case
    return f"""
//==============================================================================
// SINE-WAVE VOID
// In every z-slice the void is a circle of radius void_radius at height void_y whose
// centre is x_c(z) = x_axis + amp * Sin(2 Pi z / wavelength) (a disk translated along
// the sine, not a pipe swept normal to it).  Built as an OpenCASCADE loft through
// n_sections circles, extended z_ext beyond both faces so the cut is clean.
//==============================================================================
dy = {mm(p['dy'])};             // Void centre above the notch tip ({p['dy']:.1f} mm)
void_radius = {mm(p['radius'])};   // Void radius (D = {2 * p['radius']:.1f} mm)
x_axis = {mm(p['x_axis'])};        // Sine centreline ({p['x_axis']:.1f} mm; 14 = centered, 16 = one amplitude offset)
amp = {mm(p['amp'])};           // Sine amplitude ({p['amp']:.1f} mm)
wavelength = {mm(p['wavelength'])};    // One full period through the thickness ({p['wavelength']:.1f} mm)
n_sections = 48;
z_ext = 0.0005;

void_y = notch_tip_y + dy;
For k In {{0:n_sections}}
  zk = -z_ext + k * (extrude_z + 2 * z_ext) / n_sections;
  Circle(1000 + k) = {{x_axis + amp * Sin(2 * Pi * zk / wavelength), void_y, zk, void_radius}};
  Curve Loop(1000 + k) = {{1000 + k}};
EndFor
ThruSections(200) = {{1000:1000 + n_sections}};
BooleanDifference{{ Volume{{main_vol}}; Delete; }}{{ Volume{{200}}; Delete; }}
Printf("Sine void: x_c(z) = %g + %g sin(2 pi z / %g), y = %g, r = %g", x_axis, amp, wavelength, void_y, void_radius);
"""


def pillars_body_geo(case):
    p = case
    rows = len(p["row_z"])
    row_z = ", ".join(mm(z) for z in p["row_z"])
    row_x0 = ", ".join(mm(x) for x in p["row_x0"])
    row_n = ", ".join(str(n) for n in p["row_n"])
    return f"""
//==============================================================================
// PILLAR LAYER: lower block (with the notch) 0 <= y <= gap_y0, gap, upper block
// gap_y1 <= y <= beam_top; the blocks are joined only by cylindrical pillars (axis
// along y) on a staggered grid.  Pillars are embedded pillar_embed into both blocks
// so the union is a volume overlap, not a face contact.
//==============================================================================
gap_y0 = {mm(PILLAR_GAP_MM[0])};       // Top of the lower block (4.0 mm)
gap_y1 = {mm(PILLAR_GAP_MM[1])};      // Bottom of the upper block (4.5 mm)
beam_top = {mm(PILLAR_TOP_MM)};     // Top of the upper block = loading line (8.5 mm)
pillar_radius = {mm(PILLAR_RADIUS_MM)};  // 1 mm diameter
pillar_pitch = {mm(PILLAR_PITCH_MM)};   // x spacing within a row (1.5 mm)
pillar_embed = 0.0001;
n_rows = {rows};
row_z[] = {{{row_z}}};     // row positions through the thickness
row_x0[] = {{{row_x0}}};   // x of the first pillar in each row
row_n[] = {{{row_n}}};     // pillars per row

"""


def pillars_upper_outline_geo():
    """Upper block outline; defined before any extrusion so point 27 keeps its tag."""
    return """
// Upper block outline (loading point 27 on its top edge)
Point(41) = {      0, gap_y1, 0, lc};
Point(42) = { 0.0280, gap_y1, 0, lc};
Point(43) = { 0.0280, beam_top, 0, lc};
Point(27) = { 0.014, beam_top, 0, lc_support};  // Loading line position
Point(44) = {      0, beam_top, 0, lc};
Line(41) = {41, 42};
Line(42) = {42, 43};
Line(43) = {43, 27};     // Top right to loading line
Line(44) = {27, 44};     // Loading line to top left corner
Line(45) = {44, 41};
Curve Loop(2) = {41, 42, 43, 44, 45};
Plane Surface(2) = {2};
"""


def pillars_body_geo_tail():
    """Upper-block extrusion, pillars and the union of the three parts."""
    return """
out_top[] = Extrude {0, 0, extrude_z} {
  Surface{2};
};
top_vol = out_top[1];

pillars[] = {};
tag = 300;
For j In {0:n_rows-1}
  For i In {0:row_n[j]-1}
    Cylinder(tag) = {row_x0[j] + i * pillar_pitch, gap_y0 - pillar_embed, row_z[j],
                     0, gap_y1 - gap_y0 + 2 * pillar_embed, 0, pillar_radius};
    pillars[] += {tag};
    tag += 1;
  EndFor
EndFor
Printf("Pillars: %g", #pillars[]);
// Keep the seams: by default a union merges coplanar faces, which would erase the extruded
// support/loading edges at x = 4, 24 mm (y = 0) and x = 14 mm (top) that carry the line BCs.
Geometry.OCCUnionUnify = 0;
BooleanUnion{ Volume{main_vol}; Delete; }{ Volume{top_vol}; Volume{pillars[]}; Delete; }
"""



def size_fields_geo(case):
    """Box size fields in the meshfile_2holes layout (1 notch band, 3/7 feature, 4/5 supports)."""
    f = case["family"]
    top = mm(beam_top_mm(case))
    s = f"""
//==============================================================================
// MESH REFINEMENT FIELDS (tetrahedral elements)
//==============================================================================
// Field 1: Refine in the notch region (vertical band above the notch, full height)
Field[1] = Box;
Field[1].VIn = lc_refined;
Field[1].VOut = lc;
Field[1].XMin = 0.013;
Field[1].XMax = 0.015;
Field[1].YMin = 0;
Field[1].YMax = {top};
Field[1].ZMin = 0;
Field[1].ZMax = extrude_z;
Field[1].Thickness = 0.003;
"""
    if f in ("2holes_asym", "2holes_connected"):
        left_y = "hole_left_y" if f == "2holes_asym" else "hole_y"
        right_y = "hole_right_y" if f == "2holes_asym" else "hole_y"
        s += f"""
// Field 3: Refine between notch and LEFT hole (ligament region)
Field[3] = Box;
Field[3].VIn = lc_refined;
Field[3].VOut = lc;
Field[3].XMin = hole_left_x - hole_radius_left - 0.001;  // Extend past left hole
Field[3].XMax = notch_tip_x + 0.001;  // Extend right of notch
Field[3].YMin = notch_tip_y;  // Start at notch tip
Field[3].YMax = {left_y} + hole_radius_left + 0.001;  // Extend past hole top
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;

// Field 7: Refine between notch and RIGHT hole (ligament region)
Field[7] = Box;
Field[7].VIn = lc_refined;
Field[7].VOut = lc;
Field[7].XMin = notch_tip_x - 0.001;  // Extend left of notch
Field[7].XMax = hole_right_x + hole_radius_right + 0.001;  // Extend past right hole
Field[7].YMin = notch_tip_y;  // Start at notch tip
Field[7].YMax = {right_y} + hole_radius_right + 0.001;  // Extend past hole top
Field[7].ZMin = 0;
Field[7].ZMax = extrude_z;
Field[7].Thickness = 0.002;
"""
        field_list = "1, 3, 4, 5, 7"
    elif f == "sinewave":
        s += """
// Field 3: Refine from the notch tip to past the void, over the whole sine excursion
Field[3] = Box;
Field[3].VIn = lc_refined;
Field[3].VOut = lc;
Field[3].XMin = x_axis - amp - void_radius - 0.001;
Field[3].XMax = x_axis + amp + void_radius + 0.001;
Field[3].YMin = notch_tip_y;
Field[3].YMax = void_y + void_radius + 0.001;
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;
"""
        field_list = "1, 3, 4, 5"
    elif f == "pillars":
        s += """
// Field 3: Refine the pillar layer around the notch (x = 10 .. 18 mm)
Field[3] = Box;
Field[3].VIn = lc_refined;
Field[3].VOut = lc;
Field[3].XMin = 0.010;
Field[3].XMax = 0.018;
Field[3].YMin = gap_y0 - 0.0005;
Field[3].YMax = gap_y1 + 0.0005;
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;

// Field 8: Pillar layer over the full length (resolves every pillar)
Field[8] = Box;
Field[8].VIn = lc_pillar_layer;
Field[8].VOut = lc;
Field[8].XMin = 0;
Field[8].XMax = 0.028;
Field[8].YMin = gap_y0 - 0.00025;
Field[8].YMax = gap_y1 + 0.00025;
Field[8].ZMin = 0;
Field[8].ZMax = extrude_z;
Field[8].Thickness = 0.001;
"""
        field_list = "1, 3, 4, 5, 8"
    else:
        raise ValueError(f)
    s += f"""
// Field 4: Refine near left support line (x = 0.004)
Field[4] = Box;
Field[4].VIn = lc_support;
Field[4].VOut = lc;
Field[4].XMin = 0.002;
Field[4].XMax = 0.006;
Field[4].YMin = 0;
Field[4].YMax = 0.002;
Field[4].ZMin = 0;
Field[4].ZMax = extrude_z;
Field[4].Thickness = 0.002;

// Field 5: Refine near right support line (x = 0.024)
Field[5] = Box;
Field[5].VIn = lc_support;
Field[5].VOut = lc;
Field[5].XMin = 0.022;
Field[5].XMax = 0.026;
Field[5].YMin = 0;
Field[5].YMax = 0.002;
Field[5].ZMin = 0;
Field[5].ZMax = extrude_z;
Field[5].Thickness = 0.002;

// Field 6: Combine all fields (use minimum)
Field[6] = Min;
Field[6].FieldsList = {{{field_list}}};

Background Field = 6;

Printf("Line supports at x = 0.004 and x = 0.024 (y = 0); loading line at x = 0.014, y = {top}");
"""
    return s


def case_geo_text(case):
    """Assemble the full .geo text of one case."""
    f = case["family"]
    parts = [header_geo(case)]
    if f == "pillars":
        parts.append(f"lc_pillar_layer = {LC_PILLAR_LAYER:g};  // Mesh size over the whole pillar layer (125 microns)\n")
        parts.append(pillars_body_geo(case))           # parameters
        parts.append(notch_outline_geo(PILLAR_GAP_MM[0], with_load_point=False))
        parts.append(pillars_upper_outline_geo())
        parts.append(extrude_lower_geo())
        parts.append(pillars_body_geo_tail())
    elif f == "2holes_connected":
        parts.append(two_holes_connected_params_geo(case))
        inner_geo, inner_loops = connected_void_outline_geo(case)
        parts.append(notch_outline_geo(BEAM_HEIGHT_MM, with_load_point=True,
                                       inner_geo=inner_geo, inner_loops=inner_loops))
        parts.append(extrude_lower_geo())
    else:
        parts.append(notch_outline_geo(BEAM_HEIGHT_MM, with_load_point=True))
        parts.append(extrude_lower_geo())
        parts.append({"2holes_asym": two_holes_asym_body_geo,
                      "sinewave": sinewave_body_geo}[f](case))
    parts.append(physical_volume_geo())
    parts.append(size_fields_geo(case))
    return "".join(parts)


def write_case_geo(case):
    path = geo_path(case)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(case_geo_text(case))
    return path


# ----------------------------------------------------------------------------- checks
def notch_area_mm2(notch_round):
    """Cross-section area removed by the notch (mm^2)."""
    w = NOTCH_WIDTH_MM
    depth = NOTCH_TIP_MM[1]
    if notch_round:
        return w * (depth - w / 2) + math.pi * (w / 2) ** 2 / 2
    return w * depth


def connected_void_area_mm2(case):
    """(d) area of hole_left U hole_right U slot/throat, integrated over y (mm^2)."""
    rl, rr = case["r_left"], case["r_right"]
    area = math.pi * (rl ** 2 + rr ** 2)
    if holes_touch(case):
        h = throat_half_height_mm(case)
        band = (-h, h)
    else:
        slot = connected_slot_mm(case)
        if slot is None:
            return area
        yc = NOTCH_TIP_MM[1] + case["dy"]
        band = (slot[0] - yc, slot[1] - yc)
    span = 2 * case["dx"]                      # centre to centre
    n = 200000
    total = 0.0
    for k in range(n):                          # midpoint rule over the band
        t = band[0] + (k + 0.5) * (band[1] - band[0]) / n
        cover = (math.sqrt(max(rl ** 2 - t ** 2, 0.0)) + math.sqrt(max(rr ** 2 - t ** 2, 0.0)))
        total += max(span - cover, 0.0)
    return area + total * (band[1] - band[0]) / n


def expected_volume_mm3(case):
    """Analytic solid volume of the case (mm^3)."""
    t = BEAM_THICKNESS_MM
    notch = notch_area_mm2(NOTCH_ROUND)
    f = case["family"]
    if f == "pillars":
        lower = (BEAM_LENGTH_MM * PILLAR_GAP_MM[0] - notch) * t
        upper = BEAM_LENGTH_MM * (PILLAR_TOP_MM - PILLAR_GAP_MM[1]) * t
        n = sum(case["row_n"])
        pillars = n * math.pi * PILLAR_RADIUS_MM ** 2 * (PILLAR_GAP_MM[1] - PILLAR_GAP_MM[0])
        return lower + upper + pillars
    beam = (BEAM_LENGTH_MM * BEAM_HEIGHT_MM - notch) * t
    if f == "2holes_asym":
        return beam - 2 * math.pi * case["radius"] ** 2 * t
    if f == "2holes_connected":
        return beam - connected_void_area_mm2(case) * t
    if f == "sinewave":
        return beam - math.pi * case["radius"] ** 2 * t   # Cavalieri: every z-slice is a disk
    raise ValueError(f)


def check_geo_volume(case, gmsh):
    """Open the .geo (OCC, no meshing); return (n_volumes, occ_volume_mm3, expected_mm3)."""
    gmsh.clear()
    gmsh.open(geo_path(case))
    vols = gmsh.model.getEntities(3)
    vol = sum(gmsh.model.occ.getMass(3, tag) for _, tag in vols) * 1e9
    return len(vols), vol, expected_volume_mm3(case)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="open each .geo with gmsh and check its volume")
    args = ap.parse_args()

    for case in CASES:
        path = write_case_geo(case)
        print(f"wrote {os.path.relpath(path, TPB_DIR)}")
    if not args.check:
        return 0

    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    failures = 0
    for case in CASES:
        n_vol, vol, expected = check_geo_volume(case, gmsh)
        rel = (vol - expected) / expected
        ok = n_vol == 1 and abs(rel) < 2e-5
        failures += not ok
        print(f"{'OK ' if ok else 'BAD'} {mesh_basename(case):70s} volumes={n_vol} "
              f"V={vol:.4f} mm3 expected {expected:.4f} (rel {rel:+.1e})")
    gmsh.finalize()
    print(f"{len(CASES) - failures}/{len(CASES)} geometries pass")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
