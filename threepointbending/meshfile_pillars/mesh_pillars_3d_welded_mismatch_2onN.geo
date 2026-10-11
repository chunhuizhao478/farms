SetFactory("OpenCASCADE");

//==============================================================================
// (f) pillar structures (welded mismatch)
// Case       : welded_mismatch_2onN
// Experiment : 1-Alum_welded_mismatch_2onN_1
// Source STL : New sample stls/(f)  Pillar structures/Welded_mismatch_2onN v2.stl
// Written by threepointbending/new_samples_tools/make_new_sample_geos.py;
// parameters measured from the STL (mm in the STL, metres here).
//==============================================================================

//==============================================================================
// MESH PARAMETERS
//==============================================================================
// Global mesh size parameters
lc = 0.001;          // Global mesh size (1 mm)
lc_refined = 0.0001;   // Refined mesh size in notch/feature region (100 microns, as meshfile_2holes)
lc_support = 0.00025; // Refined mesh size at support regions (250 microns)
extrude_z = 0.008;   // Total thickness in z-direction (8 mm)
Mesh.MshFileVersion = 2.2;   // MSH 2.2 ASCII, as meshfile_height / meshfile_2holes

//==============================================================================
// NOTCH (tip at x = 14 mm, y = 1.6 mm; width 1 mm)
//==============================================================================
notch_tip_x = 0.014;
notch_tip_y = 0.0016;
// 1 = semicircular tip of radius 0.5 mm (all earlier meshfile_* sets, paper figure)
// 0 = flat tip, i.e. the rectangular 1.0 x 1.6 mm notch exactly as drawn in the STL
notch_round = 1;
lc_pillar_layer = 0.000125;  // Mesh size over the whole pillar layer (125 microns)

//==============================================================================
// PILLAR LAYER: lower block (with the notch) 0 <= y <= gap_y0, gap, upper block
// gap_y1 <= y <= beam_top; the blocks are joined only by cylindrical pillars (axis
// along y) on a staggered grid.  Pillars are embedded pillar_embed into both blocks
// so the union is a volume overlap, not a face contact.
//==============================================================================
gap_y0 = 0.004;       // Top of the lower block (4.0 mm)
gap_y1 = 0.0045;      // Bottom of the upper block (4.5 mm)
beam_top = 0.0085;     // Top of the upper block = loading line (8.5 mm)
pillar_radius = 0.0005;  // 1 mm diameter
pillar_pitch = 0.0015;   // x spacing within a row (1.5 mm)
pillar_embed = 0.0001;
n_rows = 5;
row_z[] = {0.001, 0.0025, 0.004, 0.0055, 0.007};     // row positions through the thickness
row_x0[] = {0.00125, 0.002, 0.00125, 0.002, 0.00125};   // x of the first pillar in each row
row_n[] = {18, 17, 18, 17, 18};     // pillars per row


//==============================================================================
// BEAM OUTLINE (XY plane at z = 0), extruded in z below
//==============================================================================
Point(1) = {      0,      0,   0, lc};
Point(2) = { 0.0135,      0,   0, lc};
If (notch_round == 1)
  Point(3) = { 0.0135, 0.0011,   0, lc};
  Point(4) = { 0.0140, 0.0016,   0, lc};
  Point(5) = { 0.0145, 0.0011,   0, lc};
  Point(14) = { 0.0140, 0.0011, 0, lc};  // Circle center for notch
Else
  Point(3) = { 0.0135, 0.0016,   0, lc};
  Point(4) = { 0.0140, 0.0016,   0, lc};
  Point(5) = { 0.0145, 0.0016,   0, lc};
EndIf
Point(6) = { 0.0145,      0,   0, lc};
Point(7) = { 0.0280,      0,   0, lc};
Point(8) = { 0.0280, 0.004,   0, lc};
Point(9) = {      0, 0.004,   0, lc};

// POINTS FOR LINE SUPPORTS (at center of original patches)
Point(25) = { 0.004,      0, 0, lc_support};  // Left support line position
Point(26) = { 0.024,      0, 0, lc_support};  // Right support line position

// Bottom edge: left corner -> left support -> notch -> right support -> right corner
Line(1) = {1, 25};       // Bottom left corner to left support line
Line(2) = {25, 2};       // Left support line to notch start
Line(3) = {2, 3};        // Notch vertical left
If (notch_round == 1)
  Circle(4) = {3, 14, 4};  // Notch circle left
  Circle(5) = {4, 14, 5};  // Notch circle right
Else
  Line(4) = {3, 4};        // Flat notch tip, left half
  Line(5) = {4, 5};        // Flat notch tip, right half
EndIf
Line(6) = {5, 6};        // Notch vertical right
Line(7) = {6, 26};       // Notch end to right support line
Line(8) = {26, 7};       // Right support line to right corner
Line(9) = {7, 8};        // Right edge
Line(10) = {8, 9};       // Top edge
Line(12) = {9, 1};       // Left edge

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12};
Plane Surface(1) = {1};

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

// Extrude in z (tetrahedral elements; size controlled by the background field)
out[] = Extrude {0, 0, extrude_z} {
  Surface{1};
};
main_vol = out[1];

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

// Physical Volume = MOOSE block 1 (the decks' damage block); every volume of the result
all_vols[] = Volume "*";
Printf("Total volumes found: %g", #all_vols[]);
Physical Volume("volume") = {all_vols[]};

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
Field[1].YMax = 0.0085;
Field[1].ZMin = 0;
Field[1].ZMax = extrude_z;
Field[1].Thickness = 0.003;

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
Field[6].FieldsList = {1, 3, 4, 5, 8};

Background Field = 6;

Printf("Line supports at x = 0.004 and x = 0.024 (y = 0); loading line at x = 0.014, y = 0.0085");
