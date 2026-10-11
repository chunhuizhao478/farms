SetFactory("OpenCASCADE");

//==============================================================================
// (c) 2 holes not connected (same size, asymmetric offset)
// Case       : xL3.50mm_xR2.50mm_y3.6mm_radi1.0mm
// Experiment : 3-Alum_2mmD_3.6up_2.5-3.5mm_apt
// Source STL : New sample stls/(c) 2 holes not connected (Same size, asymmetric offset)/2mm_circle_3.6up_2.5-3.5mm apart v2.stl
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
Point(8) = { 0.0280, 0.008,   0, lc};
Point(9) = {      0, 0.008,   0, lc};

// POINTS FOR LINE SUPPORTS (at center of original patches)
Point(25) = { 0.004,      0, 0, lc_support};  // Left support line position
Point(26) = { 0.024,      0, 0, lc_support};  // Right support line position
Point(27) = { 0.014, 0.008, 0, lc_support};  // Loading line position

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
Line(10) = {8, 27};      // Top right to loading line
Line(11) = {27, 9};      // Loading line to top left corner
Line(12) = {9, 1};       // Left edge

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = {1};

// Extrude in z (tetrahedral elements; size controlled by the background field)
out[] = Extrude {0, 0, extrude_z} {
  Surface{1};
};
main_vol = out[1];

//==============================================================================
// TWO CYLINDRICAL HOLES (same size, asymmetric offsets from the notch tip)
//==============================================================================
dx_left = 0.0035;     // Left hole: horizontal distance from notch tip (3.50 mm)
dx_right = 0.0025;    // Right hole: horizontal distance from notch tip (2.50 mm)
dy_left = 0.0036;    // Left hole: vertical distance from notch tip (3.6 mm)
dy_right = 0.0036;   // Right hole: vertical distance from notch tip (3.6 mm)
hole_radius_left = 0.001;    // Left hole radius (1.0 mm)
hole_radius_right = 0.001;   // Right hole radius (1.0 mm)

hole_left_x = notch_tip_x - dx_left;
hole_left_y = notch_tip_y + dy_left;
hole_right_x = notch_tip_x + dx_right;
hole_right_y = notch_tip_y + dy_right;

Cylinder(100) = {hole_left_x, hole_left_y, 0, 0, 0, extrude_z, hole_radius_left};
Cylinder(101) = {hole_right_x, hole_right_y, 0, 0, 0, extrude_z, hole_radius_right};
BooleanDifference{ Volume{main_vol}; Delete; }{ Volume{100}; Volume{101}; Delete; }
Printf("Left hole at (%g, %g), r = %g; right hole at (%g, %g), r = %g",
       hole_left_x, hole_left_y, hole_radius_left, hole_right_x, hole_right_y, hole_radius_right);

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
Field[1].YMax = 0.008;
Field[1].ZMin = 0;
Field[1].ZMax = extrude_z;
Field[1].Thickness = 0.003;

// Field 3: Refine between notch and LEFT hole (ligament region)
Field[3] = Box;
Field[3].VIn = lc_refined;
Field[3].VOut = lc;
Field[3].XMin = hole_left_x - hole_radius_left - 0.001;  // Extend past left hole
Field[3].XMax = notch_tip_x + 0.001;  // Extend right of notch
Field[3].YMin = notch_tip_y;  // Start at notch tip
Field[3].YMax = hole_left_y + hole_radius_left + 0.001;  // Extend past hole top
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
Field[7].YMax = hole_right_y + hole_radius_right + 0.001;  // Extend past hole top
Field[7].ZMin = 0;
Field[7].ZMax = extrude_z;
Field[7].Thickness = 0.002;

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
Field[6].FieldsList = {1, 3, 4, 5, 7};

Background Field = 6;

Printf("Line supports at x = 0.004 and x = 0.024 (y = 0); loading line at x = 0.014, y = 0.008");
