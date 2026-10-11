SetFactory("OpenCASCADE");

//==============================================================================
// (e) sine wave hole
// Case       : D1.0mm_y1.6mm_amp2.0mm_1ampoffset
// Experiment : Offset sine wave
// Source STL : New sample stls/(e) Sine wave hole/Alum_1mmD_1.6up_sinewavevoid_Amp2_1ampoffset.stl
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
// SINE-WAVE VOID
// In every z-slice the void is a circle of radius void_radius at height void_y whose
// centre is x_c(z) = x_axis + amp * Sin(2 Pi z / wavelength) (a disk translated along
// the sine, not a pipe swept normal to it).  Built as an OpenCASCADE loft through
// n_sections circles, extended z_ext beyond both faces so the cut is clean.
//==============================================================================
dy = 0.0016;             // Void centre above the notch tip (1.6 mm)
void_radius = 0.0005;   // Void radius (D = 1.0 mm)
x_axis = 0.016;        // Sine centreline (16.0 mm; 14 = centered, 16 = one amplitude offset)
amp = 0.002;           // Sine amplitude (2.0 mm)
wavelength = 0.008;    // One full period through the thickness (8.0 mm)
n_sections = 48;
z_ext = 0.0005;

void_y = notch_tip_y + dy;
For k In {0:n_sections}
  zk = -z_ext + k * (extrude_z + 2 * z_ext) / n_sections;
  Circle(1000 + k) = {x_axis + amp * Sin(2 * Pi * zk / wavelength), void_y, zk, void_radius};
  Curve Loop(1000 + k) = {1000 + k};
EndFor
ThruSections(200) = {1000:1000 + n_sections};
BooleanDifference{ Volume{main_vol}; Delete; }{ Volume{200}; Delete; }
Printf("Sine void: x_c(z) = %g + %g sin(2 pi z / %g), y = %g, r = %g", x_axis, amp, wavelength, void_y, void_radius);

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
Field[6].FieldsList = {1, 3, 4, 5};

Background Field = 6;

Printf("Line supports at x = 0.004 and x = 0.024 (y = 0); loading line at x = 0.014, y = 0.008");
