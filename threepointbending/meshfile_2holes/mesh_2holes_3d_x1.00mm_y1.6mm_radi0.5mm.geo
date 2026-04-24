SetFactory("OpenCASCADE");

//==============================================================================
// MESH PARAMETERS
//==============================================================================
// Global mesh size parameters
lc = 0.001;          // Global mesh size (1 mm)
lc_refined = 1e-4;   // Refined mesh size in notch region (50 microns)
lc_support = 2.5e-4; // Refined mesh size at support regions (250 microns)
extrude_z = 0.008;   // Total thickness in z-direction (8 mm)

//==============================================================================
// HOLE PARAMETERS (position relative to notch tip)
//==============================================================================
// Notch tip coordinates (reference point)
notch_tip_x = 0.0140;
notch_tip_y = 0.0016;

//==============================================================================
// TWO HOLES CONFIGURATION
// Hole 1: Left side of notch tip (negative x offset)
// Hole 2: Right side of notch tip (positive x offset)
//==============================================================================

// Horizontal offset from notch tip (absolute distance in meters)
// User can modify these values to change hole positions
dx_left = 0.001;     // Left hole: horizontal distance from notch tip (1.0 mm)
dx_right = 0.001;    // Right hole: horizontal distance from notch tip (1.0 mm)

// Vertical offset from notch tip (same for both holes by default)
dy_left = 0.0016;    // Left hole: vertical distance from notch tip (1.6 mm)
dy_right = 0.0016;   // Right hole: vertical distance from notch tip (1.6 mm)

// Hole geometry (can be different for each hole)
hole_radius_left = 0.0005;   // Left hole radius (0.5 mm)
hole_radius_right = 0.0005;  // Right hole radius (0.5 mm)

// Computed hole center coordinates
// Left hole: located at notch_tip_x - dx_left (to the left of notch)
hole_left_x = notch_tip_x - dx_left;
hole_left_y = notch_tip_y + dy_left;
hole_left_z = extrude_z / 2.0;

// Right hole: located at notch_tip_x + dx_right (to the right of notch)
hole_right_x = notch_tip_x + dx_right;
hole_right_y = notch_tip_y + dy_right;
hole_right_z = extrude_z / 2.0;

// Element type: Tetrahedral elements
// Geometry: 3D beam with notch and two cylindrical holes (one on each side of notch)
// Hole positions are parameterized by (dx_left, dy_left) and (dx_right, dy_right)
// XY-plane refinement: 50 microns in notch region, 1 mm elsewhere
// Refinement is controlled by background mesh fields, not layer distribution
//==============================================================================

// Define square corner points (2D base in XY plane at z=0)
Point(1) = {      0,      0,   0, lc};
Point(2) = { 0.0135,      0,   0, lc};
Point(3) = { 0.0135, 0.0011,   0, lc};
Point(4) = { 0.0140, 0.0016,   0, lc};
Point(5) = { 0.0145, 0.0011,   0, lc};
Point(6) = { 0.0145,      0,   0, lc};
Point(7) = { 0.0280,      0,   0, lc};
Point(8) = { 0.0280, 0.0080,   0, lc};
Point(9) = {      0, 0.0080,   0, lc};
Point(14) = { 0.0140, 0.0011, 0, lc};  // Circle center for notch

// POINTS FOR LINE SUPPORTS (at center of original patches)
// Left support line: x = 0.004 (center)
Point(25) = { 0.004,      0, 0, lc_support};  // Left support line position

// Right support line: x = 0.024 (center)
Point(26) = { 0.024,      0, 0, lc_support};  // Right support line position

// Top loading line: x = 0.014 (center)
Point(27) = { 0.014, 0.0080, 0, lc_support};  // Loading line position

// Define edges with LINE SUPPORT points (simplified from surface patches)
// Bottom edge: left corner -> left support -> notch -> right support -> right corner
Line(1) = {1, 25};       // Bottom left corner to left support line
Line(2) = {25, 2};       // Left support line to notch start
Line(3) = {2, 3};        // Notch vertical left
Circle(4) = {3, 14, 4};  // Notch circle left
Circle(5) = {4, 14, 5};  // Notch circle right
Line(6) = {5, 6};        // Notch vertical right
Line(7) = {6, 26};       // Notch end to right support line
Line(8) = {26, 7};       // Right support line to right corner

// Right edge
Line(9) = {7, 8};

// Top edge: right corner -> loading line -> left corner
Line(10) = {8, 27};      // Top right to loading line
Line(11) = {27, 9};      // Loading line to top left corner

// Left edge
Line(12) = {9, 1};

// Loop and surface for the square
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = {1};

// Extrude the surface in z-direction to create 3D volume
// Using tetrahedral elements
// Z-direction refinement will be controlled by background mesh field

out[] = Extrude {0, 0, extrude_z} {
  Surface{1};
  // No Recombine -> creates tetrahedral elements
};

Printf("Z-direction mesh: Tetrahedral elements with field-based refinement");

// Query the actual volume ID after extrusion (before creating cylinders)
vols_after_extrude[] = Volume "*";
Printf("Volume after extrusion: %g (count: %g)", vols_after_extrude[0], #vols_after_extrude[]);
main_vol = vols_after_extrude[0];

//==============================================================================
// CREATE TWO CYLINDRICAL HOLES
//==============================================================================
// Create cylinders through the entire thickness
// Cylinder(tag) = {x, y, z, dx, dy, dz, r} - creates cylinder from (x,y,z) along direction (dx,dy,dz)

// Left hole (to the left of notch tip)
Cylinder(100) = {hole_left_x, hole_left_y, 0, 0, 0, extrude_z, hole_radius_left};

// Right hole (to the right of notch tip)
Cylinder(101) = {hole_right_x, hole_right_y, 0, 0, 0, extrude_z, hole_radius_right};

Printf("=== Two Holes Configuration ===");
Printf("Left hole center: (%g, %g), radius: %g", hole_left_x, hole_left_y, hole_radius_left);
Printf("Left hole offset from notch tip: dx=%g, dy=%g", -dx_left, dy_left);
Printf("Right hole center: (%g, %g), radius: %g", hole_right_x, hole_right_y, hole_radius_right);
Printf("Right hole offset from notch tip: dx=%g, dy=%g", dx_right, dy_right);

// Subtract both cylinders from the main volume using actual volume ID
Printf("Subtracting cylinders (Volume 100, 101) from main volume (Volume %g)", main_vol);
BooleanDifference{ Volume{main_vol}; Delete; }{ Volume{100}; Volume{101}; Delete; }

// Print extrusion output for debugging
Printf("=== After Boolean Difference (two holes subtraction) ===");

// Query all volumes after Boolean operation (volume ID changes after BooleanDifference)
all_vols[] = Volume "*";
Printf("Total volumes found: %g", #all_vols[]);
If (#all_vols[] > 0)
  For i In {0:#all_vols[]-1}
    Printf("  Volume[%g] = %g", i, all_vols[i]);
  EndFor
EndIf

// Query all surfaces to find the boundary patches
// After BooleanDifference, surface IDs change, so we query by geometric location
all_surfs[] = Surface "*";
Printf("Total surfaces found: %g", #all_surfs[]);

// Define Physical Volume using the volume after Boolean operation
If (#all_vols[] > 0)
  Physical Volume("volume") = {all_vols[0]};
  Printf("Physical Volume 'volume': Volume %g", all_vols[0]);
Else
  Printf("ERROR: No volumes found!");
EndIf

//==============================================================================
// BOUNDARY DEFINITIONS
// Note: Line/point boundaries are created using MOOSE mesh generators
// Only Physical Volume is needed here
//==============================================================================

eps = 1e-6;

Printf("===================================================");
Printf("Line support positions (defined in MOOSE input file):");
Printf("  Left support line at x = 0.004");
Printf("  Right support line at x = 0.024");
Printf("  Loading line at x = 0.014, y = 0.008");

// Mesh refinement fields for tetrahedral elements
// Field 1: Refine in the notch region (XY plane)
Field[1] = Box;
Field[1].VIn = lc_refined;  // Mesh size inside the notch region (50 microns)
Field[1].VOut = lc;         // Mesh size outside (1 mm)
Field[1].XMin = 0.013;
Field[1].XMax = 0.015;
Field[1].YMin = 0;
Field[1].YMax = 0.008;
Field[1].ZMin = 0;
Field[1].ZMax = extrude_z;
Field[1].Thickness = 0.003;

// Field 3: Refine between notch and LEFT hole (ligament region)
Field[3] = Box;
Field[3].VIn = lc_refined;  // Same refined mesh size as notch region
Field[3].VOut = lc;         // Mesh size outside
Field[3].XMin = hole_left_x - hole_radius_left - 0.001;  // Extend past left hole
Field[3].XMax = notch_tip_x + 0.001;  // Extend right of notch
Field[3].YMin = notch_tip_y;  // Start at notch tip
Field[3].YMax = hole_left_y + hole_radius_left + 0.001;  // Extend past hole top
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;  // Smooth transition thickness

// Field 7: Refine between notch and RIGHT hole (ligament region)
Field[7] = Box;
Field[7].VIn = lc_refined;  // Same refined mesh size as notch region
Field[7].VOut = lc;         // Mesh size outside
Field[7].XMin = notch_tip_x - 0.001;  // Extend left of notch
Field[7].XMax = hole_right_x + hole_radius_right + 0.001;  // Extend past right hole
Field[7].YMin = notch_tip_y;  // Start at notch tip
Field[7].YMax = hole_right_y + hole_radius_right + 0.001;  // Extend past hole top
Field[7].ZMin = 0;
Field[7].ZMax = extrude_z;
Field[7].Thickness = 0.002;  // Smooth transition thickness

// Field 4: Refine near left support line (x = 0.004)
Field[4] = Box;
Field[4].VIn = lc_support;  // Mesh size at support (250 microns)
Field[4].VOut = lc;         // Mesh size outside (1 mm)
Field[4].XMin = 0.002;
Field[4].XMax = 0.006;
Field[4].YMin = 0;
Field[4].YMax = 0.002;
Field[4].ZMin = 0;
Field[4].ZMax = extrude_z;
Field[4].Thickness = 0.002;

// Field 5: Refine near right support line (x = 0.024)
Field[5] = Box;
Field[5].VIn = lc_support;  // Mesh size at support (250 microns)
Field[5].VOut = lc;         // Mesh size outside (1 mm)
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

Printf("Mesh refinement: Notch region (XY) + Both hole ligaments + Bottom supports");
Printf("Geometry: 3D beam with TWO cylindrical holes");
Printf("  Left hole at (x=%g, y=%g), radius=%g", hole_left_x, hole_left_y, hole_radius_left);
Printf("  Right hole at (x=%g, y=%g), radius=%g", hole_right_x, hole_right_y, hole_radius_right);
