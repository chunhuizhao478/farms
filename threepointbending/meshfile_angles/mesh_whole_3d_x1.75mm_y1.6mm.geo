SetFactory("OpenCASCADE");

//==============================================================================
// MESH PARAMETERS
//==============================================================================
// Global mesh size parameters
lc = 0.001;          // Global mesh size (1 mm)
lc_refined = 5e-5;   // Refined mesh size in notch region (50 microns)
lc_support = 2.5e-4; // Refined mesh size at support regions (250 microns)
extrude_z = 0.008;   // Total thickness in z-direction (8 mm)

//==============================================================================
// HOLE PARAMETERS (position relative to notch tip)
//==============================================================================
// Notch tip coordinates (reference point)
notch_tip_x = 0.0140;
notch_tip_y = 0.0016;

// Hole position offset from notch tip
// dx = hole_center_x - notch_tip_x
// dy = hole_center_y - notch_tip_y
dx = 0.00175;          // Horizontal offset (1.75 mm from notch tip)
dy = 0.0016;         // Vertical offset (1.6 mm above notch tip)

// Hole geometry
hole_radius = 0.001; // Hole radius (1 mm)

// Computed hole center coordinates
hole_center_x = notch_tip_x + dx;
hole_center_y = notch_tip_y + dy;
hole_center_z = extrude_z / 2.0;  // Center in z-direction

// Element type: Tetrahedral elements
// Geometry: 3D beam with notch and cylindrical hole above notch
// Hole position is parameterized by (dx, dy) offset from notch tip
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

// Query the actual volume ID after extrusion (before creating cylinder)
vols_after_extrude[] = Volume "*";
Printf("Volume after extrusion: %g (count: %g)", vols_after_extrude[0], #vols_after_extrude[]);
main_vol = vols_after_extrude[0];

//==============================================================================
// CREATE CYLINDRICAL HOLE
//==============================================================================
// Create a cylinder through the entire thickness
// Cylinder(tag) = {x, y, z, dx, dy, dz, r} - creates cylinder from (x,y,z) along direction (dx,dy,dz)
Cylinder(100) = {hole_center_x, hole_center_y, 0, 0, 0, extrude_z, hole_radius};

Printf("Hole center: (%g, %g), radius: %g", hole_center_x, hole_center_y, hole_radius);
Printf("Hole position relative to notch tip: dx=%g, dy=%g", dx, dy);

// Subtract the cylinder from the main volume using actual volume ID
Printf("Subtracting cylinder (Volume 100) from main volume (Volume %g)", main_vol);
BooleanDifference{ Volume{main_vol}; Delete; }{ Volume{100}; Delete; }

// Print extrusion output for debugging
Printf("=== After Boolean Difference (hole subtraction) ===");

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

// Field 3: Refine between notch and hole (ligament region)
// Region from notch tip to hole bottom edge
Field[3] = Box;
Field[3].VIn = lc_refined;  // Same refined mesh size as notch region
Field[3].VOut = lc;         // Mesh size outside
Field[3].XMin = notch_tip_x - 0.001;  // Extend left of notch
Field[3].XMax = hole_center_x + hole_radius + 0.001;  // Extend past hole
Field[3].YMin = notch_tip_y;  // Start at notch tip
Field[3].YMax = hole_center_y + hole_radius + 0.001;  // Extend past hole top
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;  // Smooth transition thickness

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
Field[6].FieldsList = {1, 3, 4, 5};

Background Field = 6;

Printf("Mesh refinement: Notch region (XY) + Bottom supports (distance-based)");
Printf("Geometry: 3D beam with cylindrical hole at (x=%g, y=%g), radius=%g", hole_center_x, hole_center_y, hole_radius);
