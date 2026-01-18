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
dx = 0.002;          // Horizontal offset (1 mm to the right of notch tip)
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

// POINTS FOR THE THREE SURFACE PATCHES (1mm wide each)
// Bottom left support patch: x = 0.0035 to 0.0045
Point(19) = { 0.0035,      0, 0, lc};  // Left boundary of left support
Point(20) = { 0.0045,      0, 0, lc};  // Right boundary of left support

// Bottom right support patch: x = 0.0235 to 0.0245
Point(21) = { 0.0235,      0, 0, lc};  // Left boundary of right support
Point(22) = { 0.0245,      0, 0, lc};  // Right boundary of right support

// Top loading patch: x = 0.0135 to 0.0145
Point(23) = { 0.0139, 0.0080, 0, lc};  // Left boundary of loading patch
Point(24) = { 0.0141, 0.0080, 0, lc};  // Right boundary of loading patch

// Define square edges with subdivisions for the three patches
Line(1) = {1, 19};       // Bottom left corner to left support patch start
Line(2) = {19, 20};      // LEFT SUPPORT PATCH (1mm wide)
Line(3) = {20, 2};       // Left support patch end to notch start
Line(4) = {2, 3};        // Notch vertical
Circle(5) = {3, 14, 4};  // Notch circle left
Circle(6) = {4, 14, 5};  // Notch circle right
Line(7) = {5, 6};        // Notch vertical
Line(8) = {6, 21};       // Notch end to right support patch start
Line(9) = {21, 22};      // RIGHT SUPPORT PATCH (1mm wide)
Line(10) = {22, 7};      // Right support patch end to right corner
Line(11) = {7, 8};       // Right edge
Line(12) = {8, 24};      // Top right to loading patch right edge
Line(13) = {24, 23};     // LOADING PATCH (1mm wide)
Line(14) = {23, 9};      // Loading patch left edge to top left corner
Line(15) = {9, 1};       // Left edge

// Loop and surface for the square
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
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

// Find surfaces by bounding box after Boolean operation
// Surface In BoundingBox {xmin, ymin, zmin, xmax, ymax, zmax}
// Use small tolerance (eps) to capture the surfaces

eps = 1e-6;

// Left support patch: x = 0.0035 to 0.0045, y = 0, z = 0 to extrude_z
left_support_surfs[] = Surface In BoundingBox {
  0.0035-eps, -eps, -eps,
  0.0045+eps, eps, extrude_z+eps
};
Printf("Left support surfaces found: %g", #left_support_surfs[]);
If (#left_support_surfs[] > 0)
  Physical Surface("bottom_left_support") = {left_support_surfs[0]};
  Printf("Physical Surface 'bottom_left_support': Surface %g", left_support_surfs[0]);
EndIf

// Right support patch: x = 0.0235 to 0.0245, y = 0, z = 0 to extrude_z
right_support_surfs[] = Surface In BoundingBox {
  0.0235-eps, -eps, -eps,
  0.0245+eps, eps, extrude_z+eps
};
Printf("Right support surfaces found: %g", #right_support_surfs[]);
If (#right_support_surfs[] > 0)
  Physical Surface("bottom_right_support") = {right_support_surfs[0]};
  Printf("Physical Surface 'bottom_right_support': Surface %g", right_support_surfs[0]);
EndIf

// Top loading patch: x = 0.0135 to 0.0145, y = 0.008, z = 0 to extrude_z
loading_surfs[] = Surface In BoundingBox {
  0.0135-eps, 0.0080-eps, -eps,
  0.0145+eps, 0.0080+eps, extrude_z+eps
};
Printf("Loading surfaces found: %g", #loading_surfs[]);
If (#loading_surfs[] > 0)
  Physical Surface("top_loading") = {loading_surfs[0]};
  Printf("Physical Surface 'top_loading': Surface %g", loading_surfs[0]);
EndIf

Printf("===================================================");

// Store surface IDs for mesh refinement fields
left_support_surf = left_support_surfs[0];
right_support_surf = right_support_surfs[0];

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
Field[3].XMax = hole_center_x + hole_radius + 0.0005;  // Extend past hole
Field[3].YMin = notch_tip_y - 0.0005;  // Start at notch tip
Field[3].YMax = hole_center_y + hole_radius + 0.0005;  // Extend past hole top
Field[3].ZMin = 0;
Field[3].ZMax = extrude_z;
Field[3].Thickness = 0.002;  // Smooth transition thickness

// Field 4: Distance from bottom left support surface
Field[4] = Distance;
Field[4].SurfacesList = {left_support_surf};

// Field 5: Threshold for smooth transition around left support
Field[5] = Threshold;
Field[5].InField = 4;
Field[5].SizeMin = lc_support;   // Mesh size at surface (250 microns)
Field[5].SizeMax = lc;           // Mesh size far from surface (1 mm)
Field[5].DistMin = 0.001;        // Start transition at 1mm from surface
Field[5].DistMax = 0.003;        // End transition at 3mm from surface

// Field 7: Distance from bottom right support surface
Field[7] = Distance;
Field[7].SurfacesList = {right_support_surf};

// Field 8: Threshold for smooth transition around right support
Field[8] = Threshold;
Field[8].InField = 7;
Field[8].SizeMin = lc_support;   // Mesh size at surface (250 microns)
Field[8].SizeMax = lc;           // Mesh size far from surface (1 mm)
Field[8].DistMin = 0.001;        // Start transition at 1mm from surface
Field[8].DistMax = 0.003;        // End transition at 3mm from surface

// Field 6: Combine all fields (use minimum)
Field[6] = Min;
Field[6].FieldsList = {1, 3, 5, 8};

Background Field = 6;

Printf("Mesh refinement: Notch region (XY) + Bottom supports (distance-based)");
Printf("Geometry: 3D beam with cylindrical hole at (x=%g, y=%g), radius=%g", hole_center_x, hole_center_y, hole_radius);
