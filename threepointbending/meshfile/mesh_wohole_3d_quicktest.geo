SetFactory("OpenCASCADE");

//==============================================================================
// MESH PARAMETERS
//==============================================================================
// Global mesh size parameters
lc = 0.001;          // Global mesh size (1 mm)
lc_refined = 1e-4;   // Refined mesh size in notch region (50 microns)
lc_support = 2.5e-4; // Refined mesh size at support regions (250 microns)
extrude_z = 0.008;   // Total thickness in z-direction (8 mm)

// Element type: Tetrahedral elements
// Z-direction refinement: Field-based (MathEval)
// - At z=0.004 (center): element size = 0.5 mm (FINE)
// - At z=0, 0.008 (edges): element size = 1.0 mm (COARSE)
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
Point(23) = { 0.0135, 0.0080, 0, lc};  // Left boundary of loading patch
Point(24) = { 0.0145, 0.0080, 0, lc};  // Right boundary of loading patch

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

// Print extrusion output for debugging
Printf("=== Extrusion Results ===");
Printf("Volume: %g", out[0]);
Printf("Total output entities: %g", #out[]);

// Identify the surfaces corresponding to our boundary condition patches
// Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
// out[0] = volume
// out[1] = top surface at z=extrude_z
// out[2] = surface from Line 1
// out[3] = surface from Line 2 (LEFT SUPPORT PATCH!)
// out[4] = surface from Line 3
// ...
// out[10] = surface from Line 9 (RIGHT SUPPORT PATCH!)
// ...
// out[14] = surface from Line 13 (LOADING PATCH!)

left_support_surf = out[3];      // From Line 2 (x = 0.0035 to 0.0045, y = 0)
right_support_surf = out[10];    // From Line 9 (x = 0.0235 to 0.0245, y = 0)
loading_patch_surf = out[14];    // From Line 13 (x = 0.0135 to 0.0145, y = 0.008)

Printf("Left support surface: %g (from Line 2: x=0.0035-0.0045, y=0)", left_support_surf);
Printf("Right support surface: %g (from Line 9: x=0.0235-0.0245, y=0)", right_support_surf);
Printf("Loading patch surface: %g (from Line 13: x=0.0135-0.0145, y=0.008)", loading_patch_surf);

// Query all volumes in the geometry
all_vols[] = Volume "*";
Printf("Total volumes found: %g", #all_vols[]);
If (#all_vols[] > 0)
  For i In {0:#all_vols[]-1}
    Printf("  Volume[%g] = %g", i, all_vols[i]);
  EndFor
EndIf

// Define Physical Volume using the first (and only) volume
// This avoids the "unknown volume" warning by using the actual volume IDs
If (#all_vols[] > 0)
  Physical Volume("volume") = {all_vols[0]};
Else
  Printf("ERROR: No volumes found!");
EndIf
Physical Surface("bottom_left_support") = {left_support_surf};
Physical Surface("bottom_right_support") = {right_support_surf};
Physical Surface("top_loading") = {loading_patch_surf};

Printf("=== Physical Groups Created Successfully ===");
If (#all_vols[] > 0)
  Printf("Physical Volume 'volume': Volume %g", all_vols[0]);
EndIf
Printf("Physical Surface 'bottom_left_support': Surface %g", left_support_surf);
Printf("Physical Surface 'bottom_right_support': Surface %g", right_support_surf);
Printf("Physical Surface 'top_loading': Surface %g", loading_patch_surf);
Printf("===================================================");

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

// Field 2: Refine towards center in Z-direction using MathEval
// Element size increases with distance from z=0.004 (center)
z_center = extrude_z / 2.0;  // Center at z=0.004
Field[2] = MathEval;
Field[2].F = Sprintf("%.6f + %.6f * Fabs(z - %.6f) / %.6f",
                     lc_refined*10,  // Minimum size at center: 0.5mm
                     lc*0.5,         // Add up to 0.5mm based on distance
                     z_center,       // Center position
                     z_center);      // Normalize by half-thickness
// At z=0.004: size = 0.5mm
// At z=0 or z=0.008: size = 0.5mm + 0.5mm = 1.0mm

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
Field[6].FieldsList = {1, 2, 5, 8};

Background Field = 6;

Printf("Mesh refinement: Notch region (XY) + Z-center refinement + Bottom supports (distance-based)");
