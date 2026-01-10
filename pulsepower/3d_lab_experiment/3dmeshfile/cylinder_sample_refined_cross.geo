SetFactory("OpenCASCADE");

// Mesh size parameters
lc_coarse = 0.001;     // Coarse mesh size for most of the domain
lc_refined = 0.0002;   // Refined mesh size in the two cross-pattern boxes

// Geometry parameters
height = 0.06;
radius_outer = 0.0125;
radius_inner = 0.0016;

// Define the main cylinder
Cylinder(1) = {0, 0, 0, 0, 0, height, radius_outer, 2*Pi};

// Define the hole cylinder
Cylinder(2) = {0, 0, 0, 0, 0, height, radius_inner, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Set coarse mesh size globally first
Mesh.CharacteristicLengthMin = lc_refined;
Mesh.CharacteristicLengthMax = lc_coarse;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 1;

// Apply coarse mesh size to all points initially
Characteristic Length{ PointsOf{ Volume{:}; } } = lc_coarse;

// Define Physical Volume for simulation
Physical Volume("Pierced Cylinder") = {1};

// Label the surfaces
Physical Surface("Bottom") = {7};      // Bottom annular face
Physical Surface("Inner") = {4};       // Inner cylindrical surface
Physical Surface("Outer") = {5};       // Outer cylindrical surface
Physical Surface("Upper") = {6};       // Top annular face

//==============================================================================
// MESH REFINEMENT FIELDS
//==============================================================================

// Box 1: Horizontal band (narrow in y-direction)
// bottom_left1 = '-1 -4e-4 0'
// top_right1 = '1 4e-4 0.06'
Field[1] = Box;
Field[1].VIn = lc_refined;     // Refined mesh inside box
Field[1].VOut = lc_coarse;     // Coarse mesh outside box
Field[1].XMin = -1.0;
Field[1].XMax = 1.0;
Field[1].YMin = -1.5e-3;
Field[1].YMax = 1.5e-3;
Field[1].ZMin = 0.0;
Field[1].ZMax = 0.06;
Field[1].Thickness = 0.001;    // Smooth transition zone

// Box 2: Vertical band (narrow in x-direction)
// bottom_left2 = '-4e-4 -1 0'
// top_right2 = '4e-4 1 0.06'
Field[2] = Box;
Field[2].VIn = lc_refined;     // Refined mesh inside box
Field[2].VOut = lc_coarse;     // Coarse mesh outside box
Field[2].XMin = -1.5e-3;
Field[2].XMax = 1.5e-3;
Field[2].YMin = -1.0;
Field[2].YMax = 1.0;
Field[2].ZMin = 0.0;
Field[2].ZMax = 0.06;
Field[2].Thickness = 0.001;    // Smooth transition zone

// Combine both refinement fields (use minimum mesh size)
Field[3] = Min;
Field[3].FieldsList = {1, 2};

// Set as background field
Background Field = 3;

Printf("Mesh refinement: Cross pattern with lc_refined = %.6f in two narrow bands", lc_refined);
Printf("  Box 1 (horizontal): x=[-1, 1], y=[-4e-4, 4e-4], z=[0, 0.06]");
Printf("  Box 2 (vertical):   x=[-4e-4, 4e-4], y=[-1, 1], z=[0, 0.06]");
Printf("  Coarse mesh elsewhere: lc_coarse = %.6f", lc_coarse);
