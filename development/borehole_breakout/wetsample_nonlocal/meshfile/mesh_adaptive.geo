SetFactory("OpenCASCADE");

lc_min = 0.0005;
lc_max = 0.001; // Maximum mesh size

// Define main cylinder
// Point(1) = {0, 0, 0, lc};
// Point(2) = {0, 0, 0.105, lc};
Cylinder(1) = {0, 0, 0, 0, 0, 0.105, 0.027, 2*Pi};

// Define parameters for the cuts
r_cut = 0.027;
h_cut = 0.005;

// Top cutting cylinder
Cylinder(2) = {0, 0, 0.105 - h_cut, 0, 0, h_cut, r_cut, 2*Pi};
// Bottom cutting cylinder
Cylinder(3) = {0, 0, 0, 0, 0, h_cut, r_cut, 2*Pi};

// Do both subtractions in one command
// NOTE: {Volume{2}; Volume{3};} means "subtract both Volume 2 AND Volume 3"
BooleanFragments{ Volume{1,2,3}; Delete; }{}

Cylinder(5) = {-0.035, 0, 0.0525, 2, 0, 0, 0.0055, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{:}; Delete; }{ Volume{5}; Delete; }

// Create a physical volume for the new, final volume
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// ------------------------------
// Add physical groups for surfaces
// You must determine these surface tags (e.g., by using the gmsh GUI)
// Example:
Physical Surface("Outer")  = {8};
Physical Surface("Inner")  = {9};
Physical Surface("Top")    = {5};
Physical Surface("Bottom") = {7};
// ------------------------------

// Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
// Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
// Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
// Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
// Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements

// --- New code: Use a Cylinder field for local mesh refinement ---
// This field imposes a fine (constant) mesh size within a cylindrical region
// around the hole (Cylinder(5)) and then transitions to a coarser mesh away.
Field[1] = Cylinder;
Field[1].XCenter = -0.035;         // Center of the cylinder (same as Cylinder(5))
Field[1].YCenter = 0;
Field[1].ZCenter = 0.0525;
Field[1].XAxis   = 1;              // Orientation: Cylinder(5) is created along the x-axis
Field[1].YAxis   = 0;
Field[1].ZAxis   = 0;
Field[1].Radius  = 0.015;         // Same radius as Cylinder(5)
Field[1].VIn     = lc_min;             // Fine mesh size inside the cylinder field
Field[1].VOut    = lc_max;           // Coarser mesh size outside the cylinder field
// Field[1].Thickness = 0.01;       // Thickness of the transition region

// Set the background field to use the cylinder field
Background Field = 1;