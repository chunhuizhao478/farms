SetFactory("OpenCASCADE");

lc = 0.001;

// Define main cylinder
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, 0.105, lc};
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
Physical Surface("Top")    = {5};
Physical Surface("Bottom") = {7};
// ------------------------------

Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements


// SetFactory("OpenCASCADE");

// // Characteristic length (mesh size)
// lc = 0.001;  // Smaller values result in finer mesh

// // Define the main cylinder with specific mesh size at base points
// Point(1) = {0, 0, 0, lc};
// Point(2) = {0, 0, 0.075, lc};
// Cylinder(1) = {0, 0, 0, 0, 0, 0.075, 0.025, 2*Pi};

// // Define the hole cylinder
// Point(7) = {-0.025, 0, 0.0375, lc};
// Point(8) = {0.025, 0, 0.0375, lc};
// Cylinder(2) = {-0.025, 0, 0.0375, 2, 0, 0, 0.00565, 2*Pi};

// // Boolean operation to subtract the hole from the main cylinder
// BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// // Define Physical Volume for simulation or further meshing
// Physical Volume("Pierced Cylinder") = {1};

// Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
// Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
// Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
// Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
// Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements