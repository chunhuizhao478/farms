SetFactory("OpenCASCADE");

lc = 0.0025;

// Define main cylinder
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, 0.075, lc};
Cylinder(1) = {0, 0, 0, 0, 0, 0.075, 0.025, 2*Pi};

// Define parameters for the cuts
r_cut = 0.025;
h_cut = 0.005;

// Top cutting cylinder
Cylinder(2) = {0, 0, 0.075 - h_cut, 0, 0, h_cut, r_cut, 2*Pi};
// Bottom cutting cylinder
Cylinder(3) = {0, 0, 0, 0, 0, h_cut, r_cut, 2*Pi};

// Do both subtractions in one command
// NOTE: {Volume{2}; Volume{3};} means "subtract both Volume 2 AND Volume 3"
BooleanFragments{ Volume{1,2,3}; Delete; }{}

// Create a physical volume for the new, final volume
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// ------------------------------
// Add physical groups for surfaces
// You must determine these surface tags (e.g., by using the gmsh GUI)
// Example:
Physical Surface("Outer")  = {1, 4, 6};
Physical Surface("Top")    = {5};
Physical Surface("Bottom") = {7};
// ------------------------------

Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements