SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc = 0.0002;  // Smaller values result in finer mesh
height = 0.06;
radius_outer = 0.0125;
radius_inner = 0.0016;

// Define the main cylinder with specific mesh size at base points
Cylinder(1) = {0, 0, 0, 0, 0, height, radius_outer, 2*Pi};

// Define the hole cylinder
Cylinder(2) = {0, 0, 0, 0, 0, height, radius_inner, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Set mesh size globally first
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 1;

// Apply mesh size to all points in the geometry
Characteristic Length{ PointsOf{ Volume{:}; } } = lc;

// Define Physical Volume for simulation or further meshing
Physical Volume("Pierced Cylinder") = {1};

// Label the surfaces - using automatic surface detection
// After BooleanDifference, we need to identify surfaces by their actual IDs
Physical Surface("Bottom") = {7};      // Bottom annular face
Physical Surface("Inner") = {4};       // Inner cylindrical surface
Physical Surface("Outer") = {5};       // Outer cylindrical surface
Physical Surface("Upper") = {6};       // Top annular face
