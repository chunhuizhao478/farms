SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc = 0.0005;  // Smaller values result in finer mesh

height = 0.001;
radius_outer = 0.0125;
radius_inner = 0.002;

// Define the main cylinder with specific mesh size at base points
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, height, lc};
Cylinder(1) = {0, 0, 0, 0, 0, height, radius_outer, 2*Pi};

// Define the hole cylinder
Point(7) = {0, 0, 0, lc};
Point(8) = {0, 0, height, lc};
Cylinder(2) = {0, 0, 0, 0, 0, height, radius_inner, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Define Physical Volume for simulation or further meshing
Physical Volume("Pierced Cylinder") = {1};

Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements

// Now label the surfaces accordingly
Physical Surface("Upper") = {6}; // Label the outer surface of the main cylinder
Physical Surface("Outer") = {5}; // Label the inner surface of the hole cylinder
Physical Surface("Inner") = {4};
Physical Surface("Bottom") = {7};

