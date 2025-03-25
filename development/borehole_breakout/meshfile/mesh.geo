SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc = 0.001;  // Smaller values result in finer mesh

// Define the main cylinder with specific mesh size at base points
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, 0.075, lc};
Cylinder(1) = {0, 0, 0, 0, 0, 0.075, 0.025, 2*Pi};

// Define the hole cylinder
Point(7) = {-0.025, 0, 0.0375, lc};
Point(8) = {0.025, 0, 0.0375, lc};
Cylinder(2) = {-0.025, 0, 0.0375, 2, 0, 0, 0.00565, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Define Physical Volume for simulation or further meshing
Physical Volume("Pierced Cylinder") = {1};

Mesh.CharacteristicLengthExtendFromBoundary = 0; // Do not extend mesh size from boundaries
Mesh.CharacteristicLengthFromPoints = 1; // Use point-specific mesh sizes
Mesh.CharacteristicLengthFromCurvature = 1; // Allow mesh size to adapt based on curvature
Mesh.CharacteristicLengthMin = lc; // Ensure the minimum mesh size is set to the characteristic length
Mesh.CharacteristicLengthMax = lc; // Set the maximum mesh size to prevent larger elements