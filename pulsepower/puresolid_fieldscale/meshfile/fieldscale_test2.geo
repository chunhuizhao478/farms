SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc = 0.0005;  // Smaller values result in finer mesh
lc2 = 0.005; // Larger mesh size for outer surfaces

height = 0.01;
radius_outer = 0.125;
radius_inner = 0.02;

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

// Define a point at the center of the cylindrical hole base
Point(3) = {0, 0, 0, lc};

// Create a mesh field that varies the mesh size radially from the center of the cylindrical hole
Field[1] = Distance;
Field[1].NodesList = {3}; // Center of the cylinder base

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc;
Field[2].LcMax = lc2;
Field[2].DistMin = radius_inner; // Start increasing mesh size at the hole radius
Field[2].DistMax = radius_outer; // Maximum distance to apply the gradient

Background Field = 2;

// Now label the surfaces accordingly
Physical Surface("Upper") = {6}; // Label the outer surface of the main cylinder
Physical Surface("Outer") = {5}; // Label the inner surface of the hole cylinder
Physical Surface("Inner") = {4};
Physical Surface("Bottom") = {7};

