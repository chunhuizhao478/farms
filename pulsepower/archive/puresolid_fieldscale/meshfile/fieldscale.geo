SetFactory("OpenCASCADE");

// Define the rectangular box dimensions
length = 10; // Length in meters
width = 10;  // Width in meters
height = 0.05; // Height in meters

// Define the cylindrical hole dimensions
holeRadius = 0.1; // Radius in meters (diameter is 0.2m)
holeHeight = height; // The hole goes through the entire thickness of the box

// Mesh size
meshSize = 0.01;
meshSize_bc = 0.1;

// Create the box
Box(1) = {0, 0, 0, length, width, height};

// Create the cylinder
Cylinder(2) = {length/2, width/2, 0, 0, 0, holeHeight, holeRadius};

// Subtract the cylinder from the box
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// Define a point at the center of the cylindrical hole base
Point(3) = {length/2, width/2, 0, meshSize};

// Create a mesh field that varies the mesh size radially from the center of the cylindrical hole
Field[1] = Distance;
Field[1].NodesList = {3}; // Center of the cylinder base

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = meshSize;
Field[2].LcMax = meshSize_bc;
Field[2].DistMin = holeRadius; // Start increasing mesh size at the hole radius
Field[2].DistMax = 8; // Maximum distance to apply the gradient

Background Field = 2;

// Mark physical surfaces
Physical Surface("Top") = {10};        // Top surface of the box
Physical Surface("Bottom") = {12};     // Bottom surface of the box
Physical Surface("Front") = {9};      // Front surface of the box
Physical Surface("Back") = {11};       // Back surface of the box
Physical Surface("Left") = {8};       // Left surface of the box
Physical Surface("Right") = {13};      // Right surface of the box
Physical Surface("Cylinder") = {7};   // Surface of the cylindrical hole

// Mark the physical volume
Physical Volume("BoxWithHole") = {1};

