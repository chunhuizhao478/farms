lc = 100;
lc2 = 2000;

// Define domain boundary points (using coarser mesh size)
Point(1) = { -52000, -40000, 0.0, lc2};
Point(2) = {  48000, -40000, 0.0, lc2};
Point(3) = {  48000,  40000, 0.0, lc2};
Point(4) = { -52000,  40000, 0.0, lc2};

// Points for the refined region (using finer mesh size)
Point(5) = { -16000,    0.0, 0.0, lc};
Point(6) = {  12000,    0.0, 0.0, lc};
Point(7) = {      0,    0.0, 0.0, lc}; // along main fault
Point(8) = { 10392.3, -6000, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,7};
Line(6) = {7,6};
Line(7) = {7,8};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

// Ensure internal points/lines are respected on the surface
Point{5,6,7,8} In Surface{2};
Line{5,6,7} In Surface{2};

// ***** Smoother Mesh Transition Setup *****

// 1. Create a Distance field from points that define the refined region
Field[1] = Distance;
Field[1].NodesList = {5, 6, 7, 8};

// 2. Create a Threshold field that smoothly transitions the mesh size
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc;
Field[2].LcMax = lc2;
Field[2].DistMin = 1000;  // Adjust as needed
Field[2].DistMax = 5000;  // Adjust as needed

// Set the Threshold field as the background field
Background Field = 2;

// ***** End of Smoother Mesh Transition Setup *****

// Define physical groups for boundaries and surfaces
Physical Curve("180") = {5,6};
Physical Curve("30")  = 7;

Physical Curve("bottom") = 1;
Physical Curve("right")  = 2;
Physical Curve("top")    = 3;
Physical Curve("left")   = 4;

Physical Surface("100") = 2;

// Choose the meshing algorithm (for example, algorithm 6 for standard triangulation)
Mesh.Algorithm = 6;
// For quadrilateral meshing, you might use algorithm 11 with recombination:
// Mesh.Algorithm = 11;
// Mesh.RecombinationAlgorithm = 2;

// For the v4 option with enlarged domain, ensure that the background field (Field[2]) smoothly transitions as desired.