// Parameters
lc = 0.01;         // global mesh size
lc_refine = 0.001; // local mesh size for the fault

// Define square corner points
Point(1) = {       0,      0,   0, lc};
Point(2) = {   0.220,      0,   0, lc_refine};
Point(3) = {   0.220,   0.05,   0, lc_refine};
Point(4) = {   0.225,   0.05,   0, lc_refine};
Point(5) = {   0.225,      0,   0, lc_refine};
Point(6) = {   0.450,      0,   0, lc};
Point(7) = {   0.450,   0.10,   0, lc};
Point(8) = {       0,   0.10,   0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

Physical Surface("domain") = {1};

// 1. Create a Distance field from points that define the refined region
Field[1] = Distance;
Field[1].NodesList = {2,3,4,5};

// 2. Create a Threshold field that smoothly transitions the mesh size
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_refine;
Field[2].LcMax = lc;
Field[2].DistMin = 0.06;  // Adjust as needed
Field[2].DistMax = 0.2;  // Adjust as needed

// Set the Threshold field as the background field
Background Field = 2;