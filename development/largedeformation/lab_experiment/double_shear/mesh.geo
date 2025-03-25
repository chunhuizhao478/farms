// -----------------------------------------------------------------------------
// Gmsh script: sample_mesh.geo
// This script creates three rectangular regions stacked vertically.
// Each rectangle is defined as a separate Plane Surface for easy subdomain
// assignment (Physical Groups).
// -----------------------------------------------------------------------------

SetFactory("OpenCASCADE"); // Use the OpenCASCADE kernel
l = 0.0005; 

// ---------------------------
// 1) Bottom rectangle
//    Coordinates: (0, 0) to (0.05, 0.004)
// ---------------------------
Point(1) = {0.0,    0.0,   0.0, l};
Point(2) = {0.05,   0.0,   0.0, l};
Point(3) = {0.05,   0.004, 0.0, l};
Point(4) = {0.0,    0.004, 0.0, l};

// ---------------------------
// 2) Middle rectangle
//    Coordinates: (0.01, 0.004) to (0.04, 0.006)
// ---------------------------
Point(5)  = {0.01,  0.004, 0.0, l};
Point(6)  = {0.04,  0.004, 0.0, l};
Point(7)  = {0.04,  0.006, 0.0, l};
Point(8)  = {0.01,  0.006, 0.0, l};

// ---------------------------
// 3) Top rectangle
//    Coordinates: (0, 0.006) to (0.05, 0.01)
// ---------------------------
Point(9)  = {0.0,   0.006, 0.0, l};
Point(10) = {0.05,  0.006, 0.0, l};
Point(11) = {0.05,  0.01,  0.0, l};
Point(12) = {0.0,   0.01,  0.0, l};

// --------------------------------------------------
// Define lines connecting the points
// --------------------------------------------------
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 6};
Line(4) = {6, 7};
Line(5) = {7, 10};
Line(6) = {10, 11};
Line(7) = {11, 12};
Line(8) = {12, 9};
Line(9) = {9, 8};
Line(10) = {8, 5};
Line(11) = {5, 4};
Line(12) = {4, 1};

// --------------------------------------------------
// Define the rectangular regions
// --------------------------------------------------
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(2) = {1};

// --------------------------------------------------
// Embed lines
// --------------------------------------------------
Line(13) = {5, 6};
Curve{13} In Surface{2};
Line(14) = {7, 8};
Curve{14} In Surface{2};
