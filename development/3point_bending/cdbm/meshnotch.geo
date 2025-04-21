// Parameters
lc = 0.005;         // global mesh size

// Define square corner points
Point(1) = {       0,      0,   0, lc};
Point(2) = {   0.220,      0,   0, lc};
Point(3) = {   0.220,   0.05,   0, lc};
Point(4) = {   0.225,   0.05,   0, lc};
Point(5) = {   0.225,      0,   0, lc};
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