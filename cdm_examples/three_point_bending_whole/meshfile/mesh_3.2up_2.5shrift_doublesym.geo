SetFactory("OpenCASCADE");

// Parameters
lc = 0.001;         // global mesh size
lc_refined = 5e-5; // refined mesh size

x0 = 0.0140; //notch tip coordinate x
y0 = 0.0016; //notch tip coordinate y

s1 = 0.0020; //shrift
u1 = 0.0032; //up
r1 = 0.001;  //radius

s2 = 0.0020; //shrift
u2 = 0.0032; //shrift
r2 = 0.001;  //radius

// Define square corner points
Point(1) = {      0,      0,   0, lc};
Point(2) = { 0.0135,      0,   0, lc};
Point(3) = { 0.0135, 0.0011,   0, lc};
Point(4) = {     x0,     y0,   0, lc};
Point(5) = { 0.0145, 0.0011,   0, lc};
Point(6) = { 0.0145,      0,   0, lc};
Point(7) = { 0.0280,      0,   0, lc};
Point(8) = { 0.0280, 0.0080,   0, lc};
Point(9) = {      0, 0.0080,   0, lc};

Point(10) = { x0 + s1 - r1  , y0 + u1       ,   0, lc};
Point(11) = { x0 + s1       , y0 + u1 + r1  ,   0, lc};
Point(12) = { x0 + s1 + r1  , y0 + u1       ,   0, lc};
Point(13) = { x0 + s1       , y0 + u1 - r1  ,   0, lc};

Point(14) = { 0.0140, 0.0011, 0, lc};
Point(15) = { x0 + s1, y0 + u1, 0, lc};

Point(16) = { 0.0040,      0, 0, lc};
Point(17) = { 0.0240,      0, 0, lc};
Point(18) = { 0.0140, 0.0080, 0, lc};

// second hole
Point(19) = { x0 - s2 + r2  , y0 + u2       ,   0, lc};
Point(20) = { x0 - s2       , y0 + u2 + r2  ,   0, lc};
Point(21) = { x0 - s2 - r2  , y0 + u2       ,   0, lc};
Point(22) = { x0 - s2       , y0 + u2 - r2  ,   0, lc};
Point(23) = { x0 - s2, y0 + u2, 0, lc};

// Define square edges
Line(1) = {1, 16};
Line(2) = {16, 2};
Line(3) = {2, 3};

Circle(4) = {3, 14, 4};
Circle(5) = {4, 14, 5};

Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 18};
Line(10) = {18, 9};
Line(11) = {9, 1};

// Loop and surface for the square
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(1) = {1};

// Define parameters for the hole
Circle(12) = {10, 15, 11};
Circle(13) = {11, 15, 12};
Circle(14) = {12, 15, 13};
Circle(15) = {13, 15, 10};

Curve Loop(2) = {12, 13, 14, 15};
Plane Surface(2) = {2};

// second hole
Circle(16) = {19, 23, 20};
Circle(17) = {20, 23, 21};
Circle(18) = {21, 23, 22};
Circle(19) = {22, 23, 19};
Curve Loop(3) = {16, 17, 18, 19};
Plane Surface(3) = {3};

// Subtract two hols from square to create a single surface
BooleanDifference { Surface{1}; Delete; } { Surface{2}; Delete; }
BooleanDifference { Surface{1}; Delete; } { Surface{3}; Delete; }

Field[1] = Box;
Field[1].VIn = lc_refined;  // Mesh size inside the fault zone
Field[1].VOut = lc;       // Mesh size outside the fault zone
Field[1].XMin = 0.010;
Field[1].XMax = 0.018;
Field[1].YMin = 0;
Field[1].YMax = 0.008;
Field[1].Thickness = 0.003;

Background Field = 1;
