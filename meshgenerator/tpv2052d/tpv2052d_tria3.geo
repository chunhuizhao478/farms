lc = 200;

Point(1) = { -40000,-40000, 0.0, lc};
Point(2) = {  40000,-40000, 0.0, lc};
Point(3) = {  40000,     0, 0.0, lc};
Point(4) = { -40000,     0, 0.0, lc};
Point(5) = { -40000, 40000, 0.0, lc};
Point(6) = {  40000, 40000, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {3,6};
Line(6) = {6,5};
Line(7) = {5,4};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1}; 

Line Loop(3) = {5,6,7,-3};
Plane Surface(4) = {3};

// Physical Curve("10")  = 3;
// Physical Curve("20")  = 7;
// Physical Surface("2") = 2;
// Physical Surface("4") = 4;

Mesh.Algorithm = 5;