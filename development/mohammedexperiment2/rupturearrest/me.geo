//use m instead of mm

//mesh size //make sure to resolve frictional length by 8 elements minimum
lc = 0.0009;

//outer domain points
Point(1) = {       0,       0, 0, lc };
Point(2) = {   0.203,       0, 0, lc };
Point(3) = {   0.203,   0.250, 0, lc };
Point(4) = {       0,   0.250, 0, lc };
Point(5) = {       0,  0.0375, 0, lc };
Point(6) = {   0.203,   0.150, 0, lc };

//inner domain points
Point(7) = { 0.142127, 0.116257, 0, lc };
Point(8) = { 0.141427, 0.115869, 0, lc };
Point(9) = { 0.141427, 0.116669, 0, lc };

//crack end point
Point(10) = { 0.112338, 0.169146, 0, lc };

//outer domain lines
Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,1};

//embeded line 1
Line(7) = {5,8};

//embeded line 2
Line(8) = {6,7};

//inner domain lines
Line(9) = {7,8};
Line(10) = {8,9};
Line(11) = {9,7};

//embeded line 3
Line(12) = {9,10};

//inner loop & surface
Line Loop(1) = {9,10,11};
Plane Surface(2) = {1};

//outer loop & surface with holes
Line Loop(3) = {1,2,3,4,5,6};
Plane Surface(4) = {3,1}; //{outer loop, inner loop}

//embed line in surface
Line{7,8,12} In Surface{4};

//define physical curves
Physical Curve("embeded1") = 7;
Physical Curve("embeded2") = 8;
Physical Curve("embeded3") = 12;

//define physical surface
Physical Surface("psurf") = 4;

Mesh.Algorithm = 5;