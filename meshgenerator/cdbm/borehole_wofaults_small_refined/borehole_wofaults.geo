//10m mesh generation
lc = 10;

//outer boundary
Point(1) = {   -1000,    -1000, 0, lc };
Point(2) = {    1000,    -1000, 0, lc };
Point(3) = {    1000,     1000, 0, lc };
Point(4) = {   -1000,     1000, 0, lc };

//circle hole
Point(5) = {       0,       0, 0, lc };
Point(6) = {    -200,       0, 0, lc };
Point(7) = {     200,       0, 0, lc };

//lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//circle
Circle(5) = {6,5,7};
Circle(6) = {7,5,6};

Curve Loop(1) = {5,6};

Curve Loop(3) = {1,2,3,4};
Plane Surface(4) = {3,1};

Mesh.Algorithm = 5;

Physical Curve("bottom") = {1};
Physical Curve("borehole") = {5,6};
Physical Surface("domain") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("left") = {4};
