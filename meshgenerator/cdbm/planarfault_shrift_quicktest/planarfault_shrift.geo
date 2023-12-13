//25m mesh generation
lc = 25;

Point(1) = {    2000,   -2000, 0, lc };
Point(2) = {    2000,    2000, 0, lc };
Point(3) = {   -2000,    2000, 0, lc };
Point(4) = {   -2000,   -2000, 0, lc };
Point(5) = {   -2000,    -100, 0, lc };
Point(6) = {    2000,    -100, 0, lc };

Line(1) = {1,6};
Line(2) = {6,2};
Line(3) = {2,3};
Line(4) = {3,5};
Line(5) = {5,4};
Line(6) = {4,1};
Line(7) = {5,6};

Line Loop(1) = {1,-7,5,6};
Plane Surface(2) = {1};

Line Loop(3) = {2,3,4,7};
Plane Surface(4) = {3};

Mesh.Algorithm = 5;