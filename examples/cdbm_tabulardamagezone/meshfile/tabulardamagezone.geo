//25m mesh generation
lc = 0.01;

Point(1) = {   -2.0,  -0.4, 0, lc };
Point(2) = {    2.0,  -0.4, 0, lc };
Point(3) = {    2.0,   0.4, 0, lc };
Point(4) = {   -2.0,   0.4, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

Mesh.Algorithm = 5;