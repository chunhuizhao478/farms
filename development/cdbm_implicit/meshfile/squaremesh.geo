//25m mesh generation
lc = 100;

Point(1) = {   5000,  -5000, 0, lc };
Point(2) = {   5000,   5000, 0, lc };
Point(3) = {  -5000,   5000, 0, lc };
Point(4) = {  -5000,  -5000, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

Mesh.Algorithm = 5;
