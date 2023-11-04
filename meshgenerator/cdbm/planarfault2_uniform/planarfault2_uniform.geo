//25m mesh generation
lc = 25;

Point(1) = {   10000,   -3000, 0, lc };
Point(2) = {   10000,    3000, 0, lc };
Point(3) = {  -10000,    3000, 0, lc };
Point(4) = {  -10000,   -3000, 0, lc };
Point(5) = {  -10000,       0, 0, lc };
Point(6) = {   10000,       0, 0, lc };

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

Field[1] = Box;
Field[1].VIn = lc;
Field[1].VOut = lc;
Field[1].XMin =  -20000;
Field[1].XMax =   20000;
Field[1].YMin =  -500;
Field[1].YMax =   500;
Field[1].Thickness = lc;

Background Field = 1;

Mesh.Algorithm = 5;