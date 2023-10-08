//25m mesh generation
lc = 100;

Point(1) = {   1000,  -625, 0, lc };
Point(2) = {   1000,     0, 0, lc };
Point(3) = {   1000,   625, 0, lc };
Point(4) = {  -4000,   625, 0, lc };
Point(5) = {  -4000,      0, 0, lc };
Point(6) = {  -4000,  -625, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};

Line Loop(1) = {1,7,5,6};
Plane Surface(2) = {1};

Line Loop(3) = {2,3,4,-7};
Plane Surface(4) = {3};

Field[1] = Box;
Field[1].VIn = lc / 8;
Field[1].VOut = lc;
Field[1].XMin =  -4000;
Field[1].XMax =   4000;
Field[1].YMin =   -300;
Field[1].YMax =    300;
Field[1].Thickness = lc;

Background Field = 1;

Mesh.Algorithm = 5;