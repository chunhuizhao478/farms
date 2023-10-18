//25m mesh generation
lc = 200;

Point(1) = {  -10000,   -5000, 0, lc };
Point(2) = {   10000,   -5000, 0, lc };
Point(3) = {   10000,    5000, 0, lc };
Point(4) = {  -10000,    5000, 0, lc };
Point(5) = {  -5000,        0, 0, lc };
Point(6) = {   5000,        0, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

Line{5} In Surface{2};

Field[1] = Box;
Field[1].VIn = lc / 8;
Field[1].VOut = lc;
Field[1].XMin =  -7000;
Field[1].XMax =   7000;
Field[1].YMin =  -2000;
Field[1].YMax =   2000;
Field[1].Thickness = lc;

Background Field = 1;

Mesh.Algorithm = 5;