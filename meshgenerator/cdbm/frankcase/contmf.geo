//25m mesh generation
lc = 20;

Point(1) = {   -1000,    -500, 0, lc };
Point(2) = {    1000,    -500, 0, lc };
Point(3) = {    1000,     500, 0, lc };
Point(4) = {   -1000,     500, 0, lc };

Point(5) = {    -600,      40, 0, lc };
Point(6) = {    -200,      40, 0, lc };

Point(7) = {     200,     -40, 0, lc };
Point(8) = {     600,     -40, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,6};
Line(6) = {7,8};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

Field[1] = Box;
Field[1].VIn = lc / 2;
Field[1].VOut = lc;
Field[1].XMin =  -1000;
Field[1].XMax =   1000;
Field[1].YMin =  -500;
Field[1].YMax =   500;
Field[1].Thickness = lc;

Background Field = 1;

Line{5} In Surface{2};
Line{6} In Surface{2};

Mesh.Algorithm = 5;