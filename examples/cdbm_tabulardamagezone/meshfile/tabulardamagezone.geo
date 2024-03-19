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

// Field[1] = Box;
// Field[1].VIn = lc / 5;
// Field[1].VOut = lc;
// Field[1].XMin =  -1.21;
// Field[1].XMax =   1.21;
// Field[1].YMin =  -0.21;
// Field[1].YMax =   0.21;
// Field[1].Thickness = lc;

// Background Field = 1;

Mesh.Algorithm = 5;