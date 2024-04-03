//0.02m mesh generation
lc = 0.2;

Point(1) = {   -20.0,  -20, 0, lc };
Point(2) = {    20.0,  -20, 0, lc };
Point(3) = {    20.0,   20, 0, lc };
Point(4) = {   -20.0,   20, 0, lc };

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