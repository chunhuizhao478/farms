//0.02m mesh generation
lc = 5;

Point(1) = {   -50.0,  -50, 0, lc };
Point(2) = {    50.0,  -50, 0, lc };
Point(3) = {    50.0, -0.1, 0, lc };
Point(4) = {   -50.0, -0.1, 0, lc };
Point(5) = {    50.0,  0.1, 0, lc };
Point(6) = {    50.0,   50, 0, lc };
Point(7) = {   -50.0,   50, 0, lc };
Point(8) = {   -50.0,  0.1, 0, lc };

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {3,5};
Line(6) = {5,8};
Line(7) = {8,4};
Line(8) = {5,6};
Line(9) = {6,7};
Line(10) = {7,8};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

Line Loop(3) = {5,6,7,-3};
Plane Surface(4) = {3};

Line Loop(5) = {8,9,10,-6};
Plane Surface(6) = {5};

Field[1] = Box;
Field[1].VIn = lc / 100;
Field[1].VOut = lc;
Field[1].XMin =  -200;
Field[1].XMax =   200;
Field[1].YMin =  -0.2;
Field[1].YMax =   0.2;
Field[1].Thickness = lc;

Field[2] = MathEval;
Field[2].F = "0.1 * y ^ 2 + 0.101";

Field[3] = Min;
Field[3].FieldsList = {1,2};
Background Field = 3;

Mesh.Algorithm = 5;