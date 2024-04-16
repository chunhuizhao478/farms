//0.02m mesh generation
lc = 2;

Point(1)  = {   -10.0,  -50, 0, lc };
Point(2)  = {    10.0,  -50, 0, lc };
Point(3)  = {   -10.0,  -45, 0, lc };
Point(4)  = {    10.0,  -45, 0, lc };
Point(5)  = {   -10.0, -0.1, 0, lc/100 };
Point(6)  = {    10.0, -0.1, 0, lc/100 };
Point(7)  = {   -10.0,  0.1, 0, lc/100 };
Point(8)  = {    10.0,  0.1, 0, lc/100 };
Point(9)  = {   -10.0,   45, 0, lc };
Point(10) = {    10.0,   45, 0, lc };
Point(11) = {   -10.0,   50, 0, lc };
Point(12) = {    10.0,   50, 0, lc };

Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {3,4};
Line(4) = {3,1};
Line(5) = {4,6};
Line(6) = {5,6};
Line(7) = {5,3};
Line(8) = {6,8};
Line(9) = {7,8};
Line(10) = {7,5};
Line(11) = {8,10};
Line(12) = {9,10};
Line(13) = {9,7};
Line(14) = {10,12};
Line(15) = {11,12};
Line(16) = {11,9};

Line Loop(1) = {1,2,-3,4};
Plane Surface(2) = {1};

Line Loop(3) = {3,5,-6,7};
Plane Surface(4) = {3};

Line Loop(5) = {6,8,-9,10};
Plane Surface(6) = {5};

Line Loop(7) = {9,11,-12,13};
Plane Surface(8) = {7};

Line Loop(9) = {12,14,-15,16};
Plane Surface(10) = {9};

Field[1] = Box;
Field[1].VIn  = lc / 100;
Field[1].VOut = lc * 4;
Field[1].XMin =  -50;
Field[1].XMax =   50;
Field[1].YMin =  -0.2;
Field[1].YMax =   0.2;
Field[1].Thickness = lc / 2;

Field[2] = MathEval;
Field[2].F = "0.1 * y ^ 2 + 0.101";

Field[3] = Min;
Field[3].FieldsList = {1,2};
Background Field = 3;

Mesh.Algorithm = 5;