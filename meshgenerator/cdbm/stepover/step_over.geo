// Gmsh project created on Sun Nov  6 10:46:12 2022
// Step-over embedded fault mesh 

lc = 200;

Point(1) = {     15000,   -6000, 0.0, lc};
Point(2) = {     15000,    5000, 0.0, lc};
Point(3) = {    -12000,    5000, 0.0, lc};
Point(4) = {    -12000,   -6000, 0.0, lc};
Point(5) = {    -5000,       0, 0.0, lc};
Point(6) = {     3000,       0, 0.0, lc};
Point(7) = {    -1000,   -1000, 0.0, lc};
Point(8) = {     8000,   -1000, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};
Line(6) = {7,8};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Line{5} In Surface{1};
Line{6} In Surface{1};

Field[1] = Box;
Field[1].VIn = lc / 8;
Field[1].VOut = lc;
Field[1].XMin =  -8000;
Field[1].XMax =  11000;
Field[1].YMin =  -3000;
Field[1].YMax =   2000;
Field[1].Thickness = lc;

Background Field = 1;