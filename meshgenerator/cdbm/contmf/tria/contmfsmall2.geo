// Gmsh project created on Sun Nov  6 10:46:12 2022
// Step-over embedded fault mesh 

lc = 400;

Point(1) = {     4000,   -1000, 0.0, lc};
Point(2) = {     4000,    1000, 0.0, lc};
Point(3) = {    -2500,    1000, 0.0, lc};
Point(4) = {    -2500,   -1000, 0.0, lc};
Point(5) = {    -1500,       0, 0.0, lc};
Point(6) = {     3000,       0, 0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {5,6};

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Line{5} In Surface{1};

Field[1] = Box;
Field[1].VIn = lc / 16;
Field[1].VOut = lc;
Field[1].XMin =  -8000;
Field[1].XMax =  11000;
Field[1].YMin =  -2000;
Field[1].YMax =   1000;
Field[1].Thickness = lc;

Background Field = 1;