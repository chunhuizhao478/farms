//25m mesh generation
lc = 25;
lc2 = 125;

Point(1) = {   10000,   -10000, 0, lc2 };
Point(2) = {   10000,    10000, 0, lc2 };
Point(3) = {  -10000,    10000, 0, lc2 };
Point(4) = {  -10000,   -10000, 0, lc2 };
Point(5) = {  -10000,        0, 0, lc };
Point(6) = {   10000,        0, 0, lc };

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

Field[1] = Distance;
Field[1].PointsList = {5,6};
Field[1].CurvesList = {7};
Field[1].Sampling = 1000;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc;
Field[2].SizeMax = lc2;
Field[2].DistMin = 3000;
Field[2].DistMax = 10000;

Field[7] = Min;
Field[7].FieldsList = {2};
Background Field = 7;

Mesh.Algorithm = 5;