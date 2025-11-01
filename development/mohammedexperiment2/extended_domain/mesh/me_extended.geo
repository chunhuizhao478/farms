//use m instead of mm

//mesh size //make sure to resolve frictional length by 8 elements minimum
lc_fault = 0.0009;
lc_far = 0.005;

//outer domain points
Point(1) = {       0,       0, 0, lc_far };
Point(2) = {   0.384,       0, 0, lc_far };
Point(3) = {   0.384,   0.300, 0, lc_far };
Point(4) = {       0,   0.300, 0, lc_far };
Point(5) = {       0,  0.036688, 0, lc_far };
Point(6) = {   0.384,   0.250, 0, lc_far };

//inner domain points
Point(7) = { 0.263167, 0.182564, 0, lc_fault };
Point(8) = { 0.262380, 0.182128, 0, lc_fault };
Point(9) = { 0.262380, 0.183028, 0, lc_fault };

//crack end point
Point(10) = { 0.233292, 0.235504, 0, lc_fault };

// Point at distance 0.010 from P9: (0.136579, 0.125415, 0.000000)
// Point at distance 0.015 from P9: (0.134155, 0.129788, 0.000000)

//outer domain lines
Line(1) = {1,2};
Line(2) = {2,6};
Line(3) = {6,3};
Line(4) = {3,4};
Line(5) = {4,5};
Line(6) = {5,1};

//embeded line 1
Line(7) = {5,8};

//embeded line 2
Line(8) = {6,7};

//inner domain lines
Line(9) = {7,8};
Line(10) = {8,9};
Line(11) = {9,7};

//embeded line 3
Line(12) = {9,10};

//inner loop (defines a hole in the outer surface)
Line Loop(1) = {9,10,11};

//outer loop & surface with inner loop excluded
Line Loop(3) = {1,2,3,4,5,6};
Plane Surface(4) = {3,-1}; //{outer loop, inner loop as hole}

//embed line in surface
Line{7,8,12} In Surface{4};

//adaptive mesh refinement near fault lines
Field[1] = Distance;
Field[1].CurvesList = {7,8,9,10,11,12};
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_fault;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.003;
Field[2].DistMax = 0.02;

Background Field = 2;

Mesh.CharacteristicLengthMin = lc_fault;
Mesh.CharacteristicLengthMax = lc_far;

//define physical curves
Physical Curve("embeded1") = 7;
Physical Curve("embeded2") = 8;
Physical Curve("embeded3") = 12;

//define physical surface
Physical Surface("psurf") = 4;

Mesh.Algorithm = 5;
