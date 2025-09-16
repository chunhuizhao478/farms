//2D algorithm: 11. Quasi-Structed Quad (Mesh.Algorithm)
//2D recombination algorithm: 2. Simple Full-Quad (Mesh.RecombinationAlgorithm)

lc = 100;

Point(1) = {  -20000, -20000, 0.0, lc};
Point(2) = {   30000, -20000, 0.0, lc};
Point(3) = {   30000,  20000, 0.0, lc};
Point(4) = {  -20000,  20000, 0.0, lc};

Point(5) = {     -10000,    0.0,  0.0, lc};
Point(6) = {      18000,    0.0,  0.0, lc};

Point(7)  = {      6000,    0.0,  0.0, lc}; //along main fault
Point(8)  = {   16392.3,  -6000,  0.0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line(5) = {5,7};
Line(6) = {7,6};
Line(7) = {7,8};

Line Loop(1) = {1,2,3,4};
Plane Surface(2) = {1};

// Transfinite Curve {5} = 161;  //28000
// Transfinite Curve {6} = 121;  //12000
// Transfinite Curve {7} = 121;

Point{5,6,7,8} In Surface{2};
Line{5,6,7} In Surface{2};

Physical Curve("180") = {5,6};
Physical Curve("30")  = 7;

Physical Curve("bottom")  = 1;
Physical Curve("right") = 2;
Physical Curve("top")  = 3;
Physical Curve("left")  = 4;

Physical Surface("100") = 2;

Mesh.Algorithm = 6;
Mesh.RecombinationAlgorithm = 2;

//v1 v2
// Mesh.Algorithm = 11;
// Mesh.RecombinationAlgorithm = 2;

//v3
//Mesh.Algorithm = 6;

//v4
//Mesh.Algorithm = 6; enlarge domain