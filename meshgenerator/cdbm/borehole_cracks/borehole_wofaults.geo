//25m mesh generation
lc = 25;

//outer boundary
Point(1) = {   -2000,    -2000, 0, lc };
Point(2) = {    2000,    -2000, 0, lc };
Point(3) = {    2000,     2000, 0, lc };
Point(4) = {   -2000,     2000, 0, lc };

//circle hole
Point(5) = {       0,       0, 0, lc };

//additional points for lines
Point(6)  = {          200,         0, 0, lc}; //1
Point(7)  = {          800,         0, 0, lc};
Point(8)  = {    173.20508,       100, 0, lc}; //2 *
Point(9)  = {    692.82032,       400, 0, lc};
Point(10) = {          100, 173.20508, 0, lc}; //3 *
Point(11) = {          400, 692.82032, 0, lc};
Point(12) = {            0,       200, 0, lc}; //4
Point(13) = {            0,       800, 0, lc};
Point(14) = {         -100, 173.20508, 0, lc}; //5 *
Point(15) = {         -400, 692.82032, 0, lc};
Point(16) = {   -173.20508,       100, 0, lc}; //6 *
Point(17) = {   -692.82032,       400, 0, lc};
Point(18) = {         -200,         0, 0, lc}; //7
Point(19) = {         -800,         0, 0, lc};
Point(20) = {   -173.20508,      -100, 0, lc}; //8 *
Point(21) = {   -692.82032,      -400, 0, lc};
Point(22) = {         -100,-173.20508, 0, lc}; //9 *
Point(23) = {         -400,-692.82032, 0, lc};
Point(24) = {            0,      -200, 0, lc}; //10
Point(25) = {            0,      -800, 0, lc};
Point(26) = {          100,-173.20508, 0, lc}; //11 *
Point(27) = {          400,-692.82032, 0, lc};
Point(28) = {    173.20508,      -100, 0, lc}; //12 *
Point(29) = {    692.82032,      -400, 0, lc};

//lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//+
Circle(5) = {6, 5, 8};
//+
Circle(6) = {8, 5, 10};
//+
Circle(7) = {10, 5, 12};
//+
Circle(8) = {12, 5, 14};
//+
Circle(9) = {14, 5, 16};
//+
Circle(10) = {16, 5, 18};
//+
Circle(11) = {18, 5, 20};
//+
Circle(12) = {20, 5, 22};
//+
Circle(13) = {22, 5, 24};
//+
Circle(14) = {24, 5, 26};
//+
Circle(15) = {26, 5, 28};
//+
Circle(16) = {28, 5, 6};

Curve Loop(1) = {5,6,7,8,9,10,11,12,13,14,15,16};

//square domain
Curve Loop(3) = {1,2,3,4};
Plane Surface(4) = {3,1};

//+
Line(17) = {6, 7};
//+
Line(18) = {8, 9};
//+
Line(19) = {10, 11};
//+
Line(20) = {12, 13};
//+
Line(21) = {14, 15};
//+
Line(22) = {16, 17};
//+
Line(23) = {18, 19};
//+
Line(24) = {20, 21};
//+
Line(25) = {22, 23};
//+
Line(26) = {24, 25};
//+
Line(27) = {26, 27};
//+
Line(28) = {28, 29};

Line{17} In Surface{4};
Line{18} In Surface{4};
Line{19} In Surface{4};
Line{20} In Surface{4};
Line{21} In Surface{4};
Line{22} In Surface{4};
Line{23} In Surface{4};
Line{24} In Surface{4};
Line{25} In Surface{4};
Line{26} In Surface{4};
Line{27} In Surface{4};
Line{28} In Surface{4};

Mesh.Algorithm = 5;

Physical Curve("bottom") = {1};
Physical Curve("borehole") = {5,6,7,8,9,10,11,12,13,14,15,16};
Physical Surface("domain") = {4};
Physical Curve("top") = {3};
Physical Curve("right") = {2};
Physical Curve("left") = {4}; //+
Physical Curve("embeded1") = {17};
Physical Curve("embeded2") = {18};
Physical Curve("embeded3") = {19};
Physical Curve("embeded4") = {20};
Physical Curve("embeded5") = {21};
Physical Curve("embeded6") = {22};
Physical Curve("embeded7") = {23};
Physical Curve("embeded8") = {24};
Physical Curve("embeded9") = {25};
Physical Curve("embeded10") = {26};
Physical Curve("embeded11") = {27};
Physical Curve("embeded12") = {28};

