// Define parameters
Lx = 8000;
Ly = 6000;
Lz = 6000;

lc1 = 100;
lc2 = 1000;

//Define points
Point(1) = {-Lx, -Ly, -Lz, lc2};
Point(2) = {+Lx, -Ly, -Lz, lc2};
Point(3) = {+Lx, -Ly,   0, lc2};
Point(4) = {-Lx, -Ly,   0, lc2};

Point(5) = {-Lx,   0, -Lz, lc1};
Point(6) = {+Lx,   0, -Lz, lc1};
Point(7) = {+Lx,   0,   0, lc1};
Point(8) = {-Lx,   0,   0, lc1};

Point(9)  = {-Lx, +Ly, -Lz, lc2};
Point(10) = {+Lx, +Ly, -Lz, lc2};
Point(11) = {+Lx, +Ly,   0, lc2};
Point(12) = {-Lx, +Ly,   0, lc2};

// Define lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {1, 5};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};

Line(9) = {5, 6};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

Line(13) = {5, 9};
Line(14) = {6, 10};
Line(15) = {7, 11};
Line(16) = {8, 12};

Line(17) = {9, 10};
Line(18) = {10, 11};
Line(19) = {11, 12};
Line(20) = {12, 9};

// Define surfaces
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {1, 6, -9, -5};
Plane Surface(2) = {2};

Line Loop(3) = {6, 10, -7, -2};
Plane Surface(3) = {3};

Line Loop(4) = {7, 11, -8, -3};
Plane Surface(4) = {4};

Line Loop(5) = {8, 12, -5, -4};
Plane Surface(5) = {5};

Line Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};

Line Loop(7) = {9, 14, -17, -13};
Plane Surface(7) = {7};

Line Loop(8) = {14, 18, -15, -10};
Plane Surface(8) = {8};

Line Loop(9) = {15, 19, -16, -11};
Plane Surface(9) = {9};

Line Loop(10) = {16, 20, -13, -12};
Plane Surface(10) = {10};

Line Loop(11) = {17, 18, 19, 20};
Plane Surface(11) = {11};

// Define volumes
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Surface Loop(2) = {6, 7, 8, 9, 10, 11};
Volume(2) = {2};

//
Field[1] = Box;
Field[1].VIn = lc1;
Field[1].VOut = lc2;
Field[1].XMin =  -6000;
Field[1].XMax =   6000;
Field[1].YMin =  -1500;
Field[1].YMax =   1500;
Field[1].ZMin =  -6000;
Field[1].ZMax =   0;
Field[1].Thickness = Ly;

Background Field = 1;

