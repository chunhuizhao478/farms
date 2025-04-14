lc = 200;
lc2 = 1000;

// Define box dimensions
xmin = -15000;
xmax = 15000;
ymin = -15000;
ymax = 0;
zmin = -15000;
zmax = 15000;
zmid = 0; // This will be the separation plane for two boxes

zuniform = 2000;

// Points for the entire volume with a separating middle surface
Point(1) = {xmin, ymin, zmin, lc2};
Point(2) = {xmax, ymin, zmin, lc2};
Point(3) = {xmax, ymax, zmin, lc2};
Point(4) = {xmin, ymax, zmin, lc2};
Point(5) = {xmin, ymin, zmid, lc};
Point(6) = {xmax, ymin, zmid, lc};
Point(7) = {xmax, ymax, zmid, lc};
Point(8) = {xmin, ymax, zmid, lc};
Point(9) = {xmin, ymin, zmax, lc2};
Point(10) = {xmax, ymin, zmax, lc2};
Point(11) = {xmax, ymax, zmax, lc2};
Point(12) = {xmin, ymax, zmax, lc2};

// Lines for the entire volume
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

// Surfaces for the entire volume with a separating middle surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};  // Bottom surface
Line Loop(2) = {5, 9, -6, -1};
Plane Surface(2) = {2};  // Side surface
Line Loop(3) = {6, 10, -7, -2};
Plane Surface(3) = {3};  // Side surface
Line Loop(4) = {7, 11, -8, -3};
Plane Surface(4) = {4};  // Side surface
Line Loop(5) = {8, 12, -5, -4};
Plane Surface(5) = {5};  // Side surface
Line Loop(6) = {9, 10, 11, 12};
Plane Surface(6) = {6};  // Middle separating surface
Line Loop(7) = {13, 17, -14, -9};
Plane Surface(7) = {7};  // Side surface
Line Loop(8) = {14, 18, -15, -10};
Plane Surface(8) = {8};  // Side surface
Line Loop(9) = {15, 19, -16, -11};
Plane Surface(9) = {9};  // Side surface
Line Loop(10) = {16, 20, -13, -12};
Plane Surface(10) = {10};  // Side surface
Line Loop(11) = {17, 18, 19, 20};
Plane Surface(11) = {11};  // Top surface

// Define the volume
Surface Loop(1) = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11};  // Include all surfaces except the middle surface
Volume(1) = {1};  // Single volume containing the middle surface

// Define mesh field to control mesh size along the normal direction (Z-axis)
Field[1] = Box;
Field[1].VIn = lc;
Field[1].VOut = lc2;
Field[1].XMin = xmin;
Field[1].XMax = xmax;
Field[1].YMin = ymin;
Field[1].YMax = ymax;
Field[1].ZMin = -zuniform;
Field[1].ZMax = zuniform;

// Apply the mesh field
Background Field = 1;

// Assign Physical Groups (optional)
// Physical Surface("BottomSurface") = {1};
// Physical Surface("TopSurface") = {7};
Physical Volume("Box") = {1};