// 2D ring with two annuli (no mesh inside the inner hole)
SetFactory("OpenCASCADE");

// Characteristic lengths
lc0 = 0.0001; // super-fine near hole
lc = 0.0001;  // fine transition
lc_strip = 0.0001; // fine mesh in crack propagation strip
lc2 = 0.01;  // coarse outer

radius_outer = 0.2;
radius_refined = 0.0925;
radius_inner = 0.09;

//---- Center point (for arc definitions + distance field) ----
Point(13) = {0.0, 0.0, 0.0, lc0};

//---- 1) Outer circle points & arcs ----
Point(1) = { radius_outer, 0.0, 0.0, lc2 };
Point(2) = { 0.0, radius_outer, 0.0, lc2 };
Point(3) = {-radius_outer, 0.0, 0.0, lc2 };
Point(4) = { 0.0, -radius_outer, 0.0, lc2 };

Circle(1) = {1, 13, 2};
Circle(2) = {2, 13, 3};
Circle(3) = {3, 13, 4};
Circle(4) = {4, 13, 1};

//---- 2) Refined‐zone circle (for splitting region) ----
Point(5) = { radius_refined, 0.0, 0.0, lc };
Point(6) = { 0.0, radius_refined, 0.0, lc };
Point(7) = {-radius_refined, 0.0, 0.0, lc };
Point(8) = { 0.0, -radius_refined, 0.0, lc };

Circle(5) = {5, 13, 6};
Circle(6) = {6, 13, 7};
Circle(7) = {7, 13, 8};
Circle(8) = {8, 13, 5};

//---- 3) Inner‐hole boundary (no mesh inside) ----
Point(9) = { radius_inner, 0.0, 0.0, lc0 };
Point(10) = { 0.0, radius_inner, 0.0, lc0 };
Point(11) = {-radius_inner, 0.0, 0.0, lc0 };
Point(12) = { 0.0, -radius_inner, 0.0, lc0 };

Circle(9) = {9, 13, 10};
Circle(10) = {10, 13, 11};
Circle(11) = {11, 13, 12};
Circle(12) = {12, 13, 9};

//---- 4) Line loops ----
Line Loop(100) = {1,2,3,4};   // outer loop
Line Loop(101) = {5,6,7,8};   // refined loop
Line Loop(102) = {9,10,11,12}; // inner-hole loop

//---- 5) Two annular surfaces ----
// 5.1 Inner annulus: radius_inner → radius_refined
Plane Surface(200) = {101, 102};
// 5.2 Outer annulus: radius_refined → radius_outer
Plane Surface(201) = {100, 101};

//---- 6) Mesh‐size fields for two box regions ----

// Field[1]: Horizontal box for horizontal crack propagation
Field[1] = Box;
Field[1].VIn = lc_strip;       // uniform fine mesh inside horizontal box
Field[1].VOut = lc2;           // coarse mesh outside
Field[1].XMin = 0;
Field[1].XMax = radius_outer;
Field[1].YMin = -4 * 4e-4;
Field[1].YMax =  4 * 4e-4;
Field[1].Thickness = 0.02;    // 5mm transition zone

// Field[2]: Vertical box for vertical crack propagation
//Field[2] = Box;
//Field[2].VIn = lc_strip;       // uniform fine mesh inside vertical box
//Field[2].VOut = lc2;           // coarse mesh outside
//Field[2].XMin = -4 * 4e-4;
//Field[2].XMax =  4 * 4e-4;
//Field[2].YMin = -radius_outer;
//Field[2].YMax = radius_outer;
//Field[2].Thickness = 0.02;    // 5mm transition zone

// Field[3]: Near hole refinement (distance-based)
Field[3] = Distance;
Field[3].NodesList = {13};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc0;          // fine near hole
Field[4].LcMax = lc2;          // coarse far from hole
Field[4].DistMin = radius_inner;
Field[4].DistMax = radius_refined;

// Field[5]: Take minimum of all fields (boxes take priority where finer)
Field[5] = Min;
Field[5].FieldsList = {1, 4};

Background Field = 5;

//---- 7) Physical groups ----
Physical Curve("OuterBoundary") = {1,2,3,4};
Physical Curve("RefinedBoundary")= {5,6,7,8};
Physical Curve("HoleBoundary") = {9,10,11,12};
Physical Surface("InnerBlock") = {200};
Physical Surface("OuterBlock") = {201};
