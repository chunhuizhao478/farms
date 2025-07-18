// 2D ring with two annuli (no mesh inside the inner hole)
SetFactory("OpenCASCADE");

// Characteristic lengths
lc0 = 0.0001; // super-fine near hole
lc = 0.0001;  // fine transition
lc2 = 0.005;  // coarse outer

radius_outer = 0.125;
radius_refined = 0.0205;
radius_inner = 0.020;

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

//---- 6) CORRECTED Mesh‐size fields ----
Field[1] = Distance;
Field[1].NodesList = {13};

// Field[2]: Outer region grading (radius_refined → radius_outer)
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc;           // size at radius_refined
Field[2].LcMax = lc2;          // size at radius_outer
Field[2].DistMin = radius_refined;
Field[2].DistMax = radius_outer;

// Field[3]: Inner region sizing (radius_inner → radius_refined)
// CORRECTED: LcMax should be lc, not lc2
Field[3] = Threshold;
Field[3].IField = 1;
Field[3].LcMin = lc0;          // size at radius_inner
Field[3].LcMax = lc2;           // size at radius_refined (CORRECTED)
Field[3].DistMin = radius_inner;
Field[3].DistMax = radius_refined;

// Field[4]: Take minimum of both fields
Field[4] = Min;
Field[4].FieldsList = {2,3};

Background Field = 4;

//---- 7) Physical groups ----
Physical Curve("OuterBoundary") = {1,2,3,4};
Physical Curve("RefinedBoundary")= {5,6,7,8};
Physical Curve("HoleBoundary") = {9,10,11,12};
Physical Surface("InnerBlock") = {200};
Physical Surface("OuterBlock") = {201};