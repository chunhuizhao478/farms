// Parameters
lc = 5e-3;
lc_refined = 1e-4;
lc_refined_2 = 1e-4;

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {0.1, 0, 0, lc};
Point(3) = {0.1, 0.04, 0, lc};
Point(4) = {0, 0.04, 0, lc};
Point(5) = {0, 0.02005, 0, lc};
Point(6) = {0.05, 0.02005, 0, lc_refined};
Point(7) = {0.05, 0.01995, 0, lc_refined};
Point(8) = {0, 0.01995, 0, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Line loop and surface
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Refine mesh near the notch
Field[1] = Box;
Field[1].VIn = lc_refined;  // Mesh size inside the fault zone
Field[1].VOut = lc;       // Mesh size outside the fault zone
Field[1].XMin = 0.0475;     // Adjusted to cover the notch area
Field[1].XMax = 0.0575;     // Adjusted to cover the notch area
Field[1].YMin = 0.0175;    // Adjusted to cover the notch area
Field[1].YMax = 0.0225;    // Adjusted to cover the notch area
Field[1].Thickness = 0.1; // Transition thickness

Background Field =1;

// Mesh 2; // Generate 2D mesh

// Adding notch line to the surface
Physical Line("Left") = {4,8};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Bottom") = {1};
Physical Surface("Domain") = {1};