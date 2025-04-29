// Parameters
lc = 5e-5;
lc_refined = 5e-5;
lc_refined_2 = 5e-5;

// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {0.001, 0, 0, lc};
Point(3) = {0.001, 0.001, 0, lc};
Point(4) = {0, 0.001, 0, lc};
Point(5) = {0, 0.000501, 0, lc};
Point(6) = {0.0005, 0.000500, 0, lc};
Point(7) = {0, 0.000499, 0, lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};

// Line loop and surface
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

// Refine mesh near the notch
Field[1] = Box;
Field[1].VIn = lc_refined;  // Mesh size inside the fault zone
Field[1].VOut = lc;       // Mesh size outside the fault zone
Field[1].XMin = 0.00045;     // Adjusted to cover the notch area
Field[1].XMax = 0.001;     // Adjusted to cover the notch area
Field[1].YMin = 0.00048;    // Adjusted to cover the notch area
Field[1].YMax = 0.00052;    // Adjusted to cover the notch area
Field[1].Thickness = 0.001; // Transition thickness

// Refine mesh near the notch
Field[2] = Box;
Field[2].VIn = lc_refined_2;  // Mesh size inside the fault zone
Field[2].VOut = lc;       // Mesh size outside the fault zone
Field[2].XMin = 0.00048;     // Adjusted to cover the notch area
Field[2].XMax = 0.00052;     // Adjusted to cover the notch area
Field[2].YMin = 0.00048;    // Adjusted to cover the notch area
Field[2].YMax = 0.00052;    // Adjusted to cover the notch area
Field[2].Thickness = 0.001; // Transition thickness

// Use the minimum of the two fields
Field[3] = Min;
Field[3].FieldsList = {1, 2};

// Set the background mesh field
Background Field = 3;

Mesh 2; // Generate 2D mesh

// Adding notch line to the surface
Physical Line("Left") = {4,7};
Physical Line("Right") = {2};
Physical Line("Top") = {3};
Physical Line("Bottom") = {1};
Physical Surface("Domain") = {1};
