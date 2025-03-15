SetFactory("OpenCASCADE");

// Define mesh sizes
lc_fault = 100;    // Fine mesh size inside the small boxes
lc = 2500;        // Coarse mesh size outside the small boxes

// Define the big square (2D) - full domain
big_xmin = -30000;
big_xmax = 30000;
big_ymin = -30000;
big_ymax = 30000;

// Define the small boxes (2D) - areas for refined mesh
small_xmin = -7500;
small_xmax = 7500;
small_ymin = -4000;
small_ymax = 4000;

// Define points for the big square
Point(1) = {big_xmin, big_ymin, 0, lc};
Point(2) = {big_xmax, big_ymin, 0, lc};
Point(3) = {big_xmax, big_ymax, 0, lc};
Point(4) = {big_xmin, big_ymax, 0, lc};

// Define lines for the big square
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create the outer surface
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define points for the lower small box (below y=0)
Point(5) = {small_xmin, small_ymin, 0, lc_fault};
Point(6) = {small_xmax, small_ymin, 0, lc_fault};
Point(7) = {small_xmax, 0, 0, lc_fault};
Point(8) = {small_xmin, 0, 0, lc_fault};

// Define lines for the lower small box
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Create line loop and surface for the lower small box
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// Define points for the upper small box (above y=0)
Point(9) = {small_xmin, 0, 0, lc_fault};
Point(10) = {small_xmax, 0, 0, lc_fault};
Point(11) = {small_xmax, small_ymax, 0, lc_fault};
Point(12) = {small_xmin, small_ymax, 0, lc_fault};

// Define lines for the upper small box
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};

// Create line loop and surface for the upper small box
Line Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3};

// Perform BooleanFragments to fragment the surfaces
BooleanFragments{ Surface{1}; Delete; }{ Surface{2, 3}; Delete; }

// Collect the new surfaces
surfaces[] = Surface{:};

// Define mesh size field: refined mesh only inside the small boxes
Field[1] = Box;
Field[1].VIn = lc_fault;    // Fine mesh size inside the small boxes
Field[1].VOut = lc;         // Coarse mesh size outside
Field[1].XMin = small_xmin;
Field[1].XMax = small_xmax;
Field[1].YMin = small_ymin;
Field[1].YMax = small_ymax;

// Apply the mesh size field
Background Field = 1;

// Assign physical groups for each box and the big domain
Physical Surface("BigDomain") = {surfaces[0]};       // The main big domain surface
Physical Surface("LowerSmallBox") = {surfaces[1]};   // Lower small box surface
Physical Surface("UpperSmallBox") = {surfaces[2]};   // Upper small box surface

// Assign physical lines for boundary conditions
Physical Line("Boundary") = {1, 2, 3, 4};
Physical Line("SharedLine") = {7, 9}; // Lines along y=0

// Generate the mesh
Mesh 2;

// Elevate the mesh order to 2
SetOrder 2;

// Optionally optimize the high-order mesh
OptimizeMesh "HighOrder";



