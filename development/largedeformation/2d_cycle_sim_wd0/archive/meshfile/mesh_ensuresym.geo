SetFactory("OpenCASCADE");

// Mesh sizes
lc_big = 200;  // Coarse mesh size for the big box
lc_small = 100;  // Fine mesh size for the small box

// Big box dimensions
big_xmin = -10000;
big_xmax = 10000;
big_ymin = -10000;
big_ymax = 10000;

// Small box dimensions
small_xmin = -4000;
small_xmax = 4000;
small_ymin = -2000;
small_ymax = 2000;

// Define points for the big box
Point(1) = {big_xmin, big_ymin, 0, lc_big};
Point(2) = {big_xmax, big_ymin, 0, lc_big};
Point(3) = {big_xmax, big_ymax, 0, lc_big};
Point(4) = {big_xmin, big_ymax, 0, lc_big};

// Define lines for the big box
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define the surface for the big box
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Define points for the small box
Point(5) = {small_xmin, small_ymin, 0, lc_small};
Point(6) = {small_xmax, small_ymin, 0, lc_small};
Point(7) = {small_xmax, small_ymax, 0, lc_small};
Point(8) = {small_xmin, small_ymax, 0, lc_small};

// Define lines for the small box
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Define the surface for the small box
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// Perform BooleanFragments to split the big box into separate domains
BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete; }

// Collect the fragmented surfaces
surfaces[] = Surface{:};

// Define physical groups
Physical Surface("BigBox") = {surfaces[0]}; // Remaining big box domain
Physical Surface("SmallBox") = {surfaces[1]}; // Small box

// Assign physical lines for boundary conditions
Physical Line("Boundary") = {1, 2, 3, 4};
Physical Line("SmallBoxEdges") = {5, 6, 7, 8};

// Apply transfinite mesh on both boxes
// Big box
Transfinite Line{1, 2, 3, 4} = 50; // Adjust divisions as needed
Transfinite Surface{surfaces[0]} = {1, 2, 3, 4};

// Small box
Transfinite Line{5, 6, 7, 8} = 40; // Adjust divisions as needed
Transfinite Surface{surfaces[1]} = {5, 6, 7, 8};

// Recombine to create quadrilateral elements
Recombine Surface{surfaces[0], surfaces[1]};

// Mesh generation
Mesh 2;

// Optional: Elevate mesh order to quadratic
SetOrder 2;

// Optional: Optimize the mesh for high-order elements
OptimizeMesh "HighOrder";