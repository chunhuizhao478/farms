SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc0 = 0.0001;  // Smaller values result in finer mesh
lc = 0.0001;  // Smaller values result in finer mesh
lc2 = 0.1; // Larger mesh size for outer surfaces

height = 0.1;
radius_outer       = 1;
radius_refinedbc   = 0.0645;
radius_inner       = 0.064;

// Define the main cylinder with specific mesh size at base points
Point(1) = {0, 0, 0, lc0};
Point(2) = {0, 0, height, lc0};
Cylinder(1) = {0, 0, 0, 0, 0, height, radius_outer, 2*Pi};
// Cylinder(2) = {0, 0, 0, 0, 0, height, radius_refinedbc, 2*Pi};

// Hole cylinder
Point(7) = {0, 0, 0, lc0};
Point(8) = {0, 0, height, lc0};
Cylinder(3) = {0, 0, 0, 0, 0, height, radius_inner, 2*Pi};

// Subtract the hole
// BooleanDifference{ Volume{1,2}; Delete; }{ Volume{3}; Delete; }
BooleanDifference{ Volume{1}; Delete; }{ Volume{3}; Delete; }

// Physical volumes (labels)
Physical Volume("Refined Cylinder") = {1};
// Physical Volume("Pierced Cylinder") = {3};

// Center point for Distance field
Point(3) = {0, 0, 0, lc0};

// 1) Distance from hole‐center
Field[1] = Distance;
Field[1].NodesList = {3};

// 2) Radial grading: lc (when r = radius_refinedbc) → lc2 (when r = radius_outer)
Field[2] = Threshold;
Field[2].IField   = 1;
Field[2].LcMin    = lc;     // at r = radius_refinedbc
Field[2].LcMax    = lc2;    // at r = radius_outer
Field[2].DistMin  = radius_refinedbc;
Field[2].DistMax  = radius_outer;

// 3) Inner blend: lc0 (at r = 0) → lc (at r = radius_refinedbc)
Field[3] = Threshold;
Field[3].IField   = 1;
Field[3].LcMin    = lc0;    // at r = 0
Field[3].LcMax    = lc;     // at r = radius_refinedbc
Field[3].DistMin  = 0.0;
Field[3].DistMax  = radius_refinedbc;

// 4) Pick the minimum at each point
Field[4] = Min;
Field[4].FieldsList = {2, 3};

// Activate the background field
Background Field = 4;

// Label the surfaces (no change here)
Physical Surface("Upper") = {6};
Physical Surface("Outer") = {5};
Physical Surface("Inner") = {4};
Physical Surface("Bottom") = {7};


