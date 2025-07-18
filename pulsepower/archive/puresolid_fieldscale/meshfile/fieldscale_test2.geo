SetFactory("OpenCASCADE");

// Characteristic length (mesh size)
lc0 = 0.0001;  // Smaller values result in finer mesh
lc = 0.0001;  // Smaller values result in finer mesh
lc2 = 0.005; // Larger mesh size for outer surfaces

height = 0.005;
radius_outer = 0.125;
radius_refinedbc = 0.0205; //0.0
radius_inner = 0.02;

// Define the main cylinder with specific mesh size at base points
Point(1) = {0, 0, 0, lc0};
Point(2) = {0, 0, height, lc0};
Cylinder(1) = {0, 0, 0, 0, 0, height, radius_outer, 2*Pi};

// Define the main cylinder with specific mesh size at base points
Cylinder(2) = {0, 0, 0, 0, 0, height, radius_refinedbc, 2*Pi};

// Define the hole cylinder
Point(7) = {0, 0, 0, lc0};
Point(8) = {0, 0, height, lc0};
Cylinder(3) = {0, 0, 0, 0, 0, height, radius_inner, 2*Pi};

// Boolean operation to subtract the hole from the main cylinder
BooleanDifference{ Volume{1,2}; Delete; }{ Volume{3}; Delete; }

// Define Physical Volume for simulation or further meshing
Physical Volume("Refined Cylinder") = {2};
Physical Volume("Pierced Cylinder") = {3};

// Define a point at the center of the cylindrical hole base
Point(3) = {0, 0, 0, lc0};

// Create a mesh field that varies the mesh size radially from the center of the cylindrical hole
Field[1] = Distance;
Field[1].NodesList = {3}; // Center of the cylinder base

//--------------------------------------
// 7) Field[2] = Threshold to grade mesh from radius_refinedbc → radius_outer
//    • If Dist ≤ radius_refinedbc, Lc = lc  (the “fine zone” outside the super‐fine inner core)
//    • If Dist ≥ radius_outer, Lc = lc2 (the coarsest zone)
//    • In between, interpolate linearly from lc → lc2.
//--------------------------------------
Field[2] = Threshold;
Field[2].IField   = 1;             // input = Dist from Field[1]
Field[2].LcMin    = lc;            // size = lc for r ≤ radius_refinedbc
Field[2].LcMax    = lc2;           // size = lc2 for r ≥ radius_outer
Field[2].DistMin  = radius_refinedbc;
Field[2].DistMax  = radius_outer;

//--------------------------------------
// 8) Field[3] = Threshold to enforce a uniform fine size (lc0) inside r ≤ radius_refinedbc
//    • If Dist ≤ radius_refinedbc, Lc = lc0
//    • If Dist ≥ radius_refinedbc, Lc = lc   (to match the inner boundary of Field[2])
//    • In other words: blend from lc0 → lc over [r = 0 … radius_refinedbc].
//--------------------------------------
Field[3] = Threshold;
Field[3].IField   = 1;              // input = Dist from Field[1]
Field[3].LcMin    = lc0;            // size = lc0 for Dist ≤ 0
Field[3].LcMax    = lc2;             // size = lc  for Dist ≥ radius_refinedbc
Field[3].DistMin  = 0.0;            // start blending at the center
Field[3].DistMax  = radius_refinedbc;

//--------------------------------------
// 9) Field[4] = Min( Field[2], Field[3] )
//    This picks the smaller of the two at each point, ensuring:
//      • For 0 ≤ r ≤ radius_refinedbc → Field[3] = lc0 is the minimum.
//      • For radius_refinedbc < r ≤ radius_outer → Field[2] (lc…lc2) is the minimum.
//      • For r ≥ radius_outer → Field[2] = lc2, Field[3] = lc, so Field[4] = lc2.
//--------------------------------------
Field[4] = Min;
Field[4].FieldsList = {2, 3};

Background Field = 4;

// Now label the surfaces accordingly
Physical Surface("Upper") = {9}; // Label the outer surface of the main cylinder
Physical Surface("Outer") = {8}; // Label the inner surface of the hole cylinder
Physical Surface("Inner") = {7};
Physical Surface("Bottom") = {10};

