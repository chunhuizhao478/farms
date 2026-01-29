// BP2-QD Benchmark Mesh with Graded Refinement
// Domain: 400 km x 400 km (matching Tandem BP1)
// Fine mesh near fault, coarsens toward far boundaries
// ALL UNITS IN METERS
//
// Usage:
//   gmsh -2 bp2_graded.geo -o bp2_graded.msh
//   gmsh -2 bp2_graded.geo -setnumber hf 200 -o bp2_200m.msh   // 200m fault mesh
//   gmsh -2 bp2_graded.geo -setnumber hf 400 -o bp2_400m.msh   // 400m fault mesh

// Mesh parameters (in meters)
DefineConstant[ hf = {400, Min 50, Max 10000, Name "Fault resolution (m)" } ];
DefineConstant[ h_far = {20000, Min 1000, Max 100000, Name "Far boundary resolution (m)" } ];
DefineConstant[ D = {400000, Min 100000, Max 1000000, Name "Domain half-width (m)" } ];
DefineConstant[ H = {400000, Min 100000, Max 1000000, Name "Domain depth (m)" } ];

// Fault geometry (BP2 specification, in meters)
Wf = 40000;      // Rate-state fault depth (m)
d1 = 15000;      // VW region depth
d2 = 16000;      // Start of transition
d3 = 18000;      // End of transition
d4 = Wf;         // End of rate-state region

// Intermediate mesh size
h_mid = (hf + h_far) / 4;

// ============================================
// SHARED FAULT POINTS (used by both domains)
// ============================================
Point(1) = {0, 0, 0, hf};           // Surface
Point(2) = {0, -d1, 0, hf};         // End of VW region
Point(3) = {0, -d2, 0, hf};         // Start of transition
Point(4) = {0, -d3, 0, hf};         // End of transition
Point(5) = {0, -d4, 0, hf};         // End of RS fault
Point(6) = {0, -H, 0, h_mid};       // Deep boundary

// ============================================
// LEFT DOMAIN (x < 0)
// ============================================
Point(11) = {-D, 0, 0, h_far};      // Top-left corner
Point(12) = {-D, -H, 0, h_far};     // Bottom-left corner

Line(1) = {11, 1};       // Top boundary left
Line(2) = {1, 2};        // Fault segment 1
Line(3) = {2, 3};        // Fault segment 2
Line(4) = {3, 4};        // Fault segment 3
Line(5) = {4, 5};        // Fault segment 4
Line(6) = {5, 6};        // Below RS fault
Line(7) = {6, 12};       // Bottom boundary left
Line(8) = {12, 11};      // Far left boundary

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// ============================================
// RIGHT DOMAIN (x > 0)
// ============================================
Point(21) = {D, 0, 0, h_far};       // Top-right corner
Point(22) = {D, -H, 0, h_far};      // Bottom-right corner

Line(21) = {1, 21};      // Top boundary right
Line(22) = {21, 22};     // Far right boundary
Line(23) = {22, 6};      // Bottom boundary right

// Fault lines (same direction as left domain for consistency)
// Use negative line numbers to reverse direction
Curve Loop(2) = {21, 22, 23, -6, -5, -4, -3, -2};
Plane Surface(2) = {2};

// ============================================
// Physical Groups (for MOOSE)
// ============================================

// Subdomains (blocks)
Physical Surface("left_block", 1) = {1};
Physical Surface("right_block", 2) = {2};

// External boundaries
Physical Curve("left_left", 101) = {8};           // Far left boundary
Physical Curve("right_right", 102) = {22};        // Far right boundary
Physical Curve("left_top", 103) = {1};            // Top left
Physical Curve("right_top", 104) = {21};          // Top right
Physical Curve("left_bottom", 105) = {7};         // Bottom left
Physical Curve("right_bottom", 106) = {23};       // Bottom right

// Fault (internal boundary - shared by both domains)
Physical Curve("fault", 200) = {2, 3, 4, 5, 6};

// ============================================
// Mesh Grading
// ============================================

// Distance field from fault
Field[1] = Distance;
Field[1].CurvesList = {2, 3, 4, 5, 6};
Field[1].Sampling = 100;

// Threshold: fine near fault, coarse far away
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hf;
Field[2].SizeMax = h_far;
Field[2].DistMin = 0;
Field[2].DistMax = 100000;  // Transition over 100 km

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;

// ============================================
// Mesh Settings
// ============================================
Mesh.Algorithm = 6;              // Frontal-Delaunay
Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = 1;
Mesh.Smoothing = 10;
Mesh.OptimizeNetgen = 1;

// Ensure coherent mesh at internal boundary
Coherence;
