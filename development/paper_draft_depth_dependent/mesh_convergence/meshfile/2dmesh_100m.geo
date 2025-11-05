SetFactory("OpenCASCADE");

// Define mesh sizes
lc_fault = 100;
lc = 1000;

// Define the big square (2D)
big_xmin = -30000;
big_xmax = 30000;
big_ymin = -30000;
big_ymax = 30000;

// Define the small box (2D)
small_xmin = -20000;
small_xmax = 20000;
small_ymin = -2000;
small_ymax = 2000;

// Define the embedded fault line (2D)
// Line from (-15000, 0) to (15000, 0)
fault_xmin = -15000;
fault_xmax = 15000;
fault_y = 0;

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

// Define points for the small box
Point(5) = {small_xmin, small_ymin, 0, lc_fault};
Point(6) = {small_xmax, small_ymin, 0, lc_fault};
Point(7) = {small_xmax, small_ymax, 0, lc_fault};
Point(8) = {small_xmin, small_ymax, 0, lc_fault};

// Define lines for the small box
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Define points for the embedded fault line
Point(9) = {fault_xmin, fault_y, 0, lc_fault};
Point(10) = {fault_xmax, fault_y, 0, lc_fault};

// Define the embedded fault line
Line(9) = {9, 10};

// Define points for the nucleation patch
// Point(13) = {nucl_xmin, nucl_ymin, 0, lc_fault};
// Point(14) = {nucl_xmax, nucl_ymin, 0, lc_fault};
// Point(15) = {nucl_xmax, nucl_ymax, 0, lc_fault};
// Point(16) = {nucl_xmin, nucl_ymax, 0, lc_fault};

// Define lines for the small box
// Line(13) = {13, 14};
// Line(14) = {14,15};
// Line(15) = {15,16};
// Line(16) = {16,13};

// Create line loops
Line Loop(1) = {1, 2, 3, 4};  // Big square
Line Loop(2) = {5, 6, 7, 8};  // Small box
// Line Loop(4) = {13,14,15,16};  // Nucleation Patch

// Create surfaces for the big square and small box
Plane Surface(1) = {1};  // Big square surface
Plane Surface(2) = {2};  // Small box surface
// Plane Surface(4) = {4};  // Small box surface

// Boolean operation to fragment all surfaces and embed the fault line
BooleanFragments{ Surface{1,2}; Delete; }{ Line{9}; }

// ====================================================================
// IMPROVED: Distance-based mesh gradation from embedded fault line
// Mesh size increases smoothly from lc_fault (100) at fault line
// to lc (2000) at domain boundaries, independent of block boundaries
// ====================================================================

// 1. Create a Distance field from the embedded fault line
Field[1] = Distance;
Field[1].CurvesList = {9};  // Embedded fault line
Field[1].Sampling = 100;  // Number of points for distance calculation

// 2. Create a Threshold field for smooth mesh size transition
// Uses distance from embedded fault line to control mesh size
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_fault;           // 100 m at fault line
Field[2].LcMax = lc;                  // 2000 m at far boundaries
Field[2].DistMin = 0;                 // Start transition immediately
Field[2].DistMax = 25000;             // Complete transition by 25 km
Field[2].Sigmoid = 1;                 // Smooth sigmoid transition (not linear)

// 3. Set the Threshold field as the background field
// This overrides point-based mesh size assignments
Background Field = 2;

// 4. Allow mesh size to vary smoothly without constraint from geometry
Mesh.MeshSizeExtendFromBoundary = 0;  // Don't extend from boundary elements
Mesh.MeshSizeFromPoints = 0;           // Don't use point-based sizes
Mesh.MeshSizeFromCurvature = 0;        // Don't adapt to curvature

// Mark all surfaces as physical surfaces
surfaces[] = Surface{:};  // Collect all existing surfaces
For i In {0:#surfaces[]-1}
    Physical Surface(Sprintf("Surface_%g", i+1)) = {surfaces[i]};
EndFor

// Print the number of surfaces created
Printf("Number of surfaces created: %g", #surfaces[]);

// // Mesh the geometry
// Mesh 2;

// // Then elevate the mesh order to 2
// SetOrder 2;

// // Optionally optimize the high-order mesh
// OptimizeMesh "HighOrder";
