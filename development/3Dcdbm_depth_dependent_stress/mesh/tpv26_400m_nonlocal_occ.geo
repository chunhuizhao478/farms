/**
 * Simple uniform TET4 mesh - Inner region only
 * No fault, no outer volume, just a uniformly meshed box
 */

SetFactory("OpenCASCADE");

// ============================================
// MESH SIZE
// ============================================
lc_uniform = 400;  // Uniform element size throughout (adjust as needed)

// ============================================
// GEOMETRY PARAMETERS
// ============================================
Fault_length = 45e3;
Fault_width = 20e3;
transition_length = 1.2e3;

// Inner box dimensions (same as your inner region)
X_min = -0.5*Fault_length - transition_length;  // -24 km
X_max =  0.5*Fault_length + transition_length;  //  24 km
Y_min = -transition_length;                      // -1.5 km
Y_max =  transition_length;                      //  1.5 km
Z_min = -Fault_width - transition_length;        // -21.5 km
Z_max = 0;                                       //   0 km

// ============================================
// CREATE SIMPLE BOX
// ============================================
Box(1) = {X_min, Y_min, Z_min, 
          X_max-X_min, Y_max-Y_min, Z_max-Z_min};

// ============================================
// UNIFORM MESH SIZING
// ============================================
// Set uniform characteristic length on all points
Characteristic Length{ PointsOf{ Volume{1}; } } = lc_uniform;

// ============================================
// MESH ALGORITHM SETTINGS
// ============================================

// 3D meshing algorithm
Mesh.Algorithm3D = 4;         // Frontal Delaunay for quality
                              // Options: 1=Delaunay, 4=Frontal, 10=HXT

// Quality optimization
Mesh.Optimize = 1;            // Enable optimization
Mesh.OptimizeNetgen = 1;      // Use Netgen optimizer for better quality
Mesh.OptimizeThreshold = 0.3; // Optimize poor quality elements

// Smoothing
Mesh.Smoothing = 10;          // Laplacian smoothing iterations

// Element order
Mesh.ElementOrder = 1;        // Linear TET4 elements

// 2D algorithm (for surfaces)
Mesh.Algorithm = 6;           // Frontal Delaunay

// Quality metric
Mesh.QualityType = 2;         // gamma (equilateral measure)

// Output
Mesh.MshFileVersion = 2.2;

// ============================================
// PHYSICAL GROUPS
// ============================================
Physical Surface(1) = Boundary{ Volume{1}; };  // All boundaries
Physical Volume(1) = {1};                      // The volume

// ============================================
// USAGE
// ============================================
// Generate mesh:
//   gmsh inner_uniform.geo -3 -o inner_mesh.msh
//
// Visualize:
//   gmsh inner_uniform.geo
//   Press '3' to mesh
//
// Adjust element size:
//   lc_uniform = 200;   // Finer (200m)
//   lc_uniform = 400;   // Default (400m)
//   lc_uniform = 800;   // Coarser (800m)
//
// Box dimensions:
//   X: 48 km (centered at x=0)
//   Y: 3 km (centered at y=0)
//   Z: 21.5 km (from surface to depth)
//
// Expected mesh:
//   Completely uniform ~400m TET4 elements
//   Optimized for equilateral quality
