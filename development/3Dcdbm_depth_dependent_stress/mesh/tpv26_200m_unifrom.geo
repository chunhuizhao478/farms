/**
 * Poroelastic fault mesh: pure tet4 everywhere.
 * - Inner region: structured uniform tet4 (transfinite brick, split by fault).
 * - Outer region: graded unstructured tet4, coarsening outward.
 * - SINGLE LAYER VERSION: Only one structured layer around fault.
 *
 * CORRECTED: Proper fault surface identification
 */

SetFactory("OpenCASCADE");

// ----------------------------------------------------
// GLOBAL & LOCAL SIZES  (COARSER, ADJUST IF NEEDED)
// ----------------------------------------------------
lc       = 2e4;    // global coarse size away from fault
lc_fault = 200;    // target fine size near fault (inner region)
Fault_length        = 45e3;
Fault_width         = 20e3;
Fault_dip           = 90*Pi/180.;
transition_length   = 800; //0.25e3;   // inner region halo around fault <- change it back if you want smaller inner volume

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

// ----------------------------------------------------
// INNER PRISM AROUND FAULT (STRUCTURED REGION)
// ----------------------------------------------------
X_inner_min = -0.5*Fault_length - 0.5*transition_length;
X_inner_max =  0.5*Fault_length + 0.5*transition_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - 0.5*transition_length;

// Compute divisions for inner region to get ~lc_fault spacing
dx_inner = X_inner_max - X_inner_min;
dy_inner = Y_inner_max - Y_inner_min;
dz_inner = Abs(Z_inner_bot - Z_inner_top);

n_x = Round(dx_inner / lc_fault);
n_z = Round(dz_inner / lc_fault);
n_x = (n_x < 1) ? 1 : n_x;
n_z = (n_z < 1) ? 1 : n_z;

// Split Y at fault (Y=0)
n_y_below = Round(Abs(Y_inner_min) / lc_fault);
n_y_above = Round(Abs(Y_inner_max) / lc_fault);
n_y_below = (n_y_below < 1) ? 1 : n_y_below;
n_y_above = (n_y_above < 1) ? 1 : n_y_above;

Printf("Inner region divisions (tet4): nx=%g, ny_below=%g, ny_above=%g, nz=%g",
       n_x, n_y_below, n_y_above, n_z);

// ----------------------------------------------------
// BUILD INNER REGION SPLIT AT Y=0 (FAULT PLANE)
// ----------------------------------------------------

// Create SINGLE box containing both regions (Y < 0 and Y > 0)
// This ensures all elements share nodes - critical for BreakMeshByBlockGenerator
Box(100) = {X_inner_min, Y_inner_min, Z_inner_bot,
            dx_inner, Y_inner_max - Y_inner_min, Abs(Z_inner_bot)}; // Entire inner region

// MOOSE will split this into blocks 100 and 200 using ParsedSubdomainMeshGenerator (based on Y coordinate)
// Then BreakMeshByBlockGenerator will duplicate nodes at Y=0 to create the fault interface

// ----------------------------------------------------
// TRANSFINITE CONSTRAINTS FOR INNER VOLUMES (TET4)
// ----------------------------------------------------

// --- Single volume with uniform meshing ---
// Total Y divisions = below + above
n_y_total = n_y_below + n_y_above;

edges_100[] = Unique(Abs(Boundary{ Surface{Boundary{ Volume{100}; }}; }));

For i In {0:#edges_100[]-1}
  e  = edges_100[i];
  bb[] = BoundingBox Curve{e};
  dx = Abs(bb[3] - bb[0]);
  dy = Abs(bb[4] - bb[1]);
  dz = Abs(bb[5] - bb[2]);

  If (dx > dy && dx > dz)
    Transfinite Curve{e} = n_x + 1;
  ElseIf (dy > dx && dy > dz)
    Transfinite Curve{e} = n_y_total + 1;  // Total Y divisions
  Else
    Transfinite Curve{e} = n_z + 1;
  EndIf
EndFor

surfs_100[] = Unique(Abs(Boundary{ Volume{100}; }));
Transfinite Surface{surfs_100[]};
Transfinite Volume{100};   // structured tet4 brick

// ----------------------------------------------------
// NOTE: Fault interface will be created by MOOSE BreakMeshByBlockGenerator
// No need to identify fault surfaces here since the volume is continuous
// ----------------------------------------------------

// ----------------------------------------------------
// OUTER REGION
// ----------------------------------------------------

Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// Subtract inner from outer (no intermediate layer)
out_outer[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{100}; };
outer_vol   = out_outer[0];
Printf("Outer volume ID: %g", outer_vol);

// ----------------------------------------------------
// MESH SIZE FIELDS  (AGGRESSIVE COARSENING IN OUTER REGION, CAPPED AT lc)
// ----------------------------------------------------

// Get all surfaces of inner volume (boundaries with outer region)
inner_surfaces[] = Unique(Abs(Boundary{ Volume{100}; }));

// Distance to inner region (for refinement near fault)
Field[1] = Distance;
Field[1].SurfacesList     = {inner_surfaces[]};
Field[1].NumPointsPerCurve = 100;

// Aggressive exponential function for rapid coarsening, capped at lc
Field[2] = MathEval;
Field[2].F = Sprintf("Min(0.5*F1 + (F1/1.0e3)^3 + %g, %g)", lc_fault, lc);

// Fast coarsening right after fault - shorter distance range
Field[6] = Threshold;
Field[6].IField  = 1;
Field[6].LcMin   = lc_fault;
Field[6].LcMax   = lc * 0.5;      // Intermediate size
Field[6].DistMin = 2*lc_fault;  // Start sooner
Field[6].DistMax = 3*lc_fault;    // Finish sooner

// Very aggressive coarsening to outer boundary - capped at lc
Field[8] = Threshold;
Field[8].IField  = 1;
Field[8].LcMin   = lc * 0.5;      // Start from intermediate
Field[8].LcMax   = lc;            // Cap at lc (20km)
Field[8].DistMin = 2*lc_fault;    // Start sooner
Field[8].DistMax = 6*lc_fault;    // Reach lc quickly

// Combine all mesh-size fields
Field[7] = Min;
Field[7].FieldsList = {2, 6, 8};
Background Field = 7;

Mesh.MeshSizeFromCurvature      = 0;
Mesh.MeshSizeExtendFromBoundary = 1;

// Enforce fine size at all inner-region points
Characteristic Length { PointsOf{ Volume{100}; } } = lc_fault;

// ----------------------------------------------------
// PHYSICAL GROUPS
// ----------------------------------------------------

// Get all unique surfaces from each region
outer_boundary_surfs[] = Unique(Abs(Boundary{ Volume{outer_vol}; }));
inner_boundary_surfs[] = Unique(Abs(Boundary{ Volume{100}; }));

Physical Surface(101) = {outer_boundary_surfs[]};           // Outer boundary (all surfaces)
Physical Surface(105) = {inner_boundary_surfs[]};           // Inner region boundary

Physical Volume(10) = {outer_vol};     // Outer volume (coarse, unstructured)
Physical Volume(12) = {100};           // Inner volume (fine, structured) - will be split by MOOSE

// ----------------------------------------------------
// FINAL MESH SETTINGS: PURE TET4
// ----------------------------------------------------

Mesh.Algorithm  = 6;   // 2D: Frontal
Mesh.Algorithm3D = 4;  // 3D: Frontal Delaunay (tet)

Mesh.RecombineAll = 0; // absolutely no hex/quad recombination

Mesh.Optimize  = 1;
Mesh.Smoothing = 5;

Mesh.MshFileVersion = 2.2;
