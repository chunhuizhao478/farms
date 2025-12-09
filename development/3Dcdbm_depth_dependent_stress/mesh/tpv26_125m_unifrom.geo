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
lc_fault = 125;    // target fine size near fault (inner region)
Fault_length        = 45e3;
Fault_width         = 20e3;
Fault_dip           = 90*Pi/180.;
transition_length   = 0.25e3;   // inner region halo around fault

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

// ----------------------------------------------------
// INNER PRISM AROUND FAULT (STRUCTURED REGION)
// ----------------------------------------------------
X_inner_min = -0.5*Fault_length;
X_inner_max =  0.5*Fault_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width;

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

// Two rectangular boxes that meet at Y=0
Box(100) = {X_inner_min, Y_inner_min, Z_inner_bot,
            dx_inner, Abs(Y_inner_min), Abs(Z_inner_bot)}; // Below fault (Y<0)

Box(101) = {X_inner_min, 0, Z_inner_bot,
            dx_inner, Y_inner_max, Abs(Z_inner_bot)};      // Above fault (Y>0)

// ----------------------------------------------------
// TRANSFINITE CONSTRAINTS FOR INNER VOLUMES (TET4)
// ----------------------------------------------------

// --- Volume 100 (below fault) ---
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
    Transfinite Curve{e} = n_y_below + 1;
  Else
    Transfinite Curve{e} = n_z + 1;
  EndIf
EndFor

surfs_100[] = Unique(Abs(Boundary{ Volume{100}; }));
Transfinite Surface{surfs_100[]};
Transfinite Volume{100};   // structured tet4 brick

// --- Volume 101 (above fault) ---
edges_101[] = Unique(Abs(Boundary{ Surface{Boundary{ Volume{101}; }}; }));

For i In {0:#edges_101[]-1}
  e  = edges_101[i];
  bb[] = BoundingBox Curve{e};
  dx = Abs(bb[3] - bb[0]);
  dy = Abs(bb[4] - bb[1]);
  dz = Abs(bb[5] - bb[2]);

  If (dx > dy && dx > dz)
    Transfinite Curve{e} = n_x + 1;
  ElseIf (dy > dx && dy > dz)
    Transfinite Curve{e} = n_y_above + 1;
  Else
    Transfinite Curve{e} = n_z + 1;
  EndIf
EndFor

surfs_101[] = Unique(Abs(Boundary{ Volume{101}; }));
Transfinite Surface{surfs_101[]};
Transfinite Volume{101};   // structured tet4 brick

// ----------------------------------------------------
// IDENTIFY FAULT SURFACES (CORRECTED)
// ----------------------------------------------------
// Get fault surfaces from each volume (they're at Y=0 plane)
surfs_100_all[] = Unique(Abs(Boundary{ Volume{100}; }));
surfs_101_all[] = Unique(Abs(Boundary{ Volume{101}; }));

// Find surfaces at Y=0 (the fault plane)
fault_100 = -1;
fault_101 = -1;

For i In {0:#surfs_100_all[]-1}
  s = surfs_100_all[i];
  bb[] = BoundingBox Surface{s};
  // Check if surface is at Y=0 (fault plane) - top surface of volume 100
  // Y coordinates should both be very close to 0
  If (Abs(bb[1]) < 1e-6 && Abs(bb[4]) < 1e-6)
    fault_100 = s;
  EndIf
EndFor

For i In {0:#surfs_101_all[]-1}
  s = surfs_101_all[i];
  bb[] = BoundingBox Surface{s};
  // Check if surface is at Y=0 (fault plane) - bottom surface of volume 101
  If (Abs(bb[1]) < 1e-6 && Abs(bb[4]) < 1e-6)
    fault_101 = s;
  EndIf
EndFor

Printf("Found fault surface IDs: %g (volume 100), %g (volume 101)", fault_100, fault_101);

// Store fault surfaces for later use
fault_surfaces[] = {fault_100, fault_101};

// ----------------------------------------------------
// OUTER REGION
// ----------------------------------------------------

Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// Subtract inner from outer (no intermediate layer)
out_outer[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{100, 101}; };
outer_vol   = out_outer[0];
Printf("Outer volume ID: %g", outer_vol);

// ----------------------------------------------------
// MESH SIZE FIELDS  (AGGRESSIVE COARSENING IN OUTER REGION, CAPPED AT lc)
// ----------------------------------------------------

// Distance to main fault (for general fault refinement)
Field[1] = Distance;
Field[1].SurfacesList     = {fault_100, fault_101};
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
Characteristic Length { PointsOf{ Volume{100, 101}; } } = lc_fault;

// ----------------------------------------------------
// PHYSICAL GROUPS
// ----------------------------------------------------

// Get all unique surfaces from each region
outer_boundary_surfs[] = Unique(Abs(Boundary{ Volume{outer_vol}; }));
inner_boundary_surfs[] = Unique(Abs(Boundary{ Volume{100, 101}; }));

// Remove fault surfaces from inner boundary (they should only be in Physical Surface 103)
inner_boundary_filtered[] = {};
For i In {0:#inner_boundary_surfs[]-1}
  If (inner_boundary_surfs[i] != fault_100 && 
      inner_boundary_surfs[i] != fault_101)
    inner_boundary_filtered[] += inner_boundary_surfs[i];
  EndIf
EndFor

Physical Surface(101) = {outer_boundary_surfs[]};           // Outer boundary (all surfaces)
Physical Surface(103) = {fault_surfaces[]};                 // Fault (main fault surfaces)
Physical Surface(105) = {inner_boundary_filtered[]};        // Inner region boundary (excluding fault)

Physical Volume(10) = {outer_vol};     // Outer volume
Physical Volume(12) = {100, 101};      // Inner volumes (structured tet4)

// ----------------------------------------------------
// FINAL MESH SETTINGS: PURE TET4
// ----------------------------------------------------

Mesh.Algorithm  = 6;   // 2D: Frontal
Mesh.Algorithm3D = 4;  // 3D: Frontal Delaunay (tet)

Mesh.RecombineAll = 0; // absolutely no hex/quad recombination

Mesh.Optimize  = 1;
Mesh.Smoothing = 5;

Mesh.MshFileVersion = 2.2;
