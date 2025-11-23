/**
 * Poroelastic fault mesh: pure tet4 everywhere.
 * - Inner region: structured uniform tet4 (transfinite brick, split by fault).
 * - Intermediate region: structured uniform tet4 (transfinite brick, split by fault).
 * - Outer region: graded unstructured tet4, coarsening outward.
 * - Same geometry and size fields as original hex version, but no recombine.
 * 
 * CORRECTED: Proper fault surface identification
 */

SetFactory("OpenCASCADE");

// ----------------------------------------------------
// GLOBAL & LOCAL SIZES  (COARSER, ADJUST IF NEEDED)
// ----------------------------------------------------
lc       = 1.25e4;    // global coarse size away from fault
lc_fault = 150;    // target fine size near fault (inner region)
lc_intermediate = 150;  // target size in intermediate region
Fault_length        = 45e3;
Fault_width         = 20e3;
Fault_dip           = 90*Pi/180.;
transition_length   = 0.3e3;   // inner region halo around fault
intermediate_length = 0.3e3;     // intermediate transition zone

// Nucleation in X,Z local coordinates
X_nucl     = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl     = 1e3;
lc_nucl    = 150;   // nucleation refinement

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

// ----------------------------------------------------
// INNER PRISM AROUND FAULT (STRUCTURED REGION)
// ----------------------------------------------------
X_inner_min = -0.5*Fault_length - transition_length;
X_inner_max =  0.5*Fault_length + transition_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - transition_length;

// Intermediate prism (transition region)
X_mid_min = -0.5*Fault_length - transition_length - intermediate_length;
X_mid_max =  0.5*Fault_length + transition_length + intermediate_length;
Y_mid_min = -transition_length - intermediate_length;
Y_mid_max =  transition_length + intermediate_length;
Z_mid_top = 0;
Z_mid_bot = -Fault_width - transition_length - intermediate_length;

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

// Compute divisions for intermediate region to get ~lc_intermediate spacing
dx_mid = X_mid_max - X_mid_min;
dy_mid = Y_mid_max - Y_mid_min;
dz_mid = Abs(Z_mid_bot - Z_mid_top);

// Number of divisions in the intermediate shell layers
n_x_mid_layers = Round(intermediate_length / lc_intermediate);
n_y_mid_layers = Round(intermediate_length / lc_intermediate);
n_z_mid_layers = Round(intermediate_length / lc_intermediate);
n_x_mid_layers = (n_x_mid_layers < 1) ? 1 : n_x_mid_layers;
n_y_mid_layers = (n_y_mid_layers < 1) ? 1 : n_y_mid_layers;
n_z_mid_layers = (n_z_mid_layers < 1) ? 1 : n_z_mid_layers;

// Total divisions for intermediate volumes
n_x_mid = n_x + 2*n_x_mid_layers;
n_z_mid = n_z + 2*n_z_mid_layers;
n_y_mid_below = n_y_below + n_y_mid_layers;
n_y_mid_above = n_y_above + n_y_mid_layers;

Printf("Intermediate region divisions (tet4): nx=%g, ny_below=%g, ny_above=%g, nz=%g",
       n_x_mid, n_y_mid_below, n_y_mid_above, n_z_mid);

// ----------------------------------------------------
// BUILD INNER REGION SPLIT AT Y=0 (FAULT PLANE)
// ----------------------------------------------------

// Two rectangular boxes that meet at Y=0
Box(100) = {X_inner_min, Y_inner_min, Z_inner_bot,
            dx_inner, Abs(Y_inner_min), Abs(Z_inner_bot)}; // Below fault (Y<0)

Box(101) = {X_inner_min, 0, Z_inner_bot,
            dx_inner, Y_inner_max, Abs(Z_inner_bot)};      // Above fault (Y>0)

// ----------------------------------------------------
// BUILD INTERMEDIATE REGION SPLIT AT Y=0 (FAULT PLANE)
// ----------------------------------------------------

// Two rectangular boxes that meet at Y=0
Box(200) = {X_mid_min, Y_mid_min, Z_mid_bot,
            dx_mid, Abs(Y_mid_min), Abs(Z_mid_bot)}; // Below fault (Y<0)

Box(201) = {X_mid_min, 0, Z_mid_bot,
            dx_mid, Y_mid_max, Abs(Z_mid_bot)};      // Above fault (Y>0)

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

// Store both fault surfaces for later use
fault_main_surfaces[] = {fault_100, fault_101};

// ----------------------------------------------------
// TRANSFINITE CONSTRAINTS FOR INTERMEDIATE VOLUMES (TET4)
// ----------------------------------------------------

// --- Volume 200 (below fault) ---
edges_200[] = Unique(Abs(Boundary{ Surface{Boundary{ Volume{200}; }}; }));

For i In {0:#edges_200[]-1}
  e  = edges_200[i];
  bb[] = BoundingBox Curve{e};
  dx = Abs(bb[3] - bb[0]);
  dy = Abs(bb[4] - bb[1]);
  dz = Abs(bb[5] - bb[2]);

  If (dx > dy && dx > dz)
    Transfinite Curve{e} = n_x_mid + 1;
  ElseIf (dy > dx && dy > dz)
    Transfinite Curve{e} = n_y_mid_below + 1;
  Else
    Transfinite Curve{e} = n_z_mid + 1;
  EndIf
EndFor

surfs_200[] = Unique(Abs(Boundary{ Volume{200}; }));
Transfinite Surface{surfs_200[]};
Transfinite Volume{200};   // structured tet4 brick

// --- Volume 201 (above fault) ---
edges_201[] = Unique(Abs(Boundary{ Surface{Boundary{ Volume{201}; }}; }));

For i In {0:#edges_201[]-1}
  e  = edges_201[i];
  bb[] = BoundingBox Curve{e};
  dx = Abs(bb[3] - bb[0]);
  dy = Abs(bb[4] - bb[1]);
  dz = Abs(bb[5] - bb[2]);

  If (dx > dy && dx > dz)
    Transfinite Curve{e} = n_x_mid + 1;
  ElseIf (dy > dx && dy > dz)
    Transfinite Curve{e} = n_y_mid_above + 1;
  Else
    Transfinite Curve{e} = n_z_mid + 1;
  EndIf
EndFor

surfs_201[] = Unique(Abs(Boundary{ Volume{201}; }));
Transfinite Surface{surfs_201[]};
Transfinite Volume{201};   // structured tet4 brick

// ----------------------------------------------------
// SUBTRACT INNER FROM INTERMEDIATE TO CREATE SHELL
// ----------------------------------------------------

out_mid[] = BooleanDifference{ Volume{200, 201}; Delete; }{ Volume{100, 101}; };
mid_vols[] = out_mid[];
For i In {0:#mid_vols[]-1}
  Printf("Intermediate volume ID[%g]: %g", i, mid_vols[i]);
EndFor

// ----------------------------------------------------
// NUCLEATION PATCH AS SURFACE IN FAULT PLANE
// ----------------------------------------------------

Rectangle(1001) = {X_nucl - R_nucl, Width_nucl - R_nucl, 0,
                   2*R_nucl, 2*R_nucl};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} {
  Surface{1001};
}
fault_nucl = 1001;

// Embed nucleation surface in the inner volumes (constrains mesh at nucleation)
Surface{fault_nucl} In Volume{100};
Surface{fault_nucl} In Volume{101};

// Fault surfaces for physical groups - nucleation plus main fault surfaces
fault_surfaces[] = {fault_nucl, fault_100, fault_101};

// ----------------------------------------------------
// OUTER REGION
// ----------------------------------------------------

Box(1) = {Xmin, Ymin, Zmin,
          Xmax - Xmin, Ymax - Ymin, -Zmin};

// Subtract intermediate from outer
out_outer[] = BooleanDifference{ Volume{1}; Delete; }{ Volume{mid_vols[]}; };
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

// Distance to nucleation patch
Field[3] = Distance;
Field[3].SurfacesList = {fault_nucl};
Field[3].NumPointsPerCurve = 50;

// Local refinement near nucleation
Field[4] = Threshold;
Field[4].IField  = 3;
Field[4].LcMin   = lc_nucl;
Field[4].LcMax   = lc_fault;
Field[4].DistMin = 0;
Field[4].DistMax = 2*R_nucl;

// Restrict nucleation refinement to nucleation surface only
Field[5] = Restrict;
Field[5].IField       = 4;
Field[5].SurfacesList = {fault_nucl};

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
Field[7].FieldsList = {2, 5, 6, 8};
Background Field = 7;

Mesh.MeshSizeFromCurvature      = 0;
Mesh.MeshSizeExtendFromBoundary = 1;

// Enforce fine size at all inner-region points
Characteristic Length { PointsOf{ Volume{100, 101}; } } = lc_fault;

// Enforce intermediate size at all intermediate-region points
Characteristic Length { PointsOf{ Volume{mid_vols[]}; } } = lc_intermediate;

// ----------------------------------------------------
// PHYSICAL GROUPS
// ----------------------------------------------------

// Get all unique surfaces from each region
outer_boundary_surfs[] = Unique(Abs(Boundary{ Volume{outer_vol}; }));
inner_boundary_surfs[] = Unique(Abs(Boundary{ Volume{100, 101}; }));
mid_boundary_surfs[] = Unique(Abs(Boundary{ Volume{mid_vols[]}; }));

// Remove fault surfaces from inner boundary (they should only be in Physical Surface 103)
inner_boundary_filtered[] = {};
For i In {0:#inner_boundary_surfs[]-1}
  If (inner_boundary_surfs[i] != fault_100 && 
      inner_boundary_surfs[i] != fault_101 && 
      inner_boundary_surfs[i] != fault_nucl)
    inner_boundary_filtered[] += inner_boundary_surfs[i];
  EndIf
EndFor

// Remove fault surfaces from intermediate boundary
mid_boundary_filtered[] = {};
For i In {0:#mid_boundary_surfs[]-1}
  If (mid_boundary_surfs[i] != fault_100 && 
      mid_boundary_surfs[i] != fault_101 && 
      mid_boundary_surfs[i] != fault_nucl)
    mid_boundary_filtered[] += mid_boundary_surfs[i];
  EndIf
EndFor

Physical Surface(101) = {outer_boundary_surfs[]};           // Outer boundary (all surfaces)
Physical Surface(103) = {fault_surfaces[]};                 // Fault (nucleation + main fault surfaces)
Physical Surface(105) = {inner_boundary_filtered[]};        // Inner region boundary (excluding fault)
Physical Surface(106) = {mid_boundary_filtered[]};          // Intermediate boundary (excluding fault)

Physical Volume(10) = {outer_vol};     // Outer volume
Physical Volume(11) = {mid_vols[]};    // Intermediate volumes (structured tet4)
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
