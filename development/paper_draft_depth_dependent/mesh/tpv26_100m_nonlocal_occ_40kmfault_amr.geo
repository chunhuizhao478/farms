/**
 * Derived from tpv26_100m.geo: adds THREE nested volumes (geometric partitions)
 * surrounding the vertical fault plane. Uses the OpenCASCADE kernel to ensure
 * the fault surface is properly integrated into the volume mesh.
 *
 * Volume structure (from innermost to outermost):
 * 1. Narrow fault zone (for nonlocal averaging):
 *    X in [-0.5*Fault_length - transition_length, 0.5*Fault_length + transition_length]
 *    Y in [-fault_zone_width, +fault_zone_width]  (±200m)
 *    Z in [0, -Fault_width - transition_length]
 *
 * 2. Intermediate transition zone:
 *    X in [-0.5*Fault_length - transition_length, 0.5*Fault_length + transition_length]
 *    Y in [-transition_length, +transition_length]  (±4km)
 *    Z in [0, -Fault_width - transition_length]
 *
 * 3. Outer coarse domain:
 *    Full domain
 */

SetFactory("OpenCASCADE"); // Required for Boolean operations

lc = 2e4;
lc_fault = 100; // fine size near fault

Fault_length = 40e3;
Fault_width = 20e3;
Fault_dip = 90*Pi/180.;
fault_zone_width = 200; // narrow zone around fault for nonlocal averaging (±200m)
transition_length = 4e3; // intermediate zone around fault (±4km)

// Nucleation in X,Z local coordinates
X_nucl = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl = 1e3;
lc_nucl = 100;

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width  *Cos(Fault_dip);
Zmin = -Xmax; // large depth extent (negative)

// Narrow fault zone dimensions (innermost box for nonlocal averaging)
X_fault_min = -0.5*Fault_length - transition_length;
X_fault_max =  0.5*Fault_length + transition_length;
Y_fault_min = -fault_zone_width;
Y_fault_max =  fault_zone_width;
Z_fault_top = 0;
Z_fault_bot = -Fault_width - transition_length;

// Intermediate transition zone dimensions
X_inner_min = -0.5*Fault_length - transition_length;
X_inner_max =  0.5*Fault_length + transition_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - transition_length;

// -------------------------------------------
// DEFINE ALL GEOMETRY USING OCC PRIMITIVES
// -------------------------------------------

// 1. Create three nested boxes (volumes 1, 2, 3 initially)
//    Box 1: Outer domain (coarsest mesh)
//    Box 2: Intermediate transition zone
//    Box 3: Narrow fault zone (UNIFORM 100m mesh for nonlocal averaging)
Box(1) = {Xmin, Ymin, Zmin, Xmax-Xmin, Ymax-Ymin, -Zmin};
Box(2) = {X_inner_min, Y_inner_min, Z_inner_bot,
                    X_inner_max-X_inner_min, Y_inner_max-Y_inner_min, Z_inner_top-Z_inner_bot};
Box(3) = {X_fault_min, Y_fault_min, Z_fault_bot,
                    X_fault_max-X_fault_min, Y_fault_max-Y_fault_min, Z_fault_top-Z_fault_bot};

// 2. Create fault plane and nucleation patch with high, non-conflicting surface tags
Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }
Rectangle(1001) = {X_nucl-R_nucl, Width_nucl-R_nucl, 0, 2*R_nucl, 2*R_nucl};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1001}; }

// 3. Fragment all three volumes by the fault surfaces (keep surfaces, delete original volumes)
BooleanFragments{ Volume{1,2,3}; Delete; }{ Surface{1000,1001}; }

// 4. Collect resulting volumes (expect 3: outer shell, intermediate zone, narrow fault zone)
vols[] = Volume{:};
// After fragmentation, we should get 3 volumes
// Note: Ordering depends on Gmsh's internal algorithm - verify in GUI if needed
// Typically ordered by size or creation sequence
shell_vol = vols[0];      // Outermost coarse domain
transition_vol = vols[1]; // Intermediate zone (±4km in Y)
fault_zone_vol = vols[2]; // Narrow fault zone (±200m in Y, UNIFORM 100m mesh)

// Assign fault surface tags directly (they are preserved: 1000 main, 1001 nucleation)
fault_main = 1000;
fault_nucl = 1001;

// -------------------------------------------
// MESH SETTINGS
// -------------------------------------------

// Distance to fault surfaces (analog of FacesList=101 in nonlocal file)
Field[1] = Distance;
Field[1].SurfacesList = {fault_main, fault_nucl};

// Smooth growth away from fault (analog of Field[2] in nonlocal file)
Field[2] = MathEval;
Field[2].F = Sprintf("0.1*F1 +(F1/2.5e3)^2 + %g", lc_fault);

// Distance to nucleation patch only (analog of Field[3] in nonlocal file)
Field[3] = Distance;
Field[3].SurfacesList = {fault_nucl};

// Threshold around nucleation (Field[4] in nonlocal)
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc_nucl;
Field[4].LcMax = lc_fault;
Field[4].DistMin = R_nucl;
Field[4].DistMax = 2*R_nucl;

// Restrict nucleation refinement strictly to fault surfaces (Field[5] in nonlocal)
Field[5] = Restrict;
Field[5].IField = 4;
Field[5].SurfacesList = {fault_main, fault_nucl};

// Propagation zone sizing transition away from fault (Field[6] in nonlocal)
Field[6] = Threshold;
Field[6].IField = 1; // based on distance to fault
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault + 0.001; // tiny offset to avoid zero interval

// UNIFORM mesh size in narrow fault zone (Field[8])
// Use a Box field to define the spatial region for uniform mesh (Y: ±200m)
Field[8] = Box;
Field[8].XMin = X_fault_min;
Field[8].XMax = X_fault_max;
Field[8].YMin = Y_fault_min;  // -200m
Field[8].YMax = Y_fault_max;  // +200m
Field[8].ZMin = Z_fault_bot;
Field[8].ZMax = Z_fault_top;
Field[8].VIn = lc_fault;   // UNIFORM 100m INSIDE the box
Field[8].VOut = lc;        // Large size OUTSIDE (will be overridden by Field[2] and Field[6])
Field[8].Thickness = 0;    // Sharp transition at box boundary

// Combine all mesh size fields
// Field[8] enforces uniform 100m in the narrow fault zone (Y: ±200m)
// Field[2] provides smooth growth outside this zone
// Field[6] manages transition from fine to coarse mesh
Field[7] = Min;
Field[7].FieldsList = {2,5,6,8}; // Field[8] now uses Box instead of Restrict
Background Field = 7;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------

// Define physical groups
Physical Surface(101) = Boundary{ Volume{shell_vol}; };        // Outer boundary (coarse)
Physical Surface(103) = {fault_main, fault_nucl};              // Fault surfaces
Physical Surface(105) = Boundary{ Volume{transition_vol}; };   // Intermediate zone boundary
Physical Surface(107) = Boundary{ Volume{fault_zone_vol}; };   // Narrow fault zone boundary

Physical Volume(10) = {shell_vol};        // Outer volume (coarse mesh)
Physical Volume(11) = {transition_vol};   // Intermediate transition zone
Physical Volume(12) = {fault_zone_vol};   // Narrow fault zone (UNIFORM 100m, for nonlocal averaging)

// Final settings
Mesh.Algorithm = 6;  // Frontal Delaunay
Mesh.MshFileVersion = 2.2;
