/**
 * Derived from tpv26_100m.geo: adds an embedded inner volume (geometric partition)
 * surrounding the vertical fault plane. Uses the OpenCASCADE kernel to ensure
 * the fault surface is properly integrated into the volume mesh.
 * The inner prism bounds:
 *   X in [-0.5*Fault_length - transition_length, 0.5*Fault_length + transition_length]
 *   Y in [-transition_length, +transition_length]
 *   Z in [0, -Fault_width - transition_length]
 * with transition_length = 2 km.
 */

SetFactory("OpenCASCADE"); // Required for Boolean operations

lc = 2e4;
lc_fault = 200; // uniform mesh size in inner domain (including fault)

Fault_length = 40e3;
Fault_width = 20e3;
Fault_dip = 90*Pi/180.;
transition_length = 4e3; // 2 km halo around fault

// Nucleation in X,Z local coordinates
X_nucl = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl = 1e3;
lc_nucl = 200; // uniform 200m mesh (same as inner domain)

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width  *Cos(Fault_dip);
Zmin = -Xmax; // large depth extent (negative)

// Inner prism around fault dimensions
X_inner_min = -0.5*Fault_length - transition_length;
X_inner_max =  0.5*Fault_length + transition_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - transition_length;

// -------------------------------------------
// DEFINE ALL GEOMETRY USING OCC PRIMITIVES
// -------------------------------------------

// 1. Create outer and inner boxes (volumes 1 & 2 initially)
Box(1) = {Xmin, Ymin, Zmin, Xmax-Xmin, Ymax-Ymin, -Zmin};
Box(2) = {X_inner_min, Y_inner_min, Z_inner_bot,
                    X_inner_max-X_inner_min, Y_inner_max-Y_inner_min, Z_inner_top-Z_inner_bot};

// 2. Create fault plane and nucleation patch with high, non-conflicting surface tags
Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }
Rectangle(1001) = {X_nucl-R_nucl, Width_nucl-R_nucl, 0, 2*R_nucl, 2*R_nucl};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1001}; }

// 3. Fragment both volumes by the two fault-related surfaces (keep surfaces, delete original volumes)
BooleanFragments{ Volume{1,2}; Delete; }{ Surface{1000,1001}; }

// 4. Collect resulting volumes (expect 2: outer shell & inner prism cut by fault)
vols[] = Volume{:};
// Assume ordering: first = outer shell, second = inner (verify in GUI if unsure)
shell_vol = vols[0];
inner_vol = vols[1];

// Assign fault surface tags directly (they are preserved: 1000 main, 1001 nucleation)
fault_main = 1000;
fault_nucl = 1001;

// -------------------------------------------
// MESH SETTINGS
// -------------------------------------------

// Distance to fault surfaces (analog of FacesList=101 in nonlocal file)
Field[1] = Distance;
Field[1].SurfacesList = {fault_main, fault_nucl};

// Smooth growth away from fault - ONLY for outer shell (analog of Field[2] in nonlocal file)
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

// UNIFORM 200m MESH IN INNER BOX (Field[8])
// Box field: enforces uniform mesh size within the inner volume
Field[8] = Box;
Field[8].VIn = lc_fault;  // 200m inside the box
Field[8].VOut = lc;       // 20km outside the box (coarse mesh)
Field[8].XMin = X_inner_min;
Field[8].XMax = X_inner_max;
Field[8].YMin = Y_inner_min;
Field[8].YMax = Y_inner_max;
Field[8].ZMin = Z_inner_bot;
Field[8].ZMax = Z_inner_top;
Field[8].Thickness = 0;  // Sharp transition at boundary

// Combine all: Min of outer mesh fields (2,5,6) and inner uniform field (8)
Field[7] = Min;
Field[7].FieldsList = {2,5,6,8};
Background Field = 7;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------

// Define physical groups
Physical Surface(101) = Boundary{ Volume{shell_vol}; };  // Outer boundary
Physical Surface(103) = {fault_main, fault_nucl};        // Fault surfaces
Physical Surface(105) = Boundary{ Volume{inner_vol}; };  // Inner volume boundary

Physical Volume(10) = {shell_vol};  // Outer volume
Physical Volume(11) = {inner_vol};  // Inner volume

// Final settings
Mesh.Algorithm = 6;  // Frontal Delaunay
Mesh.MshFileVersion = 2.2;
