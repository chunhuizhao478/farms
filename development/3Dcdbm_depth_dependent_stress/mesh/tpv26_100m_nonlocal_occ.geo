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
lc_fault = 100; // fine size near fault

Fault_length = 45e3;
Fault_width = 20e3;
Fault_dip = 90*Pi/180.;
transition_length = 1e3; // 1.5 km halo around fault

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

// 2. Create fault plane
Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }

// 3. Fragment both volumes by the fault surface (keep surface, delete original volumes)
BooleanFragments{ Volume{1,2}; Delete; }{ Surface{1000}; }

// 4. Collect resulting volumes (expect 2: outer shell & inner prism cut by fault)
vols[] = Volume{:};
// Assume ordering: first = outer shell, second = inner (verify in GUI if unsure)
shell_vol = vols[0];
inner_vol = vols[1];

// Assign fault surface tag
fault_main = 1000;

// -------------------------------------------
// MESH SETTINGS
// -------------------------------------------

// Distance to fault surface
Field[1] = Distance;
Field[1].SurfacesList = {fault_main};

// Smooth growth away from fault
Field[2] = MathEval;
Field[2].F = Sprintf("0.1*F1 +(F1/2.5e3)^2 + %g", lc_fault);

// Propagation zone sizing transition away from fault
Field[6] = Threshold;
Field[6].IField = 1; // based on distance to fault
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault + 0.001; // tiny offset to avoid zero interval

// Combine all
Field[7] = Min;
Field[7].FieldsList = {2,6};
Background Field = 7;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------

// Define physical groups
Physical Surface(101) = Boundary{ Volume{shell_vol}; };  // Outer boundary
Physical Surface(103) = {fault_main};                    // Fault surface
Physical Surface(105) = Boundary{ Volume{inner_vol}; };  // Inner volume boundary

Physical Volume(10) = {shell_vol};  // Outer volume
Physical Volume(11) = {inner_vol};  // Inner volume

// Final settings
Mesh.Algorithm = 6;  // Frontal Delaunay
Mesh.MshFileVersion = 2.2;
