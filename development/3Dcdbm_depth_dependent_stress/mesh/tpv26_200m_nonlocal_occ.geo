/**
 * Hybrid meshing approach:
 * - Inner volume: Structured HEXAHEDRAL mesh (uniform cubes)
 * - Outer volume: Adaptive TETRAHEDRAL mesh (TET4)
 * NOTE: Fault surfaces only affect outer volume meshing
 */
SetFactory("OpenCASCADE");

lc = 2e4;
lc_fault = 200;
Fault_length = 45e3;
Fault_width = 20e3;
Fault_dip = 90*Pi/180.;
transition_length = 1.5e3;

X_nucl = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl = 1e3;
lc_nucl = 200;

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax = Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

// Inner prism dimensions
X_inner_min = -0.5*Fault_length - transition_length;
X_inner_max = 0.5*Fault_length + transition_length;
Y_inner_min = -transition_length;
Y_inner_max = transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - transition_length;

// IMPORTANT: Define desired cube size for inner region
cube_size = 500; // 500m uniform cubes - adjust as needed

// -------------------------------------------
// GEOMETRY DEFINITION
// -------------------------------------------

// Create outer box
Box(1) = {Xmin, Ymin, Zmin, Xmax-Xmin, Ymax-Ymin, -Zmin};

// Create inner box (will be meshed with hex)
Box(2) = {X_inner_min, Y_inner_min, Z_inner_bot,
          X_inner_max-X_inner_min, Y_inner_max-Y_inner_min, Z_inner_top-Z_inner_bot};

// Create fault surfaces (will only cut outer volume)
Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }

Rectangle(1001) = {X_nucl-R_nucl, Width_nucl-R_nucl, 0, 2*R_nucl, 2*R_nucl};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1001}; }

// Boolean: Remove inner from outer, then fragment outer by fault
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; };
BooleanFragments{ Volume{3}; Delete; }{ Surface{1000,1001}; }

// Volumes: inner_vol=2 (untouched box), shell_vol from fragmentation
shell_vol = 3;
inner_vol = 2;
fault_main = 1000;
fault_nucl = 1001;

// -------------------------------------------
// INNER VOLUME: STRUCTURED HEX MESH
// -------------------------------------------

// Calculate number of divisions for uniform cubes
nx = Round((X_inner_max - X_inner_min) / cube_size);
ny = Round((Y_inner_max - Y_inner_min) / cube_size);
nz = Abs(Round((Z_inner_top - Z_inner_bot) / cube_size));

// Ensure at least 2 divisions per direction
If (nx < 2)
  nx = 2;
EndIf
If (ny < 2)
  ny = 2;
EndIf
If (nz < 2)
  nz = 2;
EndIf

// Apply transfinite meshing to inner volume
// Get curves of inner box and set transfinite divisions
curves_inner[] = Curve{ Volume{inner_vol}; };

// X-direction curves (4 curves)
Transfinite Curve{curves_inner[0], curves_inner[2], curves_inner[4], curves_inner[6]} = nx + 1;
// Y-direction curves (4 curves)  
Transfinite Curve{curves_inner[1], curves_inner[3], curves_inner[5], curves_inner[7]} = ny + 1;
// Z-direction curves (4 curves)
Transfinite Curve{curves_inner[8], curves_inner[9], curves_inner[10], curves_inner[11]} = nz + 1;

// Make all surfaces of inner volume transfinite
surfs_inner[] = Surface{ Volume{inner_vol}; };
Transfinite Surface{surfs_inner[]};
Recombine Surface{surfs_inner[]}; // Ensure quad faces

// Make inner volume transfinite and recombine to hex
Transfinite Volume{inner_vol};
Recombine Volume{inner_vol}; // Create hexahedral elements

// -------------------------------------------
// OUTER VOLUME: ADAPTIVE TET MESH
// -------------------------------------------

Field[1] = Distance;
Field[1].SurfacesList = {fault_main, fault_nucl};

Field[2] = MathEval;
Field[2].F = Sprintf("0.1*F1 +(F1/2.5e3)^2 + %g", lc_fault);

Field[3] = Distance;
Field[3].SurfacesList = {fault_nucl};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc_nucl;
Field[4].LcMax = lc_fault;
Field[4].DistMin = R_nucl;
Field[4].DistMax = 2*R_nucl;

Field[5] = Restrict;
Field[5].IField = 4;
Field[5].SurfacesList = {fault_main, fault_nucl};

Field[6] = Threshold;
Field[6].IField = 1;
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 2*lc_fault;
Field[6].DistMax = 2*lc_fault + 0.001;

Field[7] = Min;
Field[7].FieldsList = {2,5,6};

// Restrict background field to outer volume only
Field[8] = Restrict;
Field[8].IField = 7;
Field[8].VolumesList = {shell_vol};

Background Field = 8;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------
Physical Surface(101) = Boundary{ Volume{shell_vol}; };
Physical Surface(103) = {fault_main, fault_nucl};
Physical Surface(105) = Boundary{ Volume{inner_vol}; };
Physical Volume(10) = {shell_vol};  // Outer (TET4)
Physical Volume(11) = {inner_vol};  // Inner (HEX8)

// -------------------------------------------
// MESH SETTINGS
// -------------------------------------------
Mesh.Algorithm = 6;          // Frontal Delaunay for tets
Mesh.Algorithm3D = 1;        // Delaunay 3D
Mesh.SubdivisionAlgorithm = 1; // All hex subdivision
Mesh.RecombineAll = 0;       // Only recombine where specified
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 1;

Printf("Inner volume meshed with %g x %g x %g = %g hexahedral elements", nx, ny, nz, nx*ny*nz);
