/**
 * All-HEX: Structured cubes inner + Adaptive hex outer, NO pyramids
 */
SetFactory("OpenCASCADE");

lc = 2e4;
lc_fault = 500;
Fault_length = 45e3;
Fault_width = 20e3;
Fault_dip = 90*Pi/180.;
transition_length = 0.5e3;

Xmax = 60e3;
Xmin = -Xmax;
Ymin = -Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Ymax = Xmax + 0.5 * Fault_width * Cos(Fault_dip);
Zmin = -Xmax;

X_inner_min = -0.5*Fault_length - transition_length;
X_inner_max = 0.5*Fault_length + transition_length;
Y_inner_min = -transition_length;
Y_inner_max = transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width - transition_length;

cube_size = lc_fault;
fault_quad_size = lc_fault;

// -------------------------------------------
// INNER VOLUME - STRUCTURED CUBES
// -------------------------------------------

Lx = X_inner_max - X_inner_min;
Ly = Y_inner_max - Y_inner_min;
Lz = Z_inner_top - Z_inner_bot;

nx = Ceil(Lx / cube_size);
ny = Ceil(Ly / cube_size);
nz = Ceil(Abs(Lz) / cube_size);

Point(1) = {X_inner_min, Y_inner_min, Z_inner_bot, cube_size};
Point(2) = {X_inner_max, Y_inner_min, Z_inner_bot, cube_size};
Point(3) = {X_inner_max, Y_inner_max, Z_inner_bot, cube_size};
Point(4) = {X_inner_min, Y_inner_max, Z_inner_bot, cube_size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Transfinite Surface{1};
Recombine Surface{1};

out[] = Extrude {0, 0, Lz} {
  Surface{1};
  Layers{nz};
  Recombine;
};

inner_vol = out[1];

// -------------------------------------------
// FAULT SURFACE
// -------------------------------------------

Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }

// -------------------------------------------
// OUTER VOLUME
// -------------------------------------------

Box(100) = {Xmin, Ymin, Zmin, Xmax-Xmin, Ymax-Ymin, -Zmin};

BooleanDifference(200) = { Volume{100}; Delete; }{ Volume{inner_vol}; };
BooleanFragments{ Volume{200}; Delete; }{ Surface{1000}; }

shell_vol = 200;
fault_main = 1000;

// -------------------------------------------
// FAULT SURFACE MESHING
// -------------------------------------------

n_fault_length = Ceil(Fault_length / fault_quad_size);
n_fault_width = Ceil(Fault_width / fault_quad_size);

fault_main_curves[] = Boundary{ Surface{fault_main}; };
Transfinite Curve{fault_main_curves[0], fault_main_curves[2]} = n_fault_length + 1;
Transfinite Curve{fault_main_curves[1], fault_main_curves[3]} = n_fault_width + 1;
Transfinite Surface{fault_main};
Recombine Surface{fault_main};

// -------------------------------------------
// MESH SIZE FIELDS FOR OUTER VOLUME ONLY
// -------------------------------------------

Field[1] = Distance;
Field[1].SurfacesList = {fault_main};

Field[2] = MathEval;
Field[2].F = Sprintf("0.15*F1 + (F1/2.5e3)^2 + %g", lc_fault);

Field[3] = Threshold;
Field[3].IField = 1;
Field[3].LcMin = lc_fault;
Field[3].LcMax = lc;
Field[3].DistMin = lc_fault;
Field[3].DistMax = 2*lc_fault + 0.001;

Field[4] = Min;
Field[4].FieldsList = {2, 3};

// CRITICAL: Only apply to outer volume
Field[5] = Restrict;
Field[5].IField = 4;
Field[5].VolumesList = {shell_vol};

Background Field = 5;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------
Physical Surface(101) = Boundary{ Volume{shell_vol}; };
Physical Surface(103) = {fault_main};
Physical Surface(105) = Boundary{ Volume{inner_vol}; };
Physical Volume(10) = {shell_vol};
Physical Volume(11) = {inner_vol};

// -------------------------------------------
// CRITICAL MESH SETTINGS - ALL HEX, PROTECT INNER
// -------------------------------------------

// Step 1: Generate 2D mesh first (creates quads on all surfaces)
Mesh.Algorithm = 8;                  // Delaunay for quads
Mesh.RecombineAll = 1;               // Recombine 2D to quads

// Step 2: 3D meshing ONLY for outer volume
// The inner volume ALREADY has its structured hex from Extrude+Recombine
Mesh.Algorithm3D = 1;                // Delaunay 3D (makes tets initially)
Mesh.Recombine3DAll = 1;             // Then recombine tets to hex

// Quality and optimization
Mesh.Smoothing = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

// Output format
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll = 1;

Printf("========================================");
Printf("ALL-HEX MESH");
Printf("Inner: %g x %g x %g structured HEX cubes", nx, ny, nz);
Printf("Outer: Adaptive unstructured HEX");
Printf("NO PYRAMIDS");
Printf("========================================");
