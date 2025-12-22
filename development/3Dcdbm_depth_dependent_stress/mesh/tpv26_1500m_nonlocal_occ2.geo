/**
 * Simple uniform test mesh for slip weakening friction
 * Small domain with uniform element size throughout
 */

SetFactory("OpenCASCADE");

// Uniform element size everywhere
lc = 500; // 500m uniform mesh

// Small domain for testing
Fault_length = 10e3;   // 10 km fault
Fault_width = 5e3;     // 5 km depth
Fault_dip = 90*Pi/180.; // Vertical fault

// Small box domain
Xmax = 5e3;  // ±5 km in x
Xmin = -Xmax;
Ymin = -5e3; // ±5 km in y (perpendicular to fault)
Ymax = 5e3;
Zmin = -5e3; // 5 km depth

// Inner volume around fault (thin layer)
transition_length = 1.5e3;
X_inner_min = -0.5*Fault_length;
X_inner_max =  0.5*Fault_length;
Y_inner_min = -transition_length;
Y_inner_max =  transition_length;
Z_inner_top = 0;
Z_inner_bot = -Fault_width;

// -------------------------------------------
// GEOMETRY
// -------------------------------------------

// Outer and inner boxes
Box(1) = {Xmin, Ymin, Zmin, Xmax-Xmin, Ymax-Ymin, -Zmin};
Box(2) = {X_inner_min, Y_inner_min, Z_inner_bot,
          X_inner_max-X_inner_min, Y_inner_max-Y_inner_min, Z_inner_top-Z_inner_bot};

// Fault plane (vertical)
Rectangle(1000) = {-0.5*Fault_length, 0, 0, Fault_length, Fault_width};
Rotate{{1, 0, 0}, {0, 0, 0}, -Fault_dip} { Surface{1000}; }

// Fragment volumes by fault
BooleanFragments{ Volume{1,2}; Delete; }{ Surface{1000}; }

// Get resulting volumes
vols[] = Volume{:};
shell_vol = vols[0];
inner_vol = vols[1];

fault_main = 1000;

// -------------------------------------------
// UNIFORM HEX MESH SETTINGS
// -------------------------------------------

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Characteristic Length { PointsOf{ Volume{:}; } } = lc;

// -------------------------------------------
// PHYSICAL GROUPS
// -------------------------------------------

Physical Surface(101) = Boundary{ Volume{shell_vol}; };
Physical Surface(103) = {fault_main};
Physical Surface(105) = Boundary{ Volume{inner_vol}; };

Physical Volume(10) = {shell_vol};
Physical Volume(11) = {inner_vol};

// -------------------------------------------
// UNIFORM HEX8 MESH
// -------------------------------------------

Mesh.RecombineAll = 1;              // Create hexahedra
Mesh.Algorithm = 8;                 // Frontal-Delaunay for Quads  
Mesh.Algorithm3D = 6;               // Frontal hex
Mesh.RecombineOptimizeTopology = 5;
Mesh.Recombine3DAll = 1;
Mesh.MshFileVersion = 2.2;
