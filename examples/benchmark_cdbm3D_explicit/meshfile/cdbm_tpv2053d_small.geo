SetFactory("OpenCASCADE");

lc = 1000;
lc_fault = 100; // Fine mesh in the fault zone

Fault_length = 10e3;
Fault_width = 5e3;
Fault_thickness = 1e3;
Fault_dip = 90*Pi/180.;

// Nucleation in X,Z local coordinates
X_nucl = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl = 400;
thickness_nucl = 1000;

Xmax = 8e3;
Xmin = -Xmax;

Ymax = 0;
Ymin = -Xmax;

Zmin =  -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Zmax =   Xmax + 0.5 * Fault_width  *Cos(Fault_dip);

Box(1) = {Xmin, 0, Zmin, 2*Xmax, Ymin, 2*Xmax};

// Create a damage zone
Box(2) = {-Fault_length/2, 0-Fault_width, -Fault_thickness/2, Fault_length, Fault_width, Fault_thickness};

// Create a nucleation patch
Box(3) = {X_nucl-R_nucl/2, -Width_nucl-R_nucl/2, -thickness_nucl/2, R_nucl, R_nucl, thickness_nucl};

// Create a damage box
// damage_xdist = 80000;
// damage_ydist = 50000;
// damage_zdist = 40000;
// Box(4) = {-damage_xdist/2, -damage_ydist, -damage_zdist/2, damage_xdist, damage_ydist, damage_zdist};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2,3}; Delete; }{}

// Define mesh sizes using a progression field for smooth transition

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;
Field[1].VOut = lc/2;
Field[1].XMin = -Fault_length/2;
Field[1].XMax = Fault_length/2;
Field[1].YMin = 0-Fault_width;
Field[1].YMax = 0;
Field[1].ZMin = -Fault_thickness/2;
Field[1].ZMax = Fault_thickness/2;

Background Field = 1;

// Mark all volumes as physical volumes
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// Print the number of volumes created
Printf("Number of volumes created: %g", #volumes[]);
