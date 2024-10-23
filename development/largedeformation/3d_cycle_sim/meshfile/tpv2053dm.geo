SetFactory("OpenCASCADE");

lc = 5e3;
lc_fault = 200; //change this!!!

// Define parameters for the big box
xmin_big = -20000; //-60000
xmax_big = 20000; //60000
ymin_big = -20000; //-60000
ymax_big = 20000; //60000
zmin_big = -20000; //-60000
zmax_big = 0;

// Define dimensions for the big box
dx_big = xmax_big - xmin_big;
dy_big = ymax_big - ymin_big;
dz_big = zmax_big - zmin_big;

// Define parameters for the small box
xmin_small = -7500; //-15000
xmax_small = 7500; //15000
ymin_small = -500;
ymax_small = 500;
zmin_small = -7500; //-15000
zmax_small = 0;

// Define dimensions for the small box
dx_small = xmax_small - xmin_small;
dy_small = ymax_small - ymin_small;
dz_small = zmax_small - zmin_small;

// Create the big box
Box(1) = {xmin_big, ymin_big, zmin_big, dx_big, dy_big, dz_big};

// Create the small box
Box(2) = {xmin_small, ymin_small, zmin_small, dx_small, dy_small, dz_small};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2}; Delete; }{}

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;
Field[1].VOut = lc;
Field[1].XMin = xmin_small;
Field[1].XMax = xmax_small;
Field[1].YMin = ymin_small;
Field[1].YMax = ymax_small;
Field[1].ZMin = zmin_small;
Field[1].ZMax = zmax_small;

Background Field = 1;

// Mark all volumes as physical volumes
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// Print the number of volumes created
Printf("Number of volumes created: %g", #volumes[]);