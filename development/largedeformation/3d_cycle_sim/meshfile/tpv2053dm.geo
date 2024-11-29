SetFactory("OpenCASCADE");

lc = 2e3;
lc_fault = 100; //change this!!!

// Define parameters for the big box
xmin_big = -15000; //-60000
xmax_big = 15000; //60000
ymin_big = -15000; //-60000
ymax_big = 15000; //60000
zmin_big = -15000; //-60000
zmax_big = 0;

// Define dimensions for the big box
dx_big = xmax_big - xmin_big;
dy_big = ymax_big - ymin_big;
dz_big = zmax_big - zmin_big;

// Define parameters for the small box
xmin_small = -7500; //-15000
xmax_small = 7500; //15000
ymin_small = -2000;
ymax_small = 2000;
zmin_small = -7500; //-15000
zmax_small = 0;

// Define dimensions for the small box
dx_small = xmax_small - xmin_small;
dy_small = ymax_small - ymin_small;
dz_small = zmax_small - zmin_small;

// Create boxes
Box(1) = {xmin_big, ymin_big, zmin_big, dx_big, dy_big, dz_big};
Box(2) = {xmin_small, ymin_small, zmin_small, dx_small, dy_small, dz_small};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2}; Delete; }{}

// Get resulting volumes
volumes[] = Volume{:};

// Print volume info for debugging
Printf("Available volumes after fragmentation:");
For i In {0:#volumes[]-1}
    Printf("Volume %g", volumes[i]);
EndFor

// Assign physical volumes using correct indices
Physical Volume("InnerBox") = volumes[0];  // First fragment
Physical Volume("OuterBox") = volumes[1];  // Second fragment

// Field definition remains the same
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