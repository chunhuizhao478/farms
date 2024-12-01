SetFactory("OpenCASCADE");

lc = 2e3;
lc_fault = 100; //change this!!!

// Define parameters for the big box
xmin_big = -15000; //-60000
xmax_big = 15000; //60000
ymin_big = -10000; //-60000
ymax_big = 10000; //60000
zmin_big = -20000; //-60000
zmax_big = 0;

// Define dimensions for the big box
dx_big = xmax_big - xmin_big;
dy_big = ymax_big - ymin_big;
dz_big = zmax_big - zmin_big;

// Define parameters for the damagable box
xmin_d = -10000; //-15000
xmax_d = 10000; //15000
ymin_d = -2000;
ymax_d = 2000;
zmin_d = -15000; //-15000
zmax_d = 0;

// Define dimensions for the damagable box
dx_d = xmax_d - xmin_d;
dy_d = ymax_d - ymin_d;
dz_d = zmax_d - zmin_d;

// Define parameters for the small box
xmin_small = -8000;
xmax_small = 8000;
ymin_small = -1000;
ymax_small = 1000;
zmin_small = -12000;
zmax_small = 0;

// Define dimensions for the small box
dx_small = xmax_small - xmin_small;
dy_small = ymax_small - ymin_small;
dz_small = zmax_small - zmin_small;

// Create boxes
Box(1) = {xmin_big, ymin_big, zmin_big, dx_big, dy_big, dz_big};
Box(2) = {xmin_d, ymin_d, zmin_d, dx_d, dy_d, dz_d};
Box(3) = {xmin_small, ymin_small, zmin_small, dx_small, dy_small, dz_small};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2,3}; Delete; }{}

// Get resulting volumes
volumes[] = Volume{:};

// Print volume info for debugging
Printf("Available volumes after fragmentation:");
For i In {0:#volumes[]-1}
    Printf("Volume %g", volumes[i]);
EndFor

// Assign physical volumes using correct indices
Physical Volume("InnerBox") = volumes[0];  // First fragment
Physical Volume("DamagableBox") = volumes[1];  // Second fragment
Physical Volume("OuterBox") = volumes[2];  // Second fragment

// Field definition remains the same
Field[1] = Box;
Field[1].VIn = lc_fault;
Field[1].VOut = lc;
Field[1].XMin = xmin_d;
Field[1].XMax = xmax_d;
Field[1].YMin = ymin_d;
Field[1].YMax = ymax_d;
Field[1].ZMin = zmin_d;
Field[1].ZMax = zmax_d;

Background Field = 1;