SetFactory("OpenCASCADE");

// UNIT: m

// Define mesh sizes
lc_fault = 100; // Finer mesh size inside the fault zone
lc = 1000;     // Coarser mesh size outside

// Define the big box (outer domain)
// Later the region will be block 2, and use linear elastic material
big_xmin = -12000;
big_xmax = 12000;
big_ymin = -10000;
big_ymax = 10000;
big_zmin = -13000;
big_zmax = 0;

// Define the small fault zone box
// The region will be block 3, and use continuum damage breakage material
small_xmin = -6000;
small_xmax = 6000;
small_ymin = -1000;
small_ymax = 1000;
small_zmin = -10000;
small_zmax = -3000;

// Define the inner damage zone box
// The region will be block 1, and use continuum damage breakage material   
damage_xmin = -5000;
damage_xmax = 5000;
damage_ymin = -100;
damage_ymax = 100;
damage_zmin = -9000;
damage_zmax = -4000;

// Create the big outer box
big_box = newv;
Box(big_box) = {big_xmin, big_ymin, big_zmin, (big_xmax-big_xmin), (big_ymax-big_ymin), (big_zmax-big_zmin)};

// Create the small fault zone box
small_box = newv;
Box(small_box) = {small_xmin, small_ymin, small_zmin, (small_xmax-small_xmin), (small_ymax-small_ymin), (small_zmax-small_zmin)};

// Create the inner damage zone box
damage_box = newv;
Box(damage_box) = {damage_xmin, damage_ymin, damage_zmin, (damage_xmax-damage_xmin), (damage_ymax-damage_ymin), (damage_zmax-damage_zmin)};

// Boolean fragment to properly embed both fault zone and damage zone inside the big box
BooleanFragments{ Volume{big_box,small_box,damage_box}; Delete; }{}

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;  // Finer mesh inside the fault zone
Field[1].VOut = lc;       // Coarser mesh outside
Field[1].XMin = small_xmin;
Field[1].XMax = small_xmax;
Field[1].YMin = small_ymin;
Field[1].YMax = small_ymax;
Field[1].ZMin = small_zmin;
Field[1].ZMax = small_zmax;

// Set the background mesh size
Background Field = 1;

// Assign Physical Volumes
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// Print the number of volumes created
Printf("Number of volumes created: %g", #volumes[]);

// Mesh the geometry
// Mesh 3;

// Optional: Elevate the mesh order
// SetOrder 2;
// OptimizeMesh "HighOrder";
