SetFactory("OpenCASCADE");

// UNIT: m

// Define mesh sizes
lc_fault = 50; // Finer mesh size inside the fault zone
lc = 1000;     // Coarser mesh size outside

// Define the big box (outer domain)
// Later the region will be block 2, and use linear elastic material
big_xmin = -30000;
big_xmax = 30000;
big_ymin = -30000;
big_ymax = 30000;
big_zmin = -30000;
big_zmax = 0;

// Define the small fault zone box
// The region will be block 3, and use continuum damage breakage material
small_xmin = -2000;
small_xmax = 2000;
small_ymin = -1000;
small_ymax = 1000;
small_zmin = -10000;
small_zmax = 0;

// Define the inner damage zone box
// The region will be block 1, and use continuum damage breakage material   
damage_xmin = -1000;
damage_xmax = 1000;
damage_ymin = -100;
damage_ymax = 100;
damage_zmin = -9000;
damage_zmax = 0;

// Create the big outer box
big_box = newv;
Box(big_box) = {big_xmin, big_ymin, big_zmin, (big_xmax - big_xmin), (big_ymax - big_ymin), (big_zmax - big_zmin)};

// Create the small fault zone box
small_box = newv;
Box(small_box) = {small_xmin, small_ymin, small_zmin, (small_xmax - small_xmin), (small_ymax - small_ymin), (small_zmax - small_zmin)};

// Create the inner damage zone box
damage_box = newv;
Box(damage_box) = {damage_xmin, damage_ymin, damage_zmin, (damage_xmax - damage_xmin), (damage_ymax - damage_ymin), (damage_zmax - damage_zmin)};

// Boolean fragment to properly embed both fault zone and damage zone inside the big box
frag_volumes[] = BooleanFragments{ Volume{big_box, small_box, damage_box}; Delete; }{};

// 1. Create a Distance field from points that define the refined region
Field[1] = Distance;
Field[1].SurfacesList = {13,14,15,16,17,18};

// 2. Create a Threshold field that smoothly transitions the mesh size
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc_fault;
Field[2].LcMax = lc;
Field[2].DistMin = 1000;    // Adjust as needed
Field[2].DistMax = 5000;  // Adjust as needed

// Set the Threshold field as the background field
Background Field = 2;

// Assign Physical Volumes using the fragments created by the Boolean operation
For i In {0:#frag_volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {frag_volumes[i]};
EndFor

Printf("Number of volumes created: %g", #frag_volumes[]);

// Mesh the geometry
// Mesh 3;

