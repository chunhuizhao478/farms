SetFactory("OpenCASCADE");

// your two target sizes
lc_near = 2.5e-4;   // fine at borehole
lc_far  =   1e-3;   // coarse at outer edge
R_transition = 0.1; // 5 cm from borehole→full coarse

//-- geometry -------------------------------------------------------
Rectangle(1) = {0, 0, 0, 0.054, 0.105, 0};
Rectangle(2) = {0, 0, 0, 0.054, 0.005, 0};
Rectangle(3) = {0, 0.100, 0, 0.054, 0.005, 0};
Disk(4) = {0.027, 0.0525, 0, 0.0055, 0.0055};
BooleanDifference { Surface{1,2,3}; Delete; } { Surface{4}; Delete; }

//Point
Point(101) = {0.027, 0.0525, 0, 0}; // borehole center

//-- mesh‐size field ------------------------------------------------
// 1) distance from borehole curve (curve 5)
Field[1] = Distance;
Field[1].PointsList = {101};

// 2) threshold: size vs. distance
Field[2] = Threshold;
Field[2].IField    = 1;
Field[2].LcMin     = lc_near;     // at dist = 0
Field[2].LcMax     = lc_far;      // at dist ≥ R_transition
Field[2].DistMin   = 0;
Field[2].DistMax   = R_transition;

// 3) use that as the only mesh‐size prescription
Background Field = 2;

//-- physical tags --------------------------------------------------
Physical Curve("borehole")            = {5};
Physical Curve("top")                 = {10};
Physical Curve("bottom")              = {8};
Physical Curve("left")                = {1,6,9};
Physical Curve("right")               = {3,7,11};
Physical Surface("top_elastic_region")    = {3};
Physical Surface("bottom_elastic_region") = {2};
Physical Surface("damagable_domain")      = {4};