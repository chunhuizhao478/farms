SetFactory("OpenCASCADE");

// your size parameter
lc = 1e-3;
       
// force everything to be exactly lc
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// 1) rectangle
Rectangle(1) = {0, 0, 0, 0.054, 0.105, 0};
Rectangle(2) = {0, 0, 0, 0.054, 0.005, 0};
Rectangle(3) = {0, 0.100, 0, 0.054, 0.005, 0};

// 2) circular hole
radius_inner = 0.0055;
Disk(4) = {0.027, 0.0525, 0, radius_inner, radius_inner};

// 3) subtract the disk from the rectangle
BooleanDifference { Surface{1,2,3}; Delete; } { Surface{4}; Delete; }

// 4) physical groups
// note: after the BooleanDifference the resulting surface is 3
Physical Curve("borehole") = {5};
Physical Curve("top")      = {9};
Physical Curve("bottom")   = {6};
Physical Curve("left")     = {7};
Physical Curve("right")    = {8};

Physical Surface("top_elastic_region") = {3};
Physical Surface("bottom_elastic_region") = {2};
Physical Surface("damagable_domain") = {4};