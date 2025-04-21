SetFactory("OpenCASCADE");

// your size parameter
lc = 0.001;
       
// force everything to be exactly lc
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// 1) rectangle
Rectangle(1) = {0, 0, 0, 0.054, 0.095, 0};

// 2) circular hole
radius_inner = 0.0055;
Disk(2) = {0.027, 0.0475, 0, radius_inner, radius_inner};

// 3) subtract the disk from the rectangle
BooleanDifference { Surface{1}; Delete; } { Surface{2}; Delete; }

// 4) physical groups
// note: after the BooleanDifference the resulting surface is 3
Physical Curve("borehole") = {5};
Physical Curve("top")      = {9};
Physical Curve("bottom")   = {6};
Physical Curve("left")     = {7};
Physical Curve("right")    = {8};
Physical Surface("domain") = {3};