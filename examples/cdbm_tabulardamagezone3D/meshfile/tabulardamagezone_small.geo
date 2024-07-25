SetFactory("OpenCASCADE");

lc1 = 50;
lc2 = 0.5;

//length of large box
Lbox_length_x = 200;
Lbox_length_y = 200;
Lbox_length_z = 200;

//length of small box
Sbox_length_x = 60; 
Sbox_length_y = 10;
Sbox_length_z = 1;

// Define the larger box (Box 2)
Box(1) = {-Lbox_length_x/2, -Lbox_length_y/2, -Lbox_length_z/2, Lbox_length_x, Lbox_length_y, Lbox_length_z};

// Define the smaller box (Box 1) and length it inside the larger box
Box(2) = {-Sbox_length_x/2, -Sbox_length_y/2, -Sbox_length_z/2, Sbox_length_x, Sbox_length_y, Sbox_length_z};

// Embed the smaller box inside the larger box
BooleanFragments { Volume{1}; Delete; } { Volume{2}; Delete; }

// Define mesh sizes using Box field for smaller box
Field[1] = Box;
Field[1].VIn = lc2;
Field[1].VOut = lc1;
Field[1].XMin = -Sbox_length_x/2;
Field[1].XMax = Sbox_length_x/2;
Field[1].YMin = -Sbox_length_y/2;
Field[1].YMax = Sbox_length_y/2;
Field[1].ZMin = -Sbox_length_z/2;
Field[1].ZMax = Sbox_length_z/2;
Field[1].Thickness = 10;

Background Field = 1;

// // Mark the smaller box as a physical volume
// Physical Volume("SmallerBox") = {2};

// // If needed, mark the larger box as a physical volume
// Physical Volume("LargerBox") = {1};