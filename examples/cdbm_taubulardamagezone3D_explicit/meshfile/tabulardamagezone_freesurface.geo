SetFactory("OpenCASCADE");

lc1 = 50; lc2 = 0.05;

//length of large box
Lbox_length_x = 1000; Lbox_length_y = 1000; Lbox_length_z = 1000;

//length of small box
Sbox_length_x = 60;  Sbox_length_y = 10; Sbox_length_z = 0.2;

//length of elastic box
Ebox_length_x = 60;  Ebox_length_y = 10; Ebox_length_z = 50;

// Define the larger box (Box 2)
Box(1) = {-Lbox_length_x/2, -Lbox_length_y/2, -Lbox_length_z/2, Lbox_length_x, Lbox_length_y, Lbox_length_z};

// Define the smaller box (Box 1) and place it inside the larger box
Box(2) = {-Sbox_length_x/2, Lbox_length_y/2-Sbox_length_y, -Sbox_length_z/2, Sbox_length_x, Sbox_length_y, Sbox_length_z};

// Define the nucleation box
center_x = -Sbox_length_x/2+6; center_y = Lbox_length_y/2-Sbox_length_y/2; nucl_radius = 2;
double_nucl_radius = 2 * nucl_radius;
Box(3) = {center_x-nucl_radius,center_y-nucl_radius,-Sbox_length_z/2,double_nucl_radius,double_nucl_radius,Sbox_length_z};

Box(4) = {-Ebox_length_x/2, Lbox_length_y/2-Sbox_length_y, -Sbox_length_z/2-Ebox_length_z, Ebox_length_x, Ebox_length_y, Ebox_length_z};

Box(5) = {-Ebox_length_x/2, Lbox_length_y/2-Sbox_length_y, Sbox_length_z/2, Ebox_length_x, Ebox_length_y, Ebox_length_z};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2,3,4,5}; Delete; }{}

// Define mesh sizes using Box field for smaller box
Field[1] = Box;
Field[1].VIn = lc2;
Field[1].VOut = lc1;
Field[1].XMin = -Sbox_length_x/2;
Field[1].XMax = Sbox_length_x/2;
Field[1].YMin = Lbox_length_y/2-Sbox_length_y;
Field[1].YMax = Lbox_length_y/2;
Field[1].ZMin = -Sbox_length_z/2;
Field[1].ZMax = Sbox_length_z/2;
Field[1].Thickness = 50;

Background Field = 1;

// Mark all volumes as physical volumes
// Label each volume separately
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// Print the number of volumes created
Printf("Number of volumes created: %g", #volumes[]);


