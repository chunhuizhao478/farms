SetFactory("OpenCASCADE");

lc1 = 50;
lc2 = 0.5;

// Define the larger box (Box 2)
Box(1) = {0, 0, 0, 400, 150, 100};

// Define the smaller box (Box 1) and center it inside the larger box
Box(2) = {(400-200)/2, (150-40)/2, (100-10)/2, 200, 40, 10};

// Embed the smaller box inside the larger box
BooleanFragments { Volume{1}; Delete; } { Volume{2}; Delete; }

// Define mesh sizes using Box field for smaller box
Field[1] = Box;
Field[1].VIn = lc2;
Field[1].VOut = lc1;
Field[1].XMin = (400-200)/2;
Field[1].XMax = (400-200)/2 + 200;
Field[1].YMin = (150-40)/2;
Field[1].YMax = (150-40)/2 + 40;
Field[1].ZMin = (100-10)/2;
Field[1].ZMax = (100-10)/2 + 10;
Field[1].Thickness = 10;

Background Field = 1;