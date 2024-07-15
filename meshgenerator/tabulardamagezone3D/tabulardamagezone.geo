SetFactory("OpenCASCADE");

lc1 = 50;
lc2 = 0.5;

// Define the larger box (Box 2)
Box(1) = {0, 0, 0, 400, 150, 100};

// Define the smaller box (Box 1) and center it inside the larger box
Box(2) = {(400-200)/2, (150-40)/2, (100-0.6)/2, 200, 40, 1};

// Define mesh sizes
Field[1] = Box;
Field[1].VIn = lc2;
Field[1].VOut = lc1;
Field[1].XMin = (400-200)/2;
Field[1].XMax = (400-200)/2 + 200;
Field[1].YMin = (150-40)/2;
Field[1].YMax = (150-40)/2 + 40;
Field[1].ZMin = (100-0.6)/2;
Field[1].ZMax = (100-0.6)/2 + 0.6;

// Equivalent of propagation size on element
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc2;
Field[2].LcMax = lc1;
Field[2].DistMin = 2*lc2;
Field[2].DistMax = 2*lc2+0.001;

Field[3] = Min;
Field[3].FieldsList = {1,2};

Background Field = 3;