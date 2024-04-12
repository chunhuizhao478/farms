SetFactory("OpenCASCADE");

// Define the sizes of the boxes
Lx1 = 400; // Length of outer box in the x-direction
Ly1 = 150; // Width of outer box in the y-direction
Lz1 = 100; // Height of outer box in the z-direction

Lx2 = 100; // Length of inner box in the x-direction
Ly2 = 40;  // Width of inner box in the y-direction
Lz2 = 0.2; // Height of inner box in the z-direction

// Define the mesh sizes
OuterMeshSize = 20;
InnerMeshSize = 0.1;

// Create the outer box
OuterBox() = newv;
Box(OuterBox) = {-Lx1/2, -Ly1/2, -Lz1/2, Lx1, Ly1, Lz1};

// Create the inner box
InnerBox() = newv;
Box(InnerBox) = {-Lx2/2, -Ly2/2, -Lz2/2, Lx2, Ly2, Lz2};

// Define physical volumes
Physical Volume("OuterBox") = {OuterBox};
Physical Volume("InnerBox") = {InnerBox};

// Mesh size fields
Field[1] = Box;
Field[1].VIn = InnerMeshSize;
Field[1].VOut = OuterMeshSize;
Field[1].XMin = -Lx2/2; Field[1].XMax = Lx2/2;
Field[1].YMin = -Ly2/2; Field[1].YMax = Ly2/2;
Field[1].ZMin = -Lz2/2; Field[1].ZMax = Lz2/2;

// Apply the mesh size field
Field[2] = Min;
Field[2].FieldsList = {1};
Background Field = 2;

// Optional: Refine mesh at boundaries for smoother gradient
Field[3] = BoundaryLayer;
Field[3].EdgesList = {OuterBox, InnerBox};
Field[3].hfar = OuterMeshSize;
Field[3].hwall_n = InnerMeshSize;
Field[3].thickness = 10;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

Mesh.Algorithm = 5;