SetFactory("OpenCASCADE");

// -----------------------
// Parameters (edit these)
// -----------------------
Lx = 1.0e5;     // half-width in x (m) => domain x ∈ [-Lx, +Lx]
Lz = 4.0e4;     // depth in y (m)     => domain y ∈ [0, Lz]
Nx_half = 200;  // horizontal divisions per half (left/right)
Nz      = 160;  // vertical divisions

// -----------------------
// Points (2D: z=0 here)
// Coordinates: (x, y, 0)
// -----------------------
p1 = newp; Point(p1) = {-Lx,   0.0, 0};
p2 = newp; Point(p2) = {  0.0,  0.0, 0};
p3 = newp; Point(p3) = { Lx,    0.0, 0};
p4 = newp; Point(p4) = {-Lx,   Lz,   0};
p5 = newp; Point(p5) = {  0.0,  Lz,   0};
p6 = newp; Point(p6) = { Lx,    Lz,   0};

// -----------------------
// Lines (bottom, top, sides, fault)
// -----------------------
l1 = newl; Line(l1) = {p1, p2};      // bottom (left half)
l2 = newl; Line(l2) = {p2, p3};      // bottom (right half)
l3 = newl; Line(l3) = {p1, p4};      // left boundary
l4 = newl; Line(l4) = {p3, p6};      // right boundary
l5 = newl; Line(l5) = {p4, p5};      // top (left half)
l6 = newl; Line(l6) = {p5, p6};      // top (right half)
l7 = newl; Line(l7) = {p2, p5};      // vertical fault (internal)

// -----------------------
// Line loops & surfaces
// (counter-clockwise orientation)
// -----------------------
ll_left  = newll; Line Loop(ll_left)  = { l1,  l7, -l5, -l3 };
ll_right = newll; Line Loop(ll_right) = { l2,  l4,  l6, -l7 };

s_left  = news; Plane Surface(s_left)  = { ll_left  };
s_right = news; Plane Surface(s_right) = { ll_right };

// -----------------------
// Structured meshing (optional), quads
// -----------------------
Transfinite Line{ l3, l7, l4 } = Nz + 1;                // vertical lines
Transfinite Line{ l1, l5 }     = Nx_half + 1;           // left bottom/top
Transfinite Line{ l2, l6 }     = Nx_half + 1;           // right bottom/top

Transfinite Surface{ s_left, s_right };
Recombine Surface{ s_left, s_right };                   // comment out if you prefer triangles

// -----------------------
// Physical groups (names read by MOOSE/libMesh)
// -----------------------
Physical Surface("block_1") = { s_left  };              // left block (ID 1 in input concept)
Physical Surface("block_2") = { s_right };              // right block (ID 2)

Physical Curve("left")   = { l3 };
Physical Curve("right")  = { l4 };
Physical Curve("top")    = { l5, l6 };
Physical Curve("bottom") = { l1, l2 };
Physical Curve("fault")  = { l7 };                      // internal interface (single sideset)