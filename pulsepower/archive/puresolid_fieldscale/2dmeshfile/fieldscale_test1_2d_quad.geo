// // Clean structured quad mesh script
// Use built-in geometry kernel for better structured meshing
// SetFactory("Built-in");

// Parameters
lc0 = 0.0001;  // Finest mesh near hole
lc = 0.0001;   // Intermediate mesh
lc2 = 0.0005;  // Coarsest mesh at outer boundary

height = 0.0005;
radius_outer = 0.0125;
radius_refinedbc = 0.0025;
radius_inner = 0.0020;

// Number of divisions
nr1 = 6;   // Radial divisions: inner to middle
nr2 = 20;  // Radial divisions: middle to outer
nc = 124;   // Circumferential divisions (must be divisible by 4)

// Create center point
Point(1) = {0, 0, 0, lc0};

// Create points for inner circle (4 points + center)
Point(2) = {radius_inner, 0, 0, lc0};
Point(3) = {0, radius_inner, 0, lc0};
Point(4) = {-radius_inner, 0, 0, lc0};
Point(5) = {0, -radius_inner, 0, lc0};

// Create points for middle circle
Point(6) = {radius_refinedbc, 0, 0, lc};
Point(7) = {0, radius_refinedbc, 0, lc};
Point(8) = {-radius_refinedbc, 0, 0, lc};
Point(9) = {0, -radius_refinedbc, 0, lc};

// Create points for outer circle
Point(10) = {radius_outer, 0, 0, lc2};
Point(11) = {0, radius_outer, 0, lc2};
Point(12) = {-radius_outer, 0, 0, lc2};
Point(13) = {0, -radius_outer, 0, lc2};

// Create inner circle arcs (quarter circles)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Create middle circle arcs
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Create outer circle arcs
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};

// Create radial lines
Line(13) = {2, 6};  Line(14) = {6, 10};  // Right
Line(15) = {3, 7};  Line(16) = {7, 11};  // Top
Line(17) = {4, 8};  Line(18) = {8, 12};  // Left
Line(19) = {5, 9};  Line(20) = {9, 13};  // Bottom

// Create structured surfaces (8 surfaces total)
// Inner ring surfaces
Curve Loop(1) = {1, 15, -5, -13};   Plane Surface(1) = {1};   // Right-top
Curve Loop(2) = {2, 17, -6, -15};   Plane Surface(2) = {2};   // Left-top
Curve Loop(3) = {3, 19, -7, -17};   Plane Surface(3) = {3};   // Left-bottom
Curve Loop(4) = {4, 13, -8, -19};   Plane Surface(4) = {4};   // Right-bottom

// Outer ring surfaces
Curve Loop(5) = {5, 16, -9, -14};   Plane Surface(5) = {5};   // Right-top
Curve Loop(6) = {6, 18, -10, -16};  Plane Surface(6) = {6};   // Left-top
Curve Loop(7) = {7, 20, -11, -18};  Plane Surface(7) = {7};   // Left-bottom
Curve Loop(8) = {8, 14, -12, -20};  Plane Surface(8) = {8};   // Right-bottom

// Set transfinite properties
// Circumferential divisions
Transfinite Curve {1, 2, 3, 4} = nc/4 Using Progression 1.0;
Transfinite Curve {5, 6, 7, 8} = nc/4 Using Progression 1.0;
Transfinite Curve {9, 10, 11, 12} = nc/4 Using Progression 1.0;

// Radial divisions with progression for gradual mesh refinement
Transfinite Curve {13, 15, 17, 19} = nr1 Using Progression 1.1;  // Inner to middle
Transfinite Curve {14, 16, 18, 20} = nr2 Using Progression 1.2;  // Middle to outer

// Make all surfaces transfinite and recombine to quads
Transfinite Surface {1, 2, 3, 4, 5, 6, 7, 8};
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8};

// Mesh algorithm settings for quads
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;        // Frontal-Delaunay for quads
Mesh.RecombineAll = 1;     // Recombine triangles into quads
Mesh.SubdivisionAlgorithm = 1; // Quad subdivision

// Physical groups
Physical Curve("Inner_Hole") = {1, 2, 3, 4};
Physical Curve("Outer_Boundary") = {9, 10, 11, 12};
Physical Surface("Domain") = {1, 2, 3, 4, 5, 6, 7, 8};

// Optional: For 3D extrusion
/*
nz = 5;  // Number of layers in z-direction
Extrude {0, 0, height} {
  Surface{1, 2, 3, 4, 5, 6, 7, 8}; 
  Layers{nz}; 
  Recombine;
}
Physical Volume("Domain_3D") = {1, 2, 3, 4, 5, 6, 7, 8};
*/