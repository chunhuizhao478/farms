// Gmsh geo file: 2D square with central circular hole, horizontal fractures, and branch fractures
// Square: 20m x 20m, Hole radius: 0.5m
// Main fractures extend horizontally from the borehole
// Left branch fractures extend vertically from the main fractures
// Right branch fractures form a straight line at 45 degrees through the main fracture

// Mesh size parameters
lc_outer = 0.5;   // Mesh size at outer boundary
lc_hole = 0.01;    // Mesh size at hole boundary
lc_square = 0.01;  // Mesh size at fracture surface boundaries

// Geometry parameters
L = 20.0;                  // Square side length
R = 0.5;                   // Hole radius

// Main fracture surface parameters (user can modify these)
square_length = 5.0;       // Length of main fractures (horizontal, x-direction)
square_thickness = 0.05;   // Thickness of main fractures (vertical, y-direction)

// Branch fracture parameters (user can modify these)
branch_length = 1.0;              // Length of branch fractures
left_branch_distance = 3.0;       // Distance from borehole to left branches (along main fracture)
right_branch_distance = 2.0;      // Distance from borehole to right branches (along main fracture)
right_branch_angle = Pi/4;        // Angle of right branches from vertical (45 degrees)

// Propagation zone parameters (user can modify these)
propagation_length = 3.0;         // Length of refined mesh zone beyond fracture tips
propagation_width = 1.0;          // Width of the propagation zone (perpendicular to propagation direction)

// Derived parameters
y_half = square_thickness / 2;       // Half-thickness in y direction
x_outer_right = R + square_length;   // Outer x-coordinate of right main fracture
x_outer_left = -(R + square_length); // Outer x-coordinate of left main fracture

// Calculate x-coordinate on circle at y = ±y_half
x_arc = Sqrt(R*R - y_half*y_half);

// Branch positions
x_left_branch = -(R + left_branch_distance);    // x-coordinate of left branch center
x_right_branch = R + right_branch_distance;     // x-coordinate of right branch center
branch_half_thickness = square_thickness / 2;   // Half-thickness of branches (same as main)

// Branch vertical extents (for left branches - vertical)
y_branch_top = y_half + branch_length;          // Top of upper branches
y_branch_bottom = -(y_half + branch_length);    // Bottom of lower branches

// Right branch inclined tip offsets (45 degree angle from vertical)
right_branch_dx = branch_length * Sin(right_branch_angle);  // Horizontal offset of tip
right_branch_dy = branch_length * Cos(right_branch_angle);  // Vertical offset of tip

// Horizontal shift for lower right branch base to align with upper branch (forms straight line)
// Accounts for the main fracture thickness
right_branch_base_shift = square_thickness * Tan(right_branch_angle);

// ==================== POINTS ====================

// Square corner points (centered at origin)
Point(1) = {-L/2, -L/2, 0, lc_outer};  // Bottom-left
Point(2) = {L/2, -L/2, 0, lc_outer};   // Bottom-right
Point(3) = {L/2, L/2, 0, lc_outer};    // Top-right
Point(4) = {-L/2, L/2, 0, lc_outer};   // Top-left

// Center point for hole
Point(5) = {0, 0, 0, lc_hole};

// Points on hole boundary (top and bottom only - no mid points on fracture boundaries)
Point(7) = {0, R, 0, lc_hole};     // Top (90°)
Point(9) = {0, -R, 0, lc_hole};    // Bottom (270°)

// Points on hole boundary where main fractures attach
Point(10) = {x_arc, y_half, 0, lc_square};    // Right fracture top-inner
Point(11) = {x_arc, -y_half, 0, lc_square};   // Right fracture bottom-inner
Point(12) = {-x_arc, y_half, 0, lc_square};   // Left fracture top-inner
Point(13) = {-x_arc, -y_half, 0, lc_square};  // Left fracture bottom-inner

// Outer corner points of main fractures
Point(14) = {x_outer_right, y_half, 0, lc_square};   // Right fracture top-outer
Point(15) = {x_outer_right, -y_half, 0, lc_square};  // Right fracture bottom-outer
Point(16) = {x_outer_left, y_half, 0, lc_square};    // Left fracture top-outer
Point(17) = {x_outer_left, -y_half, 0, lc_square};   // Left fracture bottom-outer

// Left branch connection points (on main left fracture edges)
Point(18) = {x_left_branch + branch_half_thickness, y_half, 0, lc_square};   // Top edge, right of upper branch
Point(19) = {x_left_branch - branch_half_thickness, y_half, 0, lc_square};   // Top edge, left of upper branch
Point(20) = {x_left_branch + branch_half_thickness, -y_half, 0, lc_square};  // Bottom edge, right of lower branch
Point(21) = {x_left_branch - branch_half_thickness, -y_half, 0, lc_square};  // Bottom edge, left of lower branch

// Left upper branch outer corners
Point(22) = {x_left_branch + branch_half_thickness, y_branch_top, 0, lc_square};  // Top-right
Point(23) = {x_left_branch - branch_half_thickness, y_branch_top, 0, lc_square};  // Top-left

// Left lower branch outer corners
Point(24) = {x_left_branch + branch_half_thickness, y_branch_bottom, 0, lc_square};  // Bottom-right
Point(25) = {x_left_branch - branch_half_thickness, y_branch_bottom, 0, lc_square};  // Bottom-left

// Right branch connection points (on main right fracture edges)
Point(26) = {x_right_branch - branch_half_thickness, y_half, 0, lc_square};   // Top edge, left of upper branch
Point(27) = {x_right_branch + branch_half_thickness, y_half, 0, lc_square};   // Top edge, right of upper branch
Point(28) = {x_right_branch - branch_half_thickness - right_branch_base_shift, -y_half, 0, lc_square};  // Bottom edge, left of lower branch (shifted for alignment)
Point(29) = {x_right_branch + branch_half_thickness - right_branch_base_shift, -y_half, 0, lc_square};  // Bottom edge, right of lower branch (shifted for alignment)

// Right upper branch outer corners (inclined 45 degrees to the right)
Point(30) = {x_right_branch - branch_half_thickness + right_branch_dx, y_half + right_branch_dy, 0, lc_square};  // Top-left
Point(31) = {x_right_branch + branch_half_thickness + right_branch_dx, y_half + right_branch_dy, 0, lc_square};  // Top-right

// Right lower branch outer corners (inclined 45 degrees to the left, forming straight line with upper branch)
Point(32) = {x_right_branch - branch_half_thickness - right_branch_base_shift - right_branch_dx, -y_half - right_branch_dy, 0, lc_square};  // Bottom-left
Point(33) = {x_right_branch + branch_half_thickness - right_branch_base_shift - right_branch_dx, -y_half - right_branch_dy, 0, lc_square};  // Bottom-right

// ==================== LINES ====================

// Outer square edges
Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Right
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Left

// Circle arcs (all counterclockwise)
// Borehole boundary arcs (NOT shared with fractures - top and bottom portions)
Circle(5) = {10, 5, 7};   // Top-right: 10 -> 7
Circle(6) = {7, 5, 12};   // Top-left: 7 -> 12
Circle(7) = {13, 5, 9};   // Bottom-left: 13 -> 9
Circle(8) = {9, 5, 11};   // Bottom-right: 9 -> 11

// Fracture-borehole shared boundary arcs (single arc for each side, no mid-points)
Circle(9) = {11, 5, 10};  // Right fracture inner boundary: 11 -> 10
Circle(10) = {12, 5, 13}; // Left fracture inner boundary: 12 -> 13

// ----- Left main fracture edges -----
// Top edge (split by upper branch)
Line(13) = {12, 18};   // Inner segment: 12 -> 18
Line(14) = {18, 19};   // Branch connection (shared with upper branch bottom): 18 -> 19
Line(15) = {19, 16};   // Outer segment: 19 -> 16

// Outer (left) edge
Line(16) = {16, 17};   // 16 -> 17

// Bottom edge (split by lower branch)
Line(17) = {17, 21};   // Outer segment: 17 -> 21
Line(18) = {21, 20};   // Branch connection (shared with lower branch top): 21 -> 20
Line(19) = {20, 13};   // Inner segment: 20 -> 13

// ----- Left upper branch edges -----
Line(20) = {18, 22};   // Right edge (going up): 18 -> 22
Line(21) = {22, 23};   // Top edge: 22 -> 23
Line(22) = {23, 19};   // Left edge (going down): 23 -> 19

// ----- Left lower branch edges -----
Line(23) = {21, 25};   // Left edge (going down): 21 -> 25
Line(24) = {25, 24};   // Bottom edge: 25 -> 24
Line(25) = {24, 20};   // Right edge (going up): 24 -> 20

// ----- Right main fracture edges -----
// Top edge (split by upper branch)
Line(26) = {10, 26};   // Inner segment: 10 -> 26
Line(27) = {26, 27};   // Branch connection (shared with upper branch bottom): 26 -> 27
Line(28) = {27, 14};   // Outer segment: 27 -> 14

// Outer (right) edge
Line(29) = {14, 15};   // 14 -> 15

// Bottom edge (split by lower branch)
Line(30) = {15, 29};   // Outer segment: 15 -> 29
Line(31) = {29, 28};   // Branch connection (shared with lower branch top): 29 -> 28
Line(32) = {28, 11};   // Inner segment: 28 -> 11

// ----- Right upper branch edges -----
Line(33) = {27, 31};   // Right edge (going up): 27 -> 31
Line(34) = {31, 30};   // Top edge: 31 -> 30
Line(35) = {30, 26};   // Left edge (going down): 30 -> 26

// ----- Right lower branch edges -----
Line(36) = {29, 33};   // Right edge (going down): 29 -> 33
Line(37) = {33, 32};   // Bottom edge: 33 -> 32
Line(38) = {32, 28};   // Left edge (going up): 32 -> 28

// ==================== CURVE LOOPS ====================

// Outer boundary (counter-clockwise)
Curve Loop(1) = {1, 2, 3, 4};

// Inner boundary of main domain (counterclockwise around all fractures)
// Path: 10 -> 7 -> 12 -> [left fracture outer with branches] -> 13 -> 9 -> 11 -> [right fracture outer with branches] -> 10
Curve Loop(2) = {5, 6,
                 13, 20, 21, 22, 15, 16, 17, 23, 24, 25, 19,
                 7, 8,
                 -32, -38, -37, -36, -30, -29, -28, 33, 34, 35, -26};

// Left main fracture (counterclockwise)
// Path: 12 -> 18 -> 19 -> 16 -> 17 -> 21 -> 20 -> 13 -> [arc] -> 12
Curve Loop(3) = {13, 14, 15, 16, 17, 18, 19, -10};

// Left upper branch (counterclockwise)
// Path: 18 -> 22 -> 23 -> 19 -> 18
Curve Loop(4) = {20, 21, 22, -14};

// Left lower branch (counterclockwise)
// Path: 21 -> 25 -> 24 -> 20 -> 21
Curve Loop(5) = {23, 24, 25, -18};

// Right main fracture (counterclockwise)
// Path: 10 -> 26 -> 27 -> 14 -> 15 -> 29 -> 28 -> 11 -> [arc] -> 10
Curve Loop(6) = {26, 27, 28, 29, 30, 31, 32, 9};

// Right upper branch (counterclockwise)
// Path: 26 -> 27 -> 31 -> 30 -> 26 (uses -27 so shared edge with main fracture is opposite direction)
Curve Loop(7) = {-27, -35, -34, -33};

// Right lower branch (counterclockwise)
// Path: 29 -> 33 -> 32 -> 28 -> 29 (uses -31 so shared edge with main fracture is opposite direction)
Curve Loop(8) = {36, 37, 38, -31};

// ==================== SURFACES ====================

Plane Surface(1) = {1, 2};  // Main domain
Plane Surface(2) = {3};      // Left main fracture
Plane Surface(3) = {4};      // Left upper branch
Plane Surface(4) = {5};      // Left lower branch
Plane Surface(5) = {6};      // Right main fracture
Plane Surface(6) = {7};      // Right upper branch
Plane Surface(7) = {8};      // Right lower branch

// ==================== PHYSICAL GROUPS ====================

// Outer boundary curves
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};

// Borehole boundaries
Physical Curve("hole_fracture") = {9, 10};  // Curved boundaries shared by fractures and borehole
Physical Curve("hole") = {5, 6, 7, 8};      // Rest of borehole (top and bottom arcs)

// Fracture outer edges
Physical Curve("left_fracture_outer") = {16};
Physical Curve("right_fracture_outer") = {29};

// Branch outer edges
Physical Curve("branch_tops") = {21, 34};           // Top edges of upper branches
Physical Curve("branch_bottoms") = {24, 37};        // Bottom edges of lower branches
Physical Curve("branch_sides") = {20, 22, 23, 25, 33, 35, 36, 38};  // Side edges of all branches

// Physical surfaces
Physical Surface("domain") = {1};
Physical Surface("main_fractures") = {2, 5};           // Left and right main fractures
Physical Surface("branch_fractures") = {3, 4, 6, 7};   // All branch fractures

// Mesh settings
Mesh.Algorithm = 6;  // Frontal-Delaunay algorithm

// ==================== MESH REFINEMENT FIELDS ====================
// Extend refined mesh (lc_square) into propagation zones beyond fracture tips

// --- Right main fracture tip propagation zone ---
// Box extending from right fracture tip (x_outer_right) to x_outer_right + propagation_length
Field[1] = Box;
Field[1].VIn = lc_square;
Field[1].VOut = lc_outer;
Field[1].XMin = x_outer_right;
Field[1].XMax = x_outer_right + propagation_length/2;
Field[1].YMin = -propagation_width/4;
Field[1].YMax = propagation_width/4;
Field[1].ZMin = -1;
Field[1].ZMax = 1;
Field[1].Thickness = 4*propagation_width;  // Gradual transition

// --- Left main fracture tip propagation zone ---
// Box extending from left fracture tip (x_outer_left) to x_outer_left - propagation_length
Field[2] = Box;
Field[2].VIn = lc_square;
Field[2].VOut = lc_outer;
Field[2].XMin = x_outer_left - propagation_length/2;
Field[2].XMax = x_outer_left;
Field[2].YMin = -propagation_width/4;
Field[2].YMax = propagation_width/4;
Field[2].ZMin = -1;
Field[2].ZMax = 1;
Field[2].Thickness = 4*propagation_width;

// --- Left upper branch tip propagation zone ---
Field[3] = Box;
Field[3].VIn = lc_square;
Field[3].VOut = lc_outer;
Field[3].XMin = x_left_branch - propagation_width/4;
Field[3].XMax = x_left_branch + propagation_width/4;
Field[3].YMin = y_branch_top;
Field[3].YMax = y_branch_top + propagation_length/2;
Field[3].ZMin = -1;
Field[3].ZMax = 1;
Field[3].Thickness = 4*propagation_width;

// --- Left lower branch tip propagation zone ---
Field[4] = Box;
Field[4].VIn = lc_square;
Field[4].VOut = lc_outer;
Field[4].XMin = x_left_branch - propagation_width/4;
Field[4].XMax = x_left_branch + propagation_width/4;
Field[4].YMin = y_branch_bottom - propagation_length/2;
Field[4].YMax = y_branch_bottom;
Field[4].ZMin = -1;
Field[4].ZMax = 1;
Field[4].Thickness = 4*propagation_width;

// --- Right upper branch tip propagation zone (inclined) ---
Field[5] = Box;
Field[5].VIn = lc_square;
Field[5].VOut = lc_outer;
Field[5].XMin = x_right_branch + right_branch_dx - propagation_width/4;
Field[5].XMax = x_right_branch + right_branch_dx + propagation_width/4;
Field[5].YMin = y_half + right_branch_dy;
Field[5].YMax = y_half + right_branch_dy + propagation_length/2;
Field[5].ZMin = -1;
Field[5].ZMax = 1;
Field[5].Thickness = 4*propagation_width;

// --- Right lower branch tip propagation zone (inclined to the left) ---
Field[6] = Box;
Field[6].VIn = lc_square;
Field[6].VOut = lc_outer;
Field[6].XMin = x_right_branch - right_branch_base_shift - right_branch_dx - propagation_width/4;
Field[6].XMax = x_right_branch - right_branch_base_shift - right_branch_dx + propagation_width/4;
Field[6].YMin = -y_half - right_branch_dy - propagation_length/2;
Field[6].YMax = -y_half - right_branch_dy;
Field[6].ZMin = -1;
Field[6].ZMax = 1;
Field[6].Thickness = 4*propagation_width;

// --- Combine all propagation zone fields ---
Field[7] = Min;
Field[7].FieldsList = {1, 2, 3, 4, 5, 6};

// --- Set as background mesh field ---
Background Field = 7;
