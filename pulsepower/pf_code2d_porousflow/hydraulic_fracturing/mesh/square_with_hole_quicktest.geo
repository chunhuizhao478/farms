// Gmsh geo file: 2D square with central circular hole, horizontal fracture on right side, and branch fractures
// Square: 20m x 20m, Hole radius: 0.5m
// Main fracture extends horizontally from the borehole (RIGHT SIDE ONLY)
// Right branch fractures form a straight line at 45 degrees through the main fracture
// NOTE: Left main fracture and left branches have been removed

// Mesh size parameters
lc_outer = 0.5;   // Mesh size at outer boundary
lc_hole = 0.01;    // Mesh size at hole boundary
lc_square = 0.01;  // Mesh size at fracture surface boundaries

// Geometry parameters
L = 20.0;                  // Square side length
R = 0.5;                   // Hole radius

// Main fracture surface parameters (user can modify these)
square_length = 5.0;       // Length of main fracture (horizontal, x-direction)
square_thickness = 0.05;   // Thickness of main fracture (vertical, y-direction)

// Branch fracture parameters (user can modify these)
branch_length = 1.0;              // Length of branch fractures
right_branch_distance = 2.0;      // Distance from borehole to right branches (along main fracture)
right_branch_angle = Pi/4;        // Angle of right branches from vertical (45 degrees)

// Propagation zone parameters (user can modify these)
propagation_length = 3.0;         // Length of refined mesh zone beyond fracture tips
propagation_width = 1.0;          // Width of the propagation zone (perpendicular to propagation direction)

// Derived parameters
y_half = square_thickness / 2;       // Half-thickness in y direction
x_outer_right = R + square_length;   // Outer x-coordinate of right main fracture

// Calculate x-coordinate on circle at y = ±y_half
x_arc = Sqrt(R*R - y_half*y_half);

// Branch positions (right only)
x_right_branch = R + right_branch_distance;     // x-coordinate of right branch center
branch_half_thickness = square_thickness / 2;   // Half-thickness of branches (same as main)

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

// Points on hole boundary (top, left, and bottom)
Point(7) = {0, R, 0, lc_hole};     // Top (90°)
Point(8) = {-R, 0, 0, lc_hole};    // Left (180°)
Point(9) = {0, -R, 0, lc_hole};    // Bottom (270°)

// Points on hole boundary where right main fracture attaches
Point(10) = {x_arc, y_half, 0, lc_square};    // Right fracture top-inner
Point(11) = {x_arc, -y_half, 0, lc_square};   // Right fracture bottom-inner

// Outer corner points of right main fracture
Point(14) = {x_outer_right, y_half, 0, lc_square};   // Right fracture top-outer
Point(15) = {x_outer_right, -y_half, 0, lc_square};  // Right fracture bottom-outer

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
// Borehole boundary arcs (NOT shared with fracture)
Circle(5) = {10, 5, 7};   // Top-right quadrant: 10 -> 7
Circle(6) = {7, 5, 8};    // Top-left quadrant: 7 -> 8
Circle(7) = {8, 5, 9};    // Bottom-left quadrant: 8 -> 9
Circle(8) = {9, 5, 11};   // Bottom-right quadrant: 9 -> 11

// Fracture-borehole shared boundary arc (right side only)
Circle(9) = {11, 5, 10};  // Right fracture inner boundary: 11 -> 10

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

// Inner boundary of main domain (counterclockwise around hole and right fracture)
// Path: 10 -> 7 -> 8 -> 9 -> 11 -> [right fracture outer with branches] -> 10
Curve Loop(2) = {5, 6, 7, 8,
                 -32, -38, -37, -36, -30, -29, -28, 33, 34, 35, -26};

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
Physical Curve("hole_fracture") = {9};       // Curved boundary shared by right fracture and borehole
Physical Curve("hole") = {5, 6, 7, 8};       // Rest of borehole (top-right, top-left, bottom-left, bottom-right arcs)

// Fracture outer edge
Physical Curve("right_fracture_outer") = {29};

// Branch outer edges (right branches only)
Physical Curve("branch_tops") = {34};           // Top edge of upper branch
Physical Curve("branch_bottoms") = {37};        // Bottom edge of lower branch
Physical Curve("branch_sides") = {33, 35, 36, 38};  // Side edges of right branches

// Physical surfaces
Physical Surface("domain") = {1};
Physical Surface("main_fractures") = {5};              // Right main fracture only
Physical Surface("branch_fractures") = {6, 7};         // Right branch fractures only

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

// --- Right upper branch tip propagation zone (inclined) ---
Field[2] = Box;
Field[2].VIn = lc_square;
Field[2].VOut = lc_outer;
Field[2].XMin = x_right_branch + right_branch_dx - propagation_width/4;
Field[2].XMax = x_right_branch + right_branch_dx + propagation_width/4;
Field[2].YMin = y_half + right_branch_dy;
Field[2].YMax = y_half + right_branch_dy + propagation_length/2;
Field[2].ZMin = -1;
Field[2].ZMax = 1;
Field[2].Thickness = 4*propagation_width;

// --- Right lower branch tip propagation zone (inclined to the left) ---
Field[3] = Box;
Field[3].VIn = lc_square;
Field[3].VOut = lc_outer;
Field[3].XMin = x_right_branch - right_branch_base_shift - right_branch_dx - propagation_width/4;
Field[3].XMax = x_right_branch - right_branch_base_shift - right_branch_dx + propagation_width/4;
Field[3].YMin = -y_half - right_branch_dy - propagation_length/2;
Field[3].YMax = -y_half - right_branch_dy;
Field[3].ZMin = -1;
Field[3].ZMax = 1;
Field[3].Thickness = 4*propagation_width;

// --- Combine all propagation zone fields ---
Field[4] = Min;
Field[4].FieldsList = {1, 2, 3};

// --- Set as background mesh field ---
Background Field = 4;
