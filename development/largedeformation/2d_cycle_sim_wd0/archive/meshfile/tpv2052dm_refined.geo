SetFactory("OpenCASCADE");

// Define mesh sizes
lc_fault = 25;
lc = 5e3;

// Define the big square (2D)
big_xmin = -60000;
big_xmax = 60000;
big_ymin = -60000;
big_ymax = 60000;

// Define the small box (2D)
small_xmin = -15000;
small_xmax = 15000;
small_ymin = -1000;
small_ymax = 1000;

// Define the initial damage box (2D)
smalld_xmin = -10000;
smalld_xmax = 10000;
smalld_ymin = -25;
smalld_ymax = 25;

// Define the nucleation patch (2D)
nucl_xmin = -200;
nucl_xmax = 200;
nucl_ymin = -25;
nucl_ymax = 25;

// Define points for the big square
Point(1) = {big_xmin, big_ymin, 0, lc};
Point(2) = {big_xmax, big_ymin, 0, lc};
Point(3) = {big_xmax, big_ymax, 0, lc};
Point(4) = {big_xmin, big_ymax, 0, lc};

// Define lines for the big square
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Define points for the small box
Point(5) = {small_xmin, small_ymin, 0, lc_fault};
Point(6) = {small_xmax, small_ymin, 0, lc_fault};
Point(7) = {small_xmax, small_ymax, 0, lc_fault};
Point(8) = {small_xmin, small_ymax, 0, lc_fault};

// Define lines for the small box
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Define points for the initial damage strip
Point(9) = {smalld_xmin, smalld_ymin, 0, lc_fault};
Point(10) = {smalld_xmax, smalld_ymin, 0, lc_fault};
Point(11) = {smalld_xmax, smalld_ymax, 0, lc_fault};
Point(12) = {smalld_xmin, smalld_ymax, 0, lc_fault};

// Define lines for the small box
Line(9) = {9, 10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12, 9};


// Create line loops
Line Loop(1) = {1, 2, 3, 4};  // Big square
Line Loop(2) = {5, 6, 7, 8};  // Small box
Line Loop(3) = {9,10,11,12};  // Initial Damage

// Create surfaces for the big square and small box
Plane Surface(1) = {1};  // Big square surface
Plane Surface(2) = {2};  // Small box surface
Plane Surface(3) = {3};  // Small box surface

// Boolean operation to fragment all surfaces (assuming Surface IDs are 1 and 2)
BooleanFragments{ Surface{1,2,3}; Delete; }{}

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;  // Mesh size inside the fault zone
Field[1].VOut = lc;       // Mesh size outside the fault zone
Field[1].XMin = xmin_small;
Field[1].XMax = xmax_small;
Field[1].YMin = ymin_small;
Field[1].YMax = ymax_small;

Background Field = 1;

// Mark all surfaces as physical surfaces
surfaces[] = Surface{:};  // Collect all existing surfaces
For i In {0:#surfaces[]-1}
    Physical Surface(Sprintf("Surface_%g", i+1)) = {surfaces[i]};
EndFor

// Print the number of surfaces created
Printf("Number of surfaces created: %g", #surfaces[]);

// Mesh the geometry
Mesh 2;


