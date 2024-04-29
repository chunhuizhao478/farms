import gmsh
import numpy as np
import meshio
import random

# Initialize Gmsh API
gmsh.initialize()
gmsh.model.add("mesh_example")

# ser boundary domain for simulation
xmax = 800
xmin = -800
ymax = 800
ymin = -800

# set boundary domain where the faults are allowed
fxmax = 600
fxmin = -600
fymax = 600
fymin = -600

#material properties and loading 
sigmaxx = -135; sigmaxy = 70; sigmayy = -120; mus = 0.677; mud = 0.4

# Define the vertices of the main domain
lc = 200  # Coarse mesh length for the outer boundary
# Define mesh size for refined mesh
lf_refined = 8  # Uniform refined mesh size
# Random range
randint = 50

points = [
    gmsh.model.occ.addPoint(xmin, ymin, 0, lc, tag=1),
    gmsh.model.occ.addPoint(xmax, ymin, 0, lc, tag=2),
    gmsh.model.occ.addPoint(xmax, ymax, 0, lc, tag=3),
    gmsh.model.occ.addPoint(xmin, ymax, 0, lc, tag=4)
]

# Create the main square using a loop
lines = []
for i in range(len(points)):
    lines.append(gmsh.model.occ.addLine(points[i], points[(i + 1) % len(points)]))

# Close the loop and create a surface
loop = gmsh.model.occ.addCurveLoop(lines)
plane_surface = gmsh.model.occ.addPlaneSurface([loop])

# Synchronize to create the initial mesh
gmsh.model.occ.synchronize()

# Mesh the surface with coarse mesh
gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay meshing
gmsh.model.mesh.generate(2)

# Get all the elements in the mesh
element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim=2)
triangles = np.array(node_tags[0]).reshape(-1, 3) if node_tags else []

# Retrieve the node coordinates
_, node_coords, _ = gmsh.model.mesh.getNodes()
node_coords = np.array(node_coords).reshape(-1, 3)  # Reshape to [x, y, z]

#add random number
for i in range(np.shape(node_coords)[0]):

    if (node_coords[i][0] < xmax and node_coords[i][0] > xmin and node_coords[i][1] < ymax and node_coords[i][1] > ymin):
        node_coords[i][0] = node_coords[i][0] + random.SystemRandom().randint(-randint, randint)
        node_coords[i][1] = node_coords[i][1] + random.SystemRandom().randint(-randint, randint)

# Clear the existing mesh to refine
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
gmsh.model.mesh.clear()
gmsh.finalize()

#
def segmentinfobuild(x1,y1,x2,y2,segments_info_array):

    #compute fault line slope
    slope = (y2 - y1) / (x2 - x1)
    radian = np.arctan(slope)
    angle = np.degrees(radian)

    #local normal stress
    local_normal_sts = 0.5 * ( sigmaxx + sigmayy ) - 0.5 * ( sigmaxx - sigmayy ) * np.cos(2*np.radians(angle)) - sigmaxy * np.sin(2*np.radians(angle))

    #local shear stress
    local_shear_sts = -0.5 * ( sigmaxx - sigmayy ) * np.sin(2*np.radians(angle)) + sigmaxy * np.cos(2*np.radians(angle))

    #compute mu
    mu = abs(local_shear_sts / local_normal_sts)

    #compute S
    S = (mus - mu) / (mu - mud);

    #store info
    if segments_info_array is None:
        segments_info_array = [x1, y1, x2, y2, local_shear_sts, local_normal_sts, angle, mu, S]
    else:
        segments_info_array.extend([x1, y1, x2, y2, local_shear_sts, local_normal_sts, angle, mu, S])

    return segments_info_array

#
gmsh.initialize()
gmsh.model.add("mesh_example2")
points = [
    gmsh.model.occ.addPoint(xmin, ymin, 0, lf_refined),
    gmsh.model.occ.addPoint(xmax, ymin, 0, lf_refined),
    gmsh.model.occ.addPoint(xmax, ymax, 0, lf_refined),
    gmsh.model.occ.addPoint(xmin, ymax, 0, lf_refined)
]

# Create the main square using a loop
lines = []
for i in range(len(points)):
    lines.append(gmsh.model.occ.addLine(points[i], points[(i + 1) % len(points)]))

# Close the loop and create a surface
loop = gmsh.model.occ.addCurveLoop(lines)
plane_surface = gmsh.model.occ.addPlaneSurface([loop])

# Synchronize to create the initial mesh
gmsh.model.occ.synchronize()
# Add triangles within the boundaries to create a refined mesh
refined_surfaces = {}
refined_surfaces_pts_coords = {}
# build segment info array
segments_info_array = []
for i, tri in enumerate(triangles):
    pts_indices = tri - 1
    pts_coords = node_coords[pts_indices]   
    
    centroid_x = 1/3 * ( pts_coords[0][0] + pts_coords[1][0] + pts_coords[2][0] )
    centroid_y = 1/3 * ( pts_coords[0][1] + pts_coords[1][1] + pts_coords[2][1] )
    if ( centroid_x < fxmax and centroid_x > fxmin and centroid_y < fymax and centroid_y > fymin ):
    
        pids = [gmsh.model.occ.addPoint(x, y, 0, lf_refined) for x, y, _ in pts_coords]
        lns = [gmsh.model.occ.addLine(pids[j], pids[(j + 1) % 3]) for j in range(3)]
        lloop = gmsh.model.occ.addCurveLoop(lns)
        surf = gmsh.model.occ.addPlaneSurface([lloop])

        refined_surfaces[i] = surf
        refined_surfaces_pts_coords[i] = pts_coords

        #calculate material properties
        xA = pts_coords[0][0]; xB = pts_coords[1][0]; xC = pts_coords[2][0];
        yA = pts_coords[0][1]; yB = pts_coords[1][1]; yC = pts_coords[2][1];

        #
        segments_info_array = segmentinfobuild(xA,yA,xB,yB,segments_info_array)
        segments_info_array = segmentinfobuild(xB,yB,xC,yC,segments_info_array)
        segments_info_array = segmentinfobuild(xC,yC,xA,yA,segments_info_array)

#
print(np.array(segments_info_array))
segments_info_array = np.array(segments_info_array).reshape(-1,9)
np.savetxt("./segment_info_array.txt",segments_info_array, fmt="%.3f")

# Synchronize and regenerate mesh with uniform refined mesh size
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay meshing
gmsh.model.mesh.generate(2)
gmsh.fltk.run()
gmsh.write("./refinedmesh.msh")
# Finalize gmsh session
gmsh.finalize()

#
def triangle_area(x1, y1, x2, y2, x3, y3):
    return abs((x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0)

def point_in_triangle(x1, y1, x2, y2, x3, y3, px, py, tolerance=1e-10):
    # Calculate area of the full triangle
    full_area = triangle_area(x1, y1, x2, y2, x3, y3)
    
    # Calculate area of the triangle PAB
    area1 = triangle_area(px, py, x1, y1, x2, y2)
    
    # Calculate area of the triangle PBC
    area2 = triangle_area(px, py, x2, y2, x3, y3)
    
    # Calculate area of the triangle PCA
    area3 = triangle_area(px, py, x3, y3, x1, y1)
    
    # Check if sum of area1, area2, and area3 is within the tolerance of full_area
    if abs(full_area - (area1 + area2 + area3)) < tolerance:
        return True
    else:
        return False

#
#read file
m = meshio.read("./refinedmesh.msh")

#Get TRIA3 Connectivity
tria_elem_connect = m.cells_dict['triangle']

#TRIA3 Elem Num
num_elem = np.shape(tria_elem_connect)[0]

# Output refined triangle elements contained within each big surface
subdomainid = []
elementid = []

#Loop over elem
for elem_ind in range(num_elem):

    print("elem_ind: ", elem_ind)
    print("Progress: ", elem_ind/num_elem*100,"%")    

    #get connectivity of current element
    elem_connect_i = tria_elem_connect[elem_ind,:]

    #get end points coordinates for current element
    coord_data_x = m.points[elem_connect_i][:,0]
    coord_data_y = m.points[elem_connect_i][:,1]

    #get centroid point coordinate for current element
    coord_centriod_x = np.sum(coord_data_x) / 3
    coord_centriod_y = np.sum(coord_data_y) / 3  

    #loop over surfaces
    for key, data in refined_surfaces_pts_coords.items():
        
        #
        x1 = data[0][0]
        x2 = data[1][0] 
        x3 = data[2][0]
        y1 = data[0][1]
        y2 = data[1][1]
        y3 = data[2][1]

        check = point_in_triangle(x1, y1, x2, y2, x3, y3, coord_centriod_x, coord_centriod_y)
        if check:
            elementid.append(elem_ind)
            subdomainid.append(key+1)
            break
        else:
            continue

blockid = np.unique(subdomainid)

print(subdomainid)
print(elementid)
print(blockid)

np.savetxt("./subdomainid.txt",subdomainid,fmt='%i',newline=" ")
np.savetxt("./elementid.txt",elementid,fmt='%i',newline=" ")
np.savetxt("./blockid.txt",blockid,fmt='%i',newline=" ")
