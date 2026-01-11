#!/usr/bin/env python3
"""
Helper script to identify vertical line numbers in the extruded 3D mesh
for three-point bending boundary conditions.

Usage:
    python identify_lines.py mesh_wohole_3d.msh

This script reads the mesh file and identifies vertical lines at:
- Left support: x=0.004, y=0
- Right support: x=0.024, y=0
- Loading point: x=0.014, y=0.008
"""

import sys
import struct
import numpy as np

def read_msh_file(filename):
    """Read Gmsh .msh file and extract curve information"""

    print(f"Reading mesh file: {filename}")

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find $Entities section (contains geometric entity information)
    entities_start = None
    entities_end = None

    for i, line in enumerate(lines):
        if line.strip() == '$Entities':
            entities_start = i + 1
        elif line.strip() == '$EndEntities':
            entities_end = i
            break

    if entities_start is None:
        print("Error: Could not find $Entities section in mesh file")
        return

    # Parse entities section
    # Format: numPoints numCurves numSurfaces numVolumes
    entity_counts = list(map(int, lines[entities_start].split()))
    num_points, num_curves, num_surfaces, num_volumes = entity_counts

    print(f"\nEntity counts:")
    print(f"  Points: {num_points}")
    print(f"  Curves: {num_curves}")
    print(f"  Surfaces: {num_surfaces}")
    print(f"  Volumes: {num_volumes}")

    # Parse curve information
    # Format: curveTag minX minY minZ maxX maxY maxZ numPhysicalTags [physicalTags] numBoundingPoints [boundingPoints]

    curve_start = entities_start + 1 + num_points

    print(f"\n{'='*80}")
    print("VERTICAL CURVES (z direction):")
    print(f"{'='*80}")
    print(f"{'Curve ID':<10} {'X':<10} {'Y':<10} {'Z-range':<15} {'Location':<25}")
    print(f"{'-'*80}")

    target_curves = {}
    extrude_z = 0.008
    tolerance = 1e-6

    for i in range(num_curves):
        line = lines[curve_start + i].strip().split()
        curve_id = int(line[0])
        minX, minY, minZ = float(line[1]), float(line[2]), float(line[3])
        maxX, maxY, maxZ = float(line[4]), float(line[5]), float(line[6])

        # Check if this is a vertical line (same X and Y, different Z)
        if (abs(minX - maxX) < tolerance and
            abs(minY - maxY) < tolerance and
            abs(maxZ - minZ - extrude_z) < tolerance):

            x_coord = minX
            y_coord = minY

            # Identify specific locations
            location = ""
            if abs(x_coord - 0.004) < tolerance and abs(y_coord - 0.0) < tolerance:
                location = "LEFT SUPPORT ✓"
                target_curves['left_support'] = curve_id
            elif abs(x_coord - 0.024) < tolerance and abs(y_coord - 0.0) < tolerance:
                location = "RIGHT SUPPORT ✓"
                target_curves['right_support'] = curve_id
            elif abs(x_coord - 0.014) < tolerance and abs(y_coord - 0.008) < tolerance:
                location = "LOADING POINT ✓"
                target_curves['loading_point'] = curve_id

            z_range = f"[{minZ:.4f}, {maxZ:.4f}]"
            print(f"{curve_id:<10} {x_coord:<10.4f} {y_coord:<10.4f} {z_range:<15} {location:<25}")

    print(f"{'='*80}\n")

    # Print summary
    if target_curves:
        print("IDENTIFIED BOUNDARY CONDITION LINES:")
        print(f"{'-'*80}")
        if 'left_support' in target_curves:
            print(f"Left Support (x=0.004, y=0):        Curve {target_curves['left_support']}")
        if 'right_support' in target_curves:
            print(f"Right Support (x=0.024, y=0):       Curve {target_curves['right_support']}")
        if 'loading_point' in target_curves:
            print(f"Loading Point (x=0.014, y=0.008):  Curve {target_curves['loading_point']}")

        print(f"\n{'-'*80}")
        print("Add these lines to your .geo file:")
        print(f"{'-'*80}")
        if 'left_support' in target_curves:
            print(f'Physical Curve("left_support") = {{{target_curves["left_support"]}}};')
        if 'right_support' in target_curves:
            print(f'Physical Curve("right_support") = {{{target_curves["right_support"]}}};')
        if 'loading_point' in target_curves:
            print(f'Physical Curve("loading_point") = {{{target_curves["loading_point"]}}};')
        print(f"{'='*80}\n")
    else:
        print("Warning: Could not identify all target curves!")
        print("Check the coordinates and tolerance settings.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python identify_lines.py <mesh_file.msh>")
        sys.exit(1)

    mesh_file = sys.argv[1]
    read_msh_file(mesh_file)
