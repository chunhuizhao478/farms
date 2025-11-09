#!/usr/bin/env python3
"""
Parametric Study Case Generator for Pure Solid Simulations

This script generates parametric study cases by creating folders and modified
input files (elasticity.i and fracture.i) based on different combinations of
confinement pressures and domain sizes.

Author: Auto-generated
Date: 2025-11-08
"""

import os
from pathlib import Path

# ==============================================================================
# CONFIGURATION - Modify these parameters as needed
# ==============================================================================

# Base directory (where this script is located)
BASE_DIR = Path(__file__).parent

# Base input files (source templates) - these will not be deleted
BASE_ELASTICITY_FILE = BASE_DIR / "original_files" / "case_name" / "elasticity.i"
BASE_FRACTURE_FILE = BASE_DIR / "original_files" / "case_name" / "fracture.i"

# Parametric study parameters
CONFINEMENT_PRESSURES = {
    "cf1": 1e6,  # 1 MPa
    "cf5": 5e6,  # 5 MPa
    "cf10": 10e6,  # 10 MPa
}

DOMAIN_SIZES = {
    "domain1x": {
        "mesh_file": "../../../2dmeshfile/fieldscale_test1_2d.msh",
        "coord": "0.01 0.01 0",
    },
    "domain2x": {
        "mesh_file": "../../../2dmeshfile/fieldscale_test1_2d_extend2x.msh",
        "coord": "0.02 0.02 0",
    },
    "domain5x": {
        "mesh_file": "../../../2dmeshfile/fieldscale_test1_2d_extend5x.msh",
        "coord": "0.06 0.06 0",
    },
}

# Output folder naming pattern: "case_{cf_key}_{domain_key}"
FOLDER_NAME_PATTERN = "case_{cf_key}_{domain_key}"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================


def modify_elasticity_file(content, confinement_pressure, mesh_file, coord):
    """
    Modify the elasticity.i file content with specified parameters.

    Args:
        content (str): Original file content
        confinement_pressure (float): Confinement pressure value (e.g., 1e6)
        mesh_file (str): Path to mesh file
        coord (str): Coordinate string for fixed point

    Returns:
        str: Modified file content
    """
    lines = content.split("\n")
    modified_lines = []

    for i, line in enumerate(lines):
        # Modify confinement pressure
        if line.strip().startswith("confinement_pressure ="):
            modified_lines.append(f"confinement_pressure = {confinement_pressure}")
        # Modify mesh file path
        elif "file =" in line and "2dmeshfile" in line:
            # Check if this is in the Mesh section by looking at context
            # Find if we're in a FileMeshGenerator block
            context_start = max(0, i - 5)
            context = "\n".join(lines[context_start : i + 1])
            if "type = FileMeshGenerator" in context:
                indent = len(line) - len(line.lstrip())
                modified_lines.append(" " * indent + f"file =  '{mesh_file}'")
            else:
                modified_lines.append(line)
        # Modify fixed point coordinate
        elif "coord =" in line:
            # Check if this is in ExtraNodesetGenerator block
            context_start = max(0, i - 5)
            context_end = min(len(lines), i + 5)
            context = "\n".join(lines[context_start:context_end])
            if "new_boundary = corner_ptr" in context:
                indent = len(line) - len(line.lstrip())
                modified_lines.append(" " * indent + f"coord = '{coord}'")
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    return "\n".join(modified_lines)


def modify_fracture_file(content, mesh_file, coord):
    """
    Modify the fracture.i file content with specified parameters.

    Args:
        content (str): Original file content
        mesh_file (str): Path to mesh file
        coord (str): Coordinate string for fixed point

    Returns:
        str: Modified file content
    """
    lines = content.split("\n")
    modified_lines = []

    for i, line in enumerate(lines):
        # Modify mesh file path
        if "file =" in line and "2dmeshfile" in line:
            indent = len(line) - len(line.lstrip())
            modified_lines.append(" " * indent + f"file = '{mesh_file}'")
        # Modify fixed point coordinate
        elif "coord =" in line:
            # Check if this is in ExtraNodesetGenerator block
            context_start = max(0, i - 5)
            context_end = min(len(lines), i + 5)
            context = "\n".join(lines[context_start:context_end])
            if "new_boundary = corner_ptr" in context:
                indent = len(line) - len(line.lstrip())
                modified_lines.append(" " * indent + f"coord = '{coord}'")
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    return "\n".join(modified_lines)


def create_case_folder(
    case_name, cf_value, domain_config, elasticity_template, fracture_template
):
    """
    Create a case folder with modified input files.

    Args:
        case_name (str): Name of the case folder
        cf_value (float): Confinement pressure value
        domain_config (dict): Domain configuration with 'mesh_file' and 'coord'
        elasticity_template (str): Template content for elasticity.i
        fracture_template (str): Template content for fracture.i
    """
    # Create case directory
    case_dir = BASE_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)

    # Modify files
    modified_elasticity = modify_elasticity_file(
        elasticity_template,
        cf_value,
        domain_config["mesh_file"],
        domain_config["coord"],
    )

    modified_fracture = modify_fracture_file(
        fracture_template, domain_config["mesh_file"], domain_config["coord"]
    )

    # Write modified files to case directory
    with open(case_dir / "elasticity.i", "w") as f:
        f.write(modified_elasticity)

    with open(case_dir / "fracture.i", "w") as f:
        f.write(modified_fracture)

    print(f"  ✓ Created: {case_name}/")


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================


def main():
    """
    Main function to generate all parametric study cases.
    """
    print("=" * 80)
    print("Parametric Study Case Generator")
    print("=" * 80)
    print(f"\nBase directory: {BASE_DIR}")
    print(f"\nTemplate files:")
    print(f"  - Elasticity: {BASE_ELASTICITY_FILE}")
    print(f"  - Fracture:   {BASE_FRACTURE_FILE}")

    # Check if base files exist
    if not BASE_ELASTICITY_FILE.exists():
        print(f"\n✗ ERROR: Base elasticity file not found: {BASE_ELASTICITY_FILE}")
        return

    if not BASE_FRACTURE_FILE.exists():
        print(f"\n✗ ERROR: Base fracture file not found: {BASE_FRACTURE_FILE}")
        return

    # Load template files
    print("\n" + "-" * 80)
    print("Loading template files...")
    print("-" * 80)

    with open(BASE_ELASTICITY_FILE, "r") as f:
        elasticity_template = f.read()
    print(f"✓ Loaded elasticity template")

    with open(BASE_FRACTURE_FILE, "r") as f:
        fracture_template = f.read()
    print(f"✓ Loaded fracture template")

    print("\n" + "-" * 80)
    print("Generating parametric cases...")
    print("-" * 80)

    # Generate all combinations
    total_cases = 0
    for cf_key, cf_value in CONFINEMENT_PRESSURES.items():
        for domain_key, domain_config in DOMAIN_SIZES.items():
            case_name = FOLDER_NAME_PATTERN.format(cf_key=cf_key, domain_key=domain_key)

            create_case_folder(
                case_name,
                cf_value,
                domain_config,
                elasticity_template,
                fracture_template,
            )
            total_cases += 1

    print("\n" + "=" * 80)
    print(f"✓ Successfully generated {total_cases} parametric study cases!")
    print("=" * 80)

    # Print summary
    print("\nSummary:")
    print(
        f"  Confinement pressures: {list(CONFINEMENT_PRESSURES.keys())} = {list(CONFINEMENT_PRESSURES.values())}"
    )
    print(f"  Domain sizes: {list(DOMAIN_SIZES.keys())}")
    print(f"  Total cases: {total_cases}")

    print("\nGenerated folders:")
    for cf_key in CONFINEMENT_PRESSURES.keys():
        for domain_key in DOMAIN_SIZES.keys():
            case_name = FOLDER_NAME_PATTERN.format(cf_key=cf_key, domain_key=domain_key)
            print(f"  - {case_name}/")

    print("\n" + "=" * 80)
    print("Notes:")
    print(
        f"  - To modify parameters, edit the CONFIGURATION section at the top of this script"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
