#!/usr/bin/env python3
"""
Parametric Study Case Generator for Porousflow Simulations

This script generates parametric study cases by creating folders and modified
input files (elasticity.i, fracture.i, and static_solve.i) based on different
combinations of:
1. Confinement pressure variations (with varying domain sizes)
2. Initial pore pressure variations (with fixed confinement at 5 MPa and varying domain sizes)

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
BASE_STATIC_SOLVE_FILE = BASE_DIR / "original_files" / "case_name" / "static_solve.i"

# ==============================================================================
# STUDY 1: Confinement Pressure Variations
# ==============================================================================

# Confinement pressures to study
CONFINEMENT_PRESSURES = {
    "cf1": 1e6,  # 1 MPa
    "cf5": 5e6,  # 5 MPa
    "cf10": 10e6,  # 10 MPa
}

# ==============================================================================
# STUDY 2: Initial Pore Pressure Variations (Fixed Confinement at 5 MPa)
# ==============================================================================

# Initial pore pressures to study (with FIXED confinement at 5 MPa)
INITIAL_PORE_PRESSURES = {
    "pp0d0965": 0.0965e6,  # 0.0965 MPa
    "pp2": 2e6,  # 2 MPa
    "pp4": 4e6,  # 4 MPa
}

# Fixed confinement pressure for pore pressure study
FIXED_CONFINEMENT_FOR_PP_STUDY = 5e6  # 5 MPa

# ==============================================================================
# SHARED PARAMETERS: Domain Sizes (for both studies)
# ==============================================================================

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

# Folder naming patterns
CONFINEMENT_FOLDER_PATTERN = "case_{cf_key}_{domain_key}"
PORE_PRESSURE_FOLDER_PATTERN = "case_{pp_key}_{domain_key}"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================


def modify_static_solve_file(
    content, confinement_pressure, initial_pore_pressure, mesh_file, coord
):
    """
    Modify the static_solve.i file content with specified parameters.

    Args:
        content (str): Original file content
        confinement_pressure (float): Confinement pressure value (e.g., 1e6)
        initial_pore_pressure (float): Initial pore pressure value (e.g., 0.0965e6)
        mesh_file (str): Path to mesh file
        coord (str): Coordinate string for fixed point

    Returns:
        str: Modified file content
    """
    lines = content.split("\n")
    modified_lines = []

    for i, line in enumerate(lines):
        # Modify confinement pressure (line 1)
        stripped = line.strip()
        if stripped.startswith("confinement_pressure") and "=" in stripped:
            modified_lines.append(f"confinement_pressure  = {confinement_pressure}")
        # Modify initial pore pressure (line 2)
        elif stripped.startswith("initial_pore_pressure") and "=" in stripped:
            modified_lines.append(f"initial_pore_pressure = {initial_pore_pressure}")
        # Modify mesh file path
        elif "file =" in line and "2dmeshfile" in line:
            # Check if this is in the Mesh section
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
        stripped = line.strip()
        if stripped.startswith("confinement_pressure") and "=" in stripped:
            modified_lines.append(f"confinement_pressure  = {confinement_pressure}")
        # Modify mesh file path
        elif "file =" in line and "2dmeshfile" in line:
            # Check if this is in the Mesh section
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
            modified_lines.append(" " * indent + f"file =  '{mesh_file}'")
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
    case_name,
    cf_value,
    domain_config,
    static_solve_template,
    elasticity_template,
    fracture_template,
    initial_pp_value=None,
):
    """
    Create a case folder with modified input files.

    Args:
        case_name (str): Name of the case folder
        cf_value (float): Confinement pressure value
        domain_config (dict): Domain configuration with 'mesh_file' and 'coord'
        static_solve_template (str): Template content for static_solve.i
        elasticity_template (str): Template content for elasticity.i
        fracture_template (str): Template content for fracture.i
        initial_pp_value (float, optional): Initial pore pressure value. If None, uses default from template
    """
    # Create case directory
    case_dir = BASE_DIR / case_name
    case_dir.mkdir(parents=True, exist_ok=True)

    # Determine initial pore pressure
    # If initial_pp_value is provided, use it; otherwise extract from static_solve_template
    if initial_pp_value is None:
        # Extract default initial pore pressure from template
        for line in static_solve_template.split("\n"):
            if line.strip().startswith("initial_pore_pressure ="):
                # Extract the value
                initial_pp_value = float(line.split("=")[1].strip())
                break

    # Modify files
    modified_static_solve = modify_static_solve_file(
        static_solve_template,
        cf_value,
        initial_pp_value,
        domain_config["mesh_file"],
        domain_config["coord"],
    )

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
    with open(case_dir / "static_solve.i", "w") as f:
        f.write(modified_static_solve)

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
    print("Parametric Study Case Generator for Porousflow Simulations")
    print("=" * 80)
    print(f"\nBase directory: {BASE_DIR}")
    print(f"\nTemplate files:")
    print(f"  - Static Solve: {BASE_STATIC_SOLVE_FILE}")
    print(f"  - Elasticity:   {BASE_ELASTICITY_FILE}")
    print(f"  - Fracture:     {BASE_FRACTURE_FILE}")

    # Check if base files exist
    if not BASE_STATIC_SOLVE_FILE.exists():
        print(f"\n✗ ERROR: Base static_solve file not found: {BASE_STATIC_SOLVE_FILE}")
        return

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

    with open(BASE_STATIC_SOLVE_FILE, "r") as f:
        static_solve_template = f.read()
    print(f"✓ Loaded static_solve template")

    with open(BASE_ELASTICITY_FILE, "r") as f:
        elasticity_template = f.read()
    print(f"✓ Loaded elasticity template")

    with open(BASE_FRACTURE_FILE, "r") as f:
        fracture_template = f.read()
    print(f"✓ Loaded fracture template")

    # ==============================================================================
    # STUDY 1: Generate confinement pressure variation cases
    # ==============================================================================
    print("\n" + "=" * 80)
    print("STUDY 1: Generating confinement pressure variation cases...")
    print("=" * 80)

    total_cf_cases = 0
    for cf_key, cf_value in CONFINEMENT_PRESSURES.items():
        for domain_key, domain_config in DOMAIN_SIZES.items():
            case_name = CONFINEMENT_FOLDER_PATTERN.format(
                cf_key=cf_key, domain_key=domain_key
            )

            create_case_folder(
                case_name,
                cf_value,
                domain_config,
                static_solve_template,
                elasticity_template,
                fracture_template,
                initial_pp_value=None,  # Use default from template
            )
            total_cf_cases += 1

    # ==============================================================================
    # STUDY 2: Generate initial pore pressure variation cases
    # ==============================================================================
    print("\n" + "=" * 80)
    print("STUDY 2: Generating initial pore pressure variation cases...")
    print(f"         (Fixed confinement pressure: {FIXED_CONFINEMENT_FOR_PP_STUDY} Pa)")
    print("=" * 80)

    total_pp_cases = 0
    for pp_key, pp_value in INITIAL_PORE_PRESSURES.items():
        for domain_key, domain_config in DOMAIN_SIZES.items():
            case_name = PORE_PRESSURE_FOLDER_PATTERN.format(
                pp_key=pp_key, domain_key=domain_key
            )

            create_case_folder(
                case_name,
                FIXED_CONFINEMENT_FOR_PP_STUDY,  # Fixed at 5 MPa
                domain_config,
                static_solve_template,
                elasticity_template,
                fracture_template,
                initial_pp_value=pp_value,  # Vary initial pore pressure
            )
            total_pp_cases += 1

    # ==============================================================================
    # Print Summary
    # ==============================================================================
    print("\n" + "=" * 80)
    print(
        f"✓ Successfully generated {total_cf_cases + total_pp_cases} parametric study cases!"
    )
    print("=" * 80)

    print("\n" + "-" * 80)
    print("STUDY 1: Confinement Pressure Variations")
    print("-" * 80)
    print(
        f"  Confinement pressures: {list(CONFINEMENT_PRESSURES.keys())} = {list(CONFINEMENT_PRESSURES.values())}"
    )
    print(f"  Domain sizes: {list(DOMAIN_SIZES.keys())}")
    print(f"  Total cases: {total_cf_cases}")
    print("\n  Generated folders:")
    for cf_key in CONFINEMENT_PRESSURES.keys():
        for domain_key in DOMAIN_SIZES.keys():
            case_name = CONFINEMENT_FOLDER_PATTERN.format(
                cf_key=cf_key, domain_key=domain_key
            )
            print(f"    - {case_name}/")

    print("\n" + "-" * 80)
    print("STUDY 2: Initial Pore Pressure Variations")
    print("-" * 80)
    print(f"  Fixed confinement: {FIXED_CONFINEMENT_FOR_PP_STUDY} Pa (5 MPa)")
    print(
        f"  Initial pore pressures: {list(INITIAL_PORE_PRESSURES.keys())} = {list(INITIAL_PORE_PRESSURES.values())}"
    )
    print(f"  Domain sizes: {list(DOMAIN_SIZES.keys())}")
    print(f"  Total cases: {total_pp_cases}")
    print("\n  Generated folders:")
    for pp_key in INITIAL_PORE_PRESSURES.keys():
        for domain_key in DOMAIN_SIZES.keys():
            case_name = PORE_PRESSURE_FOLDER_PATTERN.format(
                pp_key=pp_key, domain_key=domain_key
            )
            print(f"    - {case_name}/")

    print("\n" + "=" * 80)
    print("Notes:")
    print(f"  - Template files location: {BASE_DIR / 'original_files' / 'case_name'}/")
    print(
        f"  - To modify parameters, edit the CONFIGURATION section at the top of this script"
    )
    print("=" * 80)


if __name__ == "__main__":
    main()
