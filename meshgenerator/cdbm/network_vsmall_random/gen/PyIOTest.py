def modify_mesh_check(mesh_check_path, 
                      element_ids_path, 
                      subdomain_ids_path,
                      blockids_path,
                      output_path):
    
    # Read element IDs data
    with open(element_ids_path, 'r') as file:
        element_ids_data = file.read().strip()

    # Read subdomain IDs data
    with open(subdomain_ids_path, 'r') as file:
        subdomain_ids_data = file.read().strip()

    # Read block IDs data
    with open(blockids_path, 'r') as file:
        blockids_data = file.read().strip()   

    # Read the mesh check file and prepare to insert the element IDs data
    with open(mesh_check_path, 'r') as file:
        mesh_check_data = file.readlines()

    # Insert the element IDs into the appropriate place in the mesh check content
    for i, line in enumerate(mesh_check_data):
        if 'element_ids =' in line:
            # Use single quotes for the data
            mesh_check_data[i] = f"element_ids = \n        '{element_ids_data}'\n"
            break
    
    # Insert the element IDs into the appropriate place in the mesh check content
    for i, line in enumerate(mesh_check_data):
        if 'subdomain_ids =' in line:
            # Use single quotes for the data
            mesh_check_data[i] = f"subdomain_ids = \n        '{subdomain_ids_data}'\n"
            break

    # Insert the element IDs into the appropriate place in the mesh check content
    for i, line in enumerate(mesh_check_data):
        if 'surrounding_blocks =' in line:
            # Use single quotes for the data
            mesh_check_data[i] = f"surrounding_blocks = \n        '{blockids_data}'\n"
            break

    # Write the updated content back to the specified output file
    with open(output_path, 'w') as file:
        file.writelines(mesh_check_data)

# Usage
modify_mesh_check('../meshcheck/meshcheck_coarse.i', './elementid.txt', './subdomainid.txt', 'blockid.txt', '../meshcheck/meshcheck_coarse_out.i')