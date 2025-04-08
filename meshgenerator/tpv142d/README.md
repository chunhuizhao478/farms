Generate TPV14 Mesh
(1) Run tpv142d_100m.geo using Gmsh

(2) Run meshio_run_tria.py to extract the element data
- "raw_results" (raw data files contains all elements above the faults)
- res_ind_elem_crack_mf.txt 
- res_ind_elem_crack_bf.txt
- res_xcoord_crack_mf.txt
- res_ycoord_crack_mf.txt
- res_xcoord_crack_bf.txt
- res_ycoord_crack_bf.txt
(3) Implement the elementIDs, subdomainIDs in "meshcheck" file

(4) Run process_coords.py to extract element data from raw_results (for postprocess)
- "preprocess_results_mf"
- "preprocess_results_bf"
- results_mf.txt                  - coordinates
- results_ids_mf.txt              - elementid -> add in vectorpostprocessor ElementMaterialSampler
- result_subdomainid_mf.txt       - subdomainid (for check)
- results_bf.txt                  - coordinates
- results_ids_bf.txt              - elementid -> add in vectorpostprocessor ElementMaterialSampler
- result_subdomainid_bf.txt       - subdomainid (for check)
(5) Implement the elementIDs, subdomainIDs in "meshcheck_sorted" file