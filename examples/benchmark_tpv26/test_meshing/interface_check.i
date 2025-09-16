#parameters

##mesh parameters
###need to check gmsh file for changing the parameters here
bottom_nodes_coord =' -60000 -60000 -60000;
                      60000 -60000 -60000;
                      60000 60000  -60000;
                     -60000 60000  -60000'

##element size
elem_size = 400 #!!! element size near the fault, need to be consistent with the mesh file

##main fault parameters
xmin_fault = -16000 #xmin of fault
xmax_fault = 12000 #xmax of fault
zmin_fault = -15000 #zmin of fault
# zmax_fault = 0 #zmax of fault

##branch fault parameters
numerical_factor1 = 0.75
numerical_factor2 = 1.1
branch_fault_length = 12000 #length of branch fault
cutoff_xmin = '${fparse 1.732 * 0.5 * elem_size * numerical_factor1}' #xmin of branch fault
cutoff_ymax = '${fparse -1 * 0.5 * elem_size * numerical_factor1}' #ymin of branch fault
cutoff_xmax = '${fparse 1.732 * 0.5 * branch_fault_length * numerical_factor2}' #xmax of branch fault
cutoff_ymin = '${fparse -1 * 0.5 * branch_fault_length * numerical_factor2}' #ymax of branch fault
zmin_branch_fault = -15000 #zmin of branch fault
##-------------------------##
##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
##-------------------------##

##model parameters##
dt = 0.0025 #time step size

end_time = 12.0 #end time for simulation

# num_steps = 40 #end_time or num_steps only one of them is needed
exodus_time_step_interval = 1 #time step interval for output
csv_time_step_interval = 40 #time step interval for csv output
checkpoint_time_step_interval = 40 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##------------------------------------------------------------------------##

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../../meshgenerator/tpv24/TPV24_100m.msh'
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y > 0'
    block_id = 100
  []
  [./new_block_2]
    type = ParsedSubdomainMeshGenerator
    input = new_block_1
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y < 0'
    block_id = 200
  []             
  # [./new_block_3]
  #   type = ParsedSubdomainMeshGenerator
  #   input = new_block_2
  #   combinatorial_geometry = '0.5774 * x + y > 0 & x >= ${cutoff_xmin} & x <= ${cutoff_xmax} & y >= ${cutoff_ymin} & y <= ${cutoff_ymax} & z >= ${zmin_branch_fault}'
  #   block_id = 300    
  # []
  # [./new_block_4]
  #   type = ParsedSubdomainMeshGenerator
  #   input = new_block_3
  #   combinatorial_geometry = '0.5774 * x + y < 0 & x >= ${cutoff_xmin} & x <= ${cutoff_xmax} & y >= ${cutoff_ymin} & y <= ${cutoff_ymax} & z >= ${zmin_branch_fault}'
  #   block_id = 400  
  # []  
  [./new_block_3]
    type = ParsedSubdomainMeshGenerator
    input = new_block_2
    combinatorial_geometry = '(- 6155.000000 * x - 10660.772721 * y + 0.000000 < 0) & (4505.772721 * x + 4505.772721 * y - 20384449.173552 < 0) & (1649.227279 * x + 6155.000000 * y + 82461.363971 < 0) & z >= ${zmin_branch_fault}'
    block_id = 300    
  []
  [./new_block_4]
    type = ParsedSubdomainMeshGenerator
    input = new_block_3
    combinatorial_geometry = '(- 6155.000000 * x - 10660.772721 * y + 0.000000 > 0) & (1649.227279 * x + 6155.000000 * y + 20384449.173552 > 0) & (4505.772721 * x + 4505.772721 * y - 82461.363971 > 0) & z >= ${zmin_branch_fault}'
    block_id = 400  
  []
  [./split_1]
    type = BreakMeshByBlockGenerator
    input = new_block_4
    split_interface = true
    block_pairs = '100 200; 300 400'
  []
  [./sidesets]
    input = split_1
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0
                0 0 -1
                0 0 1'
    new_boundary = 'left right bottom top back front'
  [] 
  [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = ${bottom_nodes_coord}
      new_boundary = corner_ptr
      input = sidesets
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        generate_output = 'stress_xx stress_yy stress_xy'
      []
    []
  []
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_y
  []
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_z
  []
[]

[Materials]
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
      use_displaced_mesh = false
  []
  [stress]
      type = ComputeLinearElasticStress
  []
  [density]
      type = GenericConstantMaterial
      prop_names = 'density'
      prop_values = '${density}'
  []
[]

[Executioner]
  type = Transient
  dt = ${dt}
  end_time = ${end_time}
  # num_steps = ${num_steps}
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    execute_on = 'timestep_end'
    time_step_interval = ${exodus_time_step_interval}
  []
[]    