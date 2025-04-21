[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole.msh'
  []
  [./elastic_region_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0.003 0 0'
    top_right = '0.005 0.001 1'
    block_id = 1
  []
  [./elastic_region_2]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_1
    bottom_left = '0.023 0 0'
    top_right = '0.025 0.001 0'
    block_id = 1
  []
  [./elastic_region_3]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_2
    bottom_left = '0.013 0.007 0'
    top_right = '0.015 0.008 1'
    block_id = 1
  []
  [./extranodeset_0]
    type = ExtraNodesetGenerator
    coord = '0.004 0 0'
    new_boundary = support1
    input = elastic_region_3
    use_closest_node=true
  []
  [./extranodeset_1]
    type = ExtraNodesetGenerator
    coord = '0.024 0 0'
    new_boundary = support2
    input = extranodeset_0
    use_closest_node=true
  []
  [./extranodeset_2]
    type = ExtraNodesetGenerator
    coord = '0.014 0.008 0;
             0.0140986 0.008 0;
             0.0139014 0.008 0'
    new_boundary = load
    input = extranodeset_1
    use_closest_node=true
  []
  displacements = 'disp_x disp_y'
[]

  [Variables]
    [disp_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
    []
  []

  [AuxVariables]
    [./strength]
      order = CONSTANT
      family = MONOMIAL
      initial_condition = 2.4e6
    [../]
    [crack_damage_aux]
      order = FIRST
      family = MONOMIAL
    [] 
  []

  [AuxKernels]
    [get_damage]
      type = MaterialRealAux
      variable = crack_damage_aux
      property = crack_damage
      block = '0'
    []
  []
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4 * t'
    []
  []

  [GlobalParams]
    displacements = 'disp_x disp_y'
  []
  
  [BCs]
    [fix_support1_x]
      type = DirichletBC
      variable = disp_x
      boundary = support1
      value = 0
    []
    [fix_support1_y]
      type = DirichletBC
      variable = disp_y
      boundary = support1
      value = 0
    []
    [fix_support2_x]
      type = DirichletBC
      variable = disp_x
      boundary = support2
      value = 0
    []
    [fix_support2_y]
      type = DirichletBC
      variable = disp_y
      boundary = support2
      value = 0
    []
    [fix_load]
      type = FunctionDirichletBC
      variable = disp_y
      boundary = load
      function = func_loading
    []
    [fix_load_x]
      type = DirichletBC
      variable = disp_x
      boundary = load
      value = 0
    []
  []

  [Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = true
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = true
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 20e9
      poissons_ratio = 0.2
    [../]
    [./damage_stress]
        type = FarmsCrackBandModel
        cracking_stress = strength
        hb = 1.414e-4 #m hb = h sqrt(3) / 2
        Gf = 90 #J/m^2
        output_properties = 'crack_damage stress elastic_strain crack_initiation_strain'
        outputs = exodus
        block = '0'
    [../]
    [./elastic_stress]
        type = ComputeFiniteStrainElasticStress
        output_properties = 'stress'
        outputs = exodus
        block = '1'
    [../]      
    [strain]
      type = ComputeFiniteStrain
      displacements = 'disp_x disp_y'
    []
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
 []
  
  [Executioner]
    type = Transient
    solve_type = NEWTON
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    petsc_options_iname  = '-ksp_type -pc_type -pc_hypre_type'
    petsc_options_value  = 'gmres hypre boomeramg'
    # petsc_options_iname = '-ksp_type -pc_type -pc_type -pc_hypre_boomeramg_agg_nl' 
    # petsc_options_value = 'flexible multigrid hypre 2'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'    
    line_search = 'none'
    # num_steps = 1
    l_max_its = 200
    nl_max_its = 100
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    dt = 0.0001
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []