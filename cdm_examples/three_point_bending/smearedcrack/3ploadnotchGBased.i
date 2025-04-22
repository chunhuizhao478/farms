[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  'meshnotch.msh'
  []
  [./elastic_region_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    top_right = '0.15 0.1 0'
    block_id = 2
  []
  [./elastic_region_2]
    type = SubdomainBoundingBoxGenerator
    input = elastic_region_1
    bottom_left = '0.30 0 0'
    top_right = '0.45 0.1 0'
    block_id = 2
  []
  [./extranodeset_1]
    type = ExtraNodesetGenerator
    coord = '   0 0 0'
    new_boundary = support1
    input = elastic_region_2
  []
  [./extranodeset_2]
    type = ExtraNodesetGenerator
    coord =  '0.45 0 0'
    new_boundary = support2
    input = extranodeset_1
  []
  [./extranodeset_3]
    type = ExtraNodesetGenerator
    coord = '0.22 0.1 0'
    new_boundary = load
    input = extranodeset_2
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
    [nonlocal_eqstrain]
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
  
  [Functions]
    [func_loading]
      type = ParsedFunction
      expression = '-1e-4 * t'
    []
  []

  [GlobalParams]
    displacements = 'disp_x disp_y'
  []

  [Physics/SolidMechanics/QuasiStatic]
    [./all]
      strain = FINITE
      add_variables = true
      generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_zx'
    [../]
  []

  [Kernels]
    # gradient based nonlocal averaging
    [react_nonlocal]
      type = Reaction
      variable = nonlocal_eqstrain
      rate = 1.0
    []
    [diffusion_nonlocal]
        type = CoefDiffusion
        variable = nonlocal_eqstrain
        coef = 6.25e-6 #5e-3
    []
    [reaction_local]
        type = ElkLocalEqstrainForce
        variable = nonlocal_eqstrain
    []    
  []
  
  [BCs]
    [fix_support_y]
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
    [apply_load]
      type = FunctionDirichletBC
      variable = disp_y
      boundary = load
      function = func_loading
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 20e9
      poissons_ratio = 0.2
    [../]
    [./elastic_stress]
        type = FarmsComputeSmearedCrackingStressGrads
        paramA = 0.99
        paramB = 100
        eqstrain_nonlocal = nonlocal_eqstrain
        cracking_stress = strength
        output_properties = 'crack_damage stress'
        outputs = exodus
    [../]
  []

  [Preconditioning]
    [smp]
      type = SMP
      full = true
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = PJFNK #in smeared cracking w/o full jacobian, this is much more efficient than NEWTON
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    line_search = 'none'
    # num_steps = 1
    l_max_its = 200
    nl_max_its = 100
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    dt = 0.01
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []