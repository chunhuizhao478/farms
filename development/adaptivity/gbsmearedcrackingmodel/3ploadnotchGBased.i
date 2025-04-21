[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  '../mesh/meshwohole.msh'
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
    [eqstrain_local_aux]
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

  [AuxKernels]
    [get_damage]
      type = MaterialRealAux
      variable = crack_damage_aux
      property = crack_damage
    []
    [get_eqstrain_local_aux]
      type = MaterialRealAux  
      variable = eqstrain_local_aux
      property = eqstrain_local
      execute_on = timestep_end
    []
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
        coef = 1e-10 #1e-4
    []
    [reaction_local]
        type = ElkLocalEqstrainForce
        variable = nonlocal_eqstrain
    []    
  []
  
  [BCs]
    [fix_right_x]
      type = DirichletBC
      variable = disp_x
      boundary = right
      value = 0
    []
    [fix_right_y]
      type = DirichletBC
      variable = disp_y
      boundary = right
      value = 0
    []
    [fix_load_x]
      type = FunctionDirichletBC
      variable = disp_x
      boundary = left
      function = func_loading
    []
    [fix_left_y]
      type = DirichletBC
      variable = disp_y
      boundary = left
      value = 0
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
        paramB = 500
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

  [Adaptivity]
    max_h_level = 3
    marker = 'combo'
    # initial_steps = 2
    # initial_marker = initial_box
    [Indicators]
        [error]
          type = GradientJumpIndicator
            variable = crack_damage_aux
        []
    []
    [Markers]
        [./combo]
            type = ComboMarker
            markers = 'error_marker'
        [../]
        [./error_marker]
            type = ErrorFractionMarker
            indicator = error
            refine = 0.9
        [../]
        # [./nonlocal_eqstrain_marker]
        #     type = ValueThresholdMarker
        #     variable = crack_damage_aux
        #     refine = 0.01
        # []         
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
    dt = 1e-5
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 10
  []