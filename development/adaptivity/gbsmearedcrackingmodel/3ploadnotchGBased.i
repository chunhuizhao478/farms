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
        scaling = 1e-5
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
        scaling = 1e-5
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
    [elastic_energy_aux]
      order = CONSTANT
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
    [get_elastic_energy]
      type = ElasticEnergyAux
      variable = elastic_energy_aux
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
        coef = 1e-9 #1e-4
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
    max_h_level = 2
    marker = 'combo'
    [Indicators]
        [error]
          type = GradientJumpIndicator
            variable = crack_damage_aux
        []
    []
    [Markers]
        [./combo]
            type = ComboMarker
            markers = 'error_marker elastic_energy_marker'
        [../]
        [./error_marker]
            type = ErrorFractionMarker
            indicator = error
            refine = 0.9
        [../]
        # peak energy: 144
        [./elastic_energy_marker]
            type = ValueThresholdMarker
            variable = elastic_energy_aux
            refine = 100
        []         
    []
  []
  
  [Executioner]
    type = Transient
    solve_type = PJFNK #in smeared cracking w/o full jacobian, this is much more efficient than NEWTON
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol' 
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9'
    automatic_scaling = true
    line_search = 'bt'
    # num_steps = 1
    l_max_its = 50
    nl_max_its = 100
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    dt = 0.0001
    # verbose = true
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []