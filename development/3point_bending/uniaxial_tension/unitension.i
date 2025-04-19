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
  []

  [AuxVariables]
    [./strength]
      order = CONSTANT
      family = MONOMIAL
      initial_condition = 6.43e6
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

  [UserObjects]
    [eqstrain_averaging]
        type = ElkRadialAverage
        length_scale = 2e-5
        prop_name = eqstrain_local
        radius = 4e-5
        weights = BAZANT
        execute_on = LINEAR
    []
  []
  
  [BCs]
    [fix_support_x]
      type = DirichletBC
      variable = disp_x
      boundary = right
      value = 0
    []
    [fix_support_y]
      type = DirichletBC
      variable = disp_y
      boundary = right
      value = 0
    []
    [apply_load]
      type = FunctionDirichletBC
      variable = disp_x
      boundary = left
      function = func_loading
    []
  []
  
  [Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 40e6
      poissons_ratio = 0.25
    [../]
    [./elastic_stress]
        type = ComputeSmearedCrackingStressDebugIBased
        model = NONLOCAL
        damage_evolution_law_span = 2.0
        cracking_stress = strength
        softening_models = abrupt_softening
        cracked_elasticity_type = FULL
        max_cracks = 1
        output_properties = 'crack_damage stress elastic_strain crack_initiation_strain'
        outputs = exodus
    [../]
    # [./elastic_stress]
    #     type = ElkComputeSmearedCrackingStressModifiedMazars
    #     cracking_stress = strength
    #     paramA = 0.99
    #     paramB = 10000
    #     model = 'NONLOCAL'
    #     output_properties = 'crack_damage stress eqstrain_local'
    #     outputs = exodus
    # [../]
    # [./elastic_stress]
    #     type = ComputeSmearedCrackingStress
    #     cracking_stress = 6.43e6
    #     softening_models = abrupt_softening
    #     cracked_elasticity_type = FULL
    #     max_cracks = 1
    #     outputs = exodus
    # [../]
    [./abrupt_softening]
      type = AbruptSoftening
    [../]
    [./exponential_softening]
      type = ExponentialSoftening
    [../]
    [nonlocal_eqstrain]
      type = ElkNonlocalEqstrain
      average_UO = eqstrain_averaging
      output_properties = 'eqstrain_nonlocal'
      outputs = exodus
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
    solve_type = Newton
    petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    petsc_options_value = '101                asm      lu'
    line_search = 'none'
    # num_steps = 1
    l_max_its = 200
    nl_max_its = 10
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    l_tol = 1e-5
    start_time = 0.0
    end_time = 100
    # dt = 0.01
    [./TimeIntegrator]
      type = ImplicitEuler
    [../]
    [TimeStepper]
      type = IterationAdaptiveDT
      cutback_factor_at_failure = 0.5
      growth_factor = 1.5
      optimal_iterations = 5
      linear_iteration_ratio = 200
      dt = 0.01
    []
  []
  
  [Outputs]
    exodus = true
    time_step_interval = 1
  []