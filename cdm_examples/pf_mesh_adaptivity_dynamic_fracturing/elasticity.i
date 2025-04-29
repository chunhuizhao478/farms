E = 3.2e10
nu = 0.2
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'

Gc = 3
l = 2e-4 # N * h, N: number of elements, h: element size

[Adaptivity]
  max_h_level = 3
  marker = 'combo'
  cycles_per_step = 2
  [Markers]
      [./combo]
          type = ComboMarker
          markers = 'damage_marker strain_energy_marker'
      [../]
      [damage_marker]
        type = ValueThresholdMarker
        variable = d
        refine = 0.1
      []
      [strain_energy_marker]
        type = ValueThresholdMarker
        variable = psie_active
        refine = '${fparse 0.4*0.5*Gc/l}'
      []      
  []
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc=${Gc};l=${l}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = psie_active
    source_variable = psie_active
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file =  'mesh.msh'
  []
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [fy]
  []
  [d]
    family = LAGRANGE
    order = FIRST
  []
  #err measurement of active strain energy
  [eng_err]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
[]

[AuxKernels]
  [./error_measure]
    type = ErrorPsiMeasure
    variable = eng_err
  [../]
    [accel_x]
      type = NewmarkAccelAux
      variable = accel_x
      displacement = disp_x
      velocity = vel_x
      beta = 0.25
      execute_on = 'TIMESTEP_END'
  []
  [vel_x]
      type = NewmarkVelAux
      variable = vel_x
      acceleration = accel_x
      gamma = 0.5
      execute_on = 'TIMESTEP_END'
  []
  [accel_y]
      type = NewmarkAccelAux
      variable = accel_y
      displacement = disp_y
      velocity = vel_y
      beta = 0.25
      execute_on = 'TIMESTEP_END'
  []
  [vel_y]
      type = NewmarkVelAux
      variable = vel_y
      acceleration = accel_y
      gamma = 0.5
      execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    save_in = fy
  []
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    acceleration = accel_x
    velocity = vel_x
    beta = 0.25
    gamma = 0.5
    eta = 0
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    acceleration = accel_y
    velocity = vel_y
    beta = 0.25
    gamma = 0.5
    eta = 0
  []
[]

[Functions]
  [func_loading]
    type = ParsedFunction
    expression = '1e-7 * t'
  []
[]

[BCs]
  [load_top]
    type = NeumannBC
    variable = disp_y
    boundary = 3
    value = 1e6
  []
  [load_bottom]
    type = NeumannBC
    variable = disp_y
    boundary = 4
    value = -1e6
  [] 
[]

[Materials]
  [bulk]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [degradation]
    type = PowerDegradationFunction
    property_name = g
    expression = (1-d)^p*(1-eta)+eta
    phase_field = d
    parameter_names = 'p eta '
    parameter_values = '2 1e-6'
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = None
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [density]
    type = ADGenericConstantMaterial
    prop_names = 'density'
    prop_values = '2450'
  []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10

  dt = 1e-6
  end_time = 1000

  fixed_point_max_its = 5
  accept_on_max_fixed_point_iteration = true
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10

  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
[]

[Postprocessors]
  [Fy]
    type = NodalSum
    variable = fy
    boundary = 1
  []
[]

[Outputs]
  exodus = true
  print_linear_residuals = false
  csv = true
[]
