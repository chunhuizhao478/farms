# Validation test for ConfigurationalForceUserObject
# Simple penny-shaped crack under Mode I (tensile) loading
#
# Geometry: 3D cube with a small central crack perpendicular to z-axis
# Loading: Uniform tension in z-direction
# Expected: Configurational force primarily in z-direction (Mode I opening)
#           Force magnitude should be positive (crack wants to propagate)
#           Fx and Fy should be near zero due to symmetry

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 10
    ny = 10
    nz = 10
    xmin = -1.0
    xmax = 1.0
    ymin = -1.0
    ymax = 1.0
    zmin = -1.0
    zmax = 1.0
  []

  # Create a crack plane in the center (z=0)
  [crack_plane]
    type = ParsedGenerateSideset
    input = gen
    combinatorial_geometry = 'abs(z) < 0.05 & (x*x + y*y) < 0.25'
    new_sideset_name = crack_surface
    normal = '0 0 1'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxVariables]
  # Configurational force components
  [force_x]
    order = FIRST
    family = LAGRANGE
  []
  [force_y]
    order = FIRST
    family = LAGRANGE
  []
  [force_z]
    order = FIRST
    family = LAGRANGE
  []
  [force_mag]
    order = FIRST
    family = LAGRANGE
  []

  # Stress and energy for debugging
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [energy_density]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [TensorMechanics]
    # Standard tensor mechanics for small strain
  []
[]

[AuxKernels]
  [force_x_aux]
    type = ConfigurationalForceAux
    variable = force_x
    configurational_force_uo = config_force
    component = x
    execute_on = 'TIMESTEP_END'
  []

  [force_y_aux]
    type = ConfigurationalForceAux
    variable = force_y
    configurational_force_uo = config_force
    component = y
    execute_on = 'TIMESTEP_END'
  []

  [force_z_aux]
    type = ConfigurationalForceAux
    variable = force_z
    configurational_force_uo = config_force
    component = z
    execute_on = 'TIMESTEP_END'
  []

  [force_mag_aux]
    type = ConfigurationalForceAux
    variable = force_mag
    configurational_force_uo = config_force
    component = magnitude
    execute_on = 'TIMESTEP_END'
  []

  [stress_zz_aux]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []

  [energy_aux]
    type = MaterialRealAux
    property = elastic_energy
    variable = energy_density
  []
[]

[BCs]
  # Fix bottom in z
  [fix_z_bottom]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  []

  # Apply tensile displacement on top
  [pull_z_top]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'top'
    function = '0.01*t'  # Linear ramp: 0.01 m displacement
  []

  # Symmetry conditions to prevent rigid body motion
  [fix_x_corner]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0.0
  []

  [fix_y_corner]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
[]

[Materials]
  # Simple linear elastic material
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9  # 1 GPa (typical soft material)
    poissons_ratio = 0.3
  []

  [strain]
    type = ComputeSmallStrain
  []

  [stress]
    type = ComputeLinearElasticStress
  []

  # Compute elastic energy density for configurational force
  [elastic_energy]
    type = ElasticEnergyMaterial
    outputs = exodus
  []
[]

[UserObjects]
  [config_force]
    type = ConfigurationalForceUserObject

    # Material properties
    energy_density = elastic_energy
    stress = stress

    # Use simple formulation (no displacement gradient)
    use_displacement_gradient = false

    # Optional: only compute near crack
    # crack_front_boundaries = 'crack_surface'

    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  # Maximum configurational force magnitude
  [max_force]
    type = NodalExtremeValue
    variable = force_mag
    value_type = max
  []

  # Average force near crack (z≈0)
  [avg_force_crack]
    type = ElementAverageValue
    variable = force_mag
  []

  # Force components at center node (should be origin)
  [force_x_center]
    type = PointValue
    variable = force_x
    point = '0 0 0'
  []

  [force_y_center]
    type = PointValue
    variable = force_y
    point = '0 0 0'
  []

  [force_z_center]
    type = PointValue
    variable = force_z
    point = '0 0 0'
  []

  # Check energy and stress
  [max_energy]
    type = ElementExtremeValue
    variable = energy_density
  []

  [max_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'

  # Small timesteps for quasi-static loading
  dt = 0.1
  end_time = 1.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  l_max_its = 100
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false

  [console]
    type = Console
    max_rows = 10
  []
[]
