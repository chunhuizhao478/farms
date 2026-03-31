# Simple validation test for ConfigurationalForceUserObject
# 3D bar under uniaxial tension (no actual crack for now)
# This tests the basic functionality of computing Eshelby stress and nodal forces
#
# Expected results:
# - Forces should be non-zero in tension direction
# - Energy density should increase with loading
# - Eshelby stress Σ = ΨI - σ should be computed correctly

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 5
    ny = 5
    nz = 10
    xmin = 0
    xmax = 0.5
    ymin = 0
    ymax = 0.5
    zmin = 0
    zmax = 1.0
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
  []
  [force_y]
  []
  [force_z]
  []
  [force_mag]
  []

  # For debugging
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_energy_aux]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [TensorMechanics]
    # Uses GlobalParams for displacements
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

  [stress_xx_aux]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
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
    property = strain_energy_density
    variable = strain_energy_aux
  []
[]

[BCs]
  # Fix bottom face
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom'
    value = 0.0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  []

  # Apply displacement on top
  [top_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'top'
    function = '0.005*t'  # 0.5% strain at t=1
  []
[]

[Materials]
  # Elasticity tensor
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 10e9  # 10 GPa
    poissons_ratio = 0.3
  []

  # Strain
  [strain]
    type = ComputeSmallStrain
  []

  # Stress
  [stress]
    type = ComputeLinearElasticStress
  []

  # Strain energy density Ψ = 0.5 * σ : ε
  [strain_energy]
    type = StrainEnergyDensity
  []
[]

[UserObjects]
  [config_force]
    type = ConfigurationalForceUserObject

    # Material properties
    energy_density = strain_energy_density
    stress = stress

    # Simple formulation
    use_displacement_gradient = false

    execute_on = 'TIMESTEP_END'
  []
[]

[Postprocessors]
  # Track configurational forces
  [max_force_mag]
    type = NodalExtremeValue
    variable = force_mag
    value_type = max
  []

  [avg_force_z]
    type = ElementAverageValue
    variable = force_z
  []

  # Track mechanical quantities
  [max_stress_zz]
    type = ElementExtremeValue
    variable = stress_zz
    value_type = max
  []

  [max_energy]
    type = ElementExtremeValue
    variable = strain_energy_aux
    value_type = max
  []

  # Reaction force on top
  [react_z]
    type = NodalSum
    variable = disp_z
    boundary = 'top'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dt = 0.2
  end_time = 1.0

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]
