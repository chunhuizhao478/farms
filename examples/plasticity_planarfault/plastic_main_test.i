#sample test geometry
[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmin = -5000
    xmax = 5000
    ymin = -5000
    ymax = 5000
    elem_type = TRI3
  []
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [./stsdev_x]
    type = StressDivergenceTensors
    variable = disp_x
    displacements = 'disp_x disp_y'
    component = 0
    use_displaced_mesh = false
  []
  [./stsdev_y]
    type = StressDivergenceTensors
    variable = disp_y
    displacements = 'disp_x disp_y'
    component = 1
    use_displaced_mesh = false
  []
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
[]

[Materials]
  [strain]
    type = ComputeIncrementalSmallStrain
    eigenstrain_names = eigenstrain_all
    displacements = 'disp_x disp_y'
  []
  [eigenstrain_all]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'func_initial_stress_xx           func_initial_strike_shear_stress 0  
                      func_initial_strike_shear_stress func_initial_stress_yy           0      
                      0 0 func_initial_stress_zz'
    eigenstrain_name = eigenstrain_all
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    lambda = 32.04e9
    shear_modulus = 32.04e9
    use_displaced_mesh = false
  []
  [./admissible]
    type = ComputeMultipleInelasticStress
    inelastic_models = cdp
    perform_finite_strain_rotations = false
  [../]
  [./cdp]
    type = CappedDruckerPragerStressUpdate
    DP_model = dp
    tensile_strength = ts
    compressive_strength = cs
    yield_function_tol = 1E-8
    tip_smoother = 4
    smoothing_tol = 1E-5
  [../]
[]

[Functions]
  #this function is used in medimum
  [func_initial_stress_xy_const]
    type = ConstantFunction
    value = 70e6
  []
  [./func_initial_stress_yy]
    type = ConstantFunction
    value = 120e6
  []
  #In problems with inelasticity, the sigma11 is important
  #This is different from pure tpv205 
  [./func_initial_stress_xx]
    type = ConstantFunction
    value = 135e6
  []
  [./func_initial_stress_zz]
    type = ConstantFunction
    value = 63.75e6
  []
[]

[UserObjects]
  [./ts]
    type = TensorMechanicsHardeningConstant
    value = 1000
  [../]
  [./cs]
    type = TensorMechanicsHardeningConstant
    value = 1000
  [../]
  [./mc_coh]
    type = TensorMechanicsHardeningConstant
    value = 4.3e6
  [../]
  [./mc_phi]
    type = TensorMechanicsHardeningConstant
    value = 30.5
    convert_to_radians = true
  [../]
  [./mc_psi]
    type = TensorMechanicsHardeningConstant
    value = 0
  [../]
  [./dp]
    type = TensorMechanicsPlasticDruckerPrager
    mc_cohesion = mc_coh
    mc_friction_angle = mc_phi
    mc_dilation_angle = mc_psi
    mc_interpolation_scheme = outer_tip
    internal_constraint_tolerance = 1 # irrelevant here
    yield_function_tolerance = 1      # irrelevant here
  [../]
[]

[Executioner]
  type = Transient
  dt = 1
  end_time = 4.0
  num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  interval = 20
[]