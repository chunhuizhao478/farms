# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #
# [Note]: This serves as a test file, to run the full problem, please extend the domain size by modifying nx, ny, xmin, xmax, ymin, ymax

  [Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -8000
      xmax = 8000
      ymin = -2000
      ymax = 2000
      zmin = -10000
      zmax = 0
      nx = 80
      ny = 20
      nz = 100
      subdomain_ids = 1
    []
    [rename_boundary]
      type = RenameBoundaryGenerator
      input = msh
      old_boundary = 'top bottom left right front back'
      new_boundary = 'back front left right top bottom'
    []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = rename_boundary
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y > 0 & z >= -15000'
      block_id = 2
    []
    [./new_block_2]
      type = ParsedSubdomainMeshGenerator
      input = new_block_1
      combinatorial_geometry = 'x >= -15000 & x <= 15000 & y < 0 & z >= -15000'
      block_id = 3
    []       
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_2
        split_interface = true
        block_pairs = '2 3'
    []       
  []

  [GlobalParams]
    #primary variables
    displacements = 'disp_x disp_y disp_z'
    #damping ratio
    q = 0.4
    #characteristic length (m)
    Dc = 0.4
    #initial normal stress (Pa)
    T2_o = 120e6
    #initial shear stress along dip direction (Pa)
    T3_o = 0.0
    #dynamic friction coefficient
    mu_d = 0.525
    #element edge length (m)
    len = 200
  []

  [AuxVariables]
    [./resid_x]
      order = FIRST
      family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    []
    [./resid_z]
      order = FIRST
      family = LAGRANGE
    []
    [./resid_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_z]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_z]
      order = FIRST
      family = LAGRANGE
    []
    [./vel_slipweakening_x]
      order = FIRST
      family = LAGRANGE
    []
    [./vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_z]
      order = FIRST
      family = LAGRANGE
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_shear_stress]
        order = CONSTANT
        family = MONOMIAL
    []
    [./jump_x_rate]
        order = CONSTANT
        family = MONOMIAL
    []
    [./jump_y_rate]
      order = CONSTANT
      family = MONOMIAL
    []
    [./jump_z_rate]
      order = CONSTANT
      family = MONOMIAL
    []
  []

  [Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
      boundary = 'Block2_Block3'
      strain = SMALL
      generate_output='traction_x traction_y traction_z jump_x jump_y jump_z'
    [../]
  []

  [Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          add_variables = true
          generate_output = 'stress_xx stress_yy stress_xy'
          extra_vector_tags = 'restore_tag'
        []
      []
    []
  []

  [Problem]
    extra_tag_vectors = 'restore_tag'
  []

  [Functions]
    [func_static_friction_coeff_mus]
      type = PiecewiseConstant
      axis=x
      x = '-1000e3 -15e3 15e3'
      y = '10000 0.677 10000.0'
      direction = left
    []
    [func_initial_strike_shear_stress]
      type = InitialShearStressTPV2053d
    []
  []

  [AuxKernels]
    [Displacment_x]
      type = ProjectionAux
      variable = disp_slipweakening_x
      v = disp_x
      execute_on = 'TIMESTEP_END'
    []
    [Displacement_y]
      type = ProjectionAux
      variable = disp_slipweakening_y
      v = disp_y
      execute_on = 'TIMESTEP_END'
    []
    [Displacement_z]
      type = ProjectionAux
      variable = disp_slipweakening_z
      v = disp_z
      execute_on = 'TIMESTEP_END'
    []
    [Vel_x]
      type = CompVarRate
      variable = vel_slipweakening_x
      coupled = disp_x
      execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
      type = CompVarRate
      variable = vel_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_slipweakening_z
      coupled = disp_z
      execute_on = 'TIMESTEP_END'
    []
    [jump_x_rate]
      type = FDCompVarRate
      variable = jump_x_rate
      coupled = jump_x
      execute_on = 'TIMESTEP_END'
    []
    [jump_y_rate]
      type = FDCompVarRate
      variable = jump_y_rate
      coupled = jump_y
      execute_on = 'TIMESTEP_END'
    []
    [jump_z_rate]
      type = FDCompVarRate
      variable = jump_z_rate
      coupled = jump_z
      execute_on = 'TIMESTEP_END'
    []
    [restore_x]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_x'
      variable = 'resid_x'
      execute_on = 'TIMESTEP_END'
    []
    [restore_y]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_y'
      variable = 'resid_y'
      execute_on = 'TIMESTEP_END'
    []
    [restore_z]
      type = TagVectorAux
      vector_tag = 'restore_tag'
      v = 'disp_z'
      variable = 'resid_z'
      execute_on = 'TIMESTEP_END'
    []
    [StaticFricCoeff]
      type = FunctionAux
      variable = mu_s
      function = func_static_friction_coeff_mus
      execute_on = 'TIMESTEP_BEGIN'
    []
    [StrikeShearStress]
      type = FunctionAux
      variable = ini_shear_stress
      function = func_initial_strike_shear_stress
      execute_on = 'TIMESTEP_BEGIN'
    []
  []

  [Kernels]
    # [./inertia_x]
    #   type = InertialForce
    #   use_displaced_mesh = false
    #   variable = disp_x
    # []
    # [./inertia_y]
    #   type = InertialForce
    #   use_displaced_mesh = false
    #   variable = disp_y
    # []
    # [./inertia_z]
    #   type = InertialForce
    #   use_displaced_mesh = false
    #   variable = disp_z
    # []
    [massmatrix]
      type = MassMatrix
      density = 2670
      matrix_tags = 'system'
      variable = disp_x
    []
    [massmatrix_y]
      type = MassMatrix
      density = 2670
      matrix_tags = 'system'
      variable = disp_y
    []
    [massmatrix_z]
      type = MassMatrix
      density = 2670
      matrix_tags = 'system'
      variable = disp_z
    []  
    [./Reactionx]
      type = StiffPropDamping
      variable = 'disp_x'
      component = '0'
    []
    [./Reactiony]
      type = StiffPropDamping
      variable = 'disp_y'
      component = '1'
    []
    [./Reactionz]
      type = StiffPropDamping
      variable = 'disp_z'
      component = '2'
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
        prop_names = density
        prop_values = 2670
    []
    [./czm_mat]
        type = SlipWeakeningFrictionczm3d
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        disp_slipweakening_z     = disp_slipweakening_z
        vel_slipweakening_x      = vel_slipweakening_x
        vel_slipweakening_y      = vel_slipweakening_y
        vel_slipweakening_z      = vel_slipweakening_z
        reaction_slipweakening_x = resid_x
        reaction_slipweakening_y = resid_y
        reaction_slipweakening_z = resid_z
        mu_s = mu_s
        ini_shear_sts = ini_shear_stress
        boundary = 'Block2_Block3'
    [../]
  []

  [UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
  []

  [Executioner]
    type = Transient
    dt = 0.0025
    end_time = 12
    # num_steps = 1
    # [TimeIntegrator]
    #   type = CentralDifference
    #   solve_type = lumped
    #   use_constant_mass = true
    # []
    [TimeIntegrator]
      type = DirectCentralDifference
      mass_matrix_tag = 'system'
      use_constant_mass = true
    []
  []

  [Outputs]
    exodus = true
    time_step_interval = 40
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z traction_x traction_y traction_z jump_x jump_y jump_z jump_x_rate jump_y_rate jump_z_rate mu_s ini_shear_stress'
  []

  [BCs]
    ##non-reflecting bc
    #top surface is free surface
  #   [./dashpot_top_x]
  #     type = NonReflectDashpotBC3d
  #     component = 0
  #     variable = disp_x
  #     disp_x = disp_x
  #     disp_y = disp_y
  #     disp_z = disp_z
  #     p_wave_speed = 6000
  #     shear_wave_speed = 3464
  #     boundary = top
  # []
  # [./dashpot_top_y]
  #     type = NonReflectDashpotBC3d
  #     component = 1
  #     variable = disp_y
  #     disp_x = disp_x
  #     disp_y = disp_y
  #     disp_z = disp_z
  #     p_wave_speed = 6000
  #     shear_wave_speed = 3464
  #     boundary = top
  # []
  # [./dashpot_top_z]
  #     type = NonReflectDashpotBC3d
  #     component = 2
  #     variable = disp_z
  #     disp_x = disp_x
  #     disp_y = disp_y
  #     disp_z = disp_z
  #     p_wave_speed = 6000
  #     shear_wave_speed = 3464
  #     boundary = top
  # []
  [./dashpot_bottom_x]
      type = NonReflectDashpotBC3d
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_bottom_y]
      type = NonReflectDashpotBC3d
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_bottom_z]
      type = NonReflectDashpotBC3d
      component = 2
      variable = disp_z
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_left_x]
      type = NonReflectDashpotBC3d
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_left_y]
      type = NonReflectDashpotBC3d
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_left_z]
      type = NonReflectDashpotBC3d
      component = 2
      variable = disp_z
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_right_x]
      type = NonReflectDashpotBC3d
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
  [./dashpot_right_y]
      type = NonReflectDashpotBC3d
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
  [./dashpot_right_z]
      type = NonReflectDashpotBC3d
      component = 2
      variable = disp_z
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
  [./dashpot_front_x]
    type = NonReflectDashpotBC3d
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    p_wave_speed = 6000
    shear_wave_speed = 3464
    boundary = front
  []
  [./dashpot_front_y]
      type = NonReflectDashpotBC3d
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = front
  []
  [./dashpot_front_z]
      type = NonReflectDashpotBC3d
      component = 2
      variable = disp_z
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = front
  []
  [./dashpot_back_x]
    type = NonReflectDashpotBC3d
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    p_wave_speed = 6000
    shear_wave_speed = 3464
    boundary = back
  []
  [./dashpot_back_y]
      type = NonReflectDashpotBC3d
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = back
  []
  [./dashpot_back_z]
      type = NonReflectDashpotBC3d
      component = 2
      variable = disp_z
      disp_x = disp_x
      disp_y = disp_y
      disp_z = disp_z
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = back
  []
[]
