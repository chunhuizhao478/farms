# Verification of Benchmark Problem TPV205-3D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #

# [Version 3] #
# This file serves to: #
# Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces #

[Mesh]
    [./msh]
      type = GeneratedMeshGenerator
      dim = 3
      nx = 150
      ny = 150
      nz = 150
      xmin = -15000
      xmax = 15000
      ymin = -15000
      ymax = 15000
      zmin = -15000
      zmax = 15000
    []
    [./new_block]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'y<0'
      block_id = 1
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = new_block
      split_interface = true
    []
    [interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = split
        primary_block = 0
        paired_block = 1
        new_boundary = 'interface'
    []
    [secondary_interface]
        type = SideSetsBetweenSubdomainsGenerator
        input = interface
        primary_block = 1
        paired_block = 0
        new_boundary = 'secondary_interface'
    []
  []
    
  [GlobalParams]
    displacements = 'disp_x disp_y disp_z'
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
    T3_o = 0.0
    mu_d = 0.525
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
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
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
  []
  
  [Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
      boundary = 'Block0_Block1'
      strain = SMALL
    [../]
  []
    
    
  [Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          add_variables = true
        [../]
      [../]
    [../]
  []
  
  [AuxKernels]
    [Vel_x]
      type = CompVarRate
      variable = vel_slipweakening_x
      coupled = disp_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
      type = CompVarRate
      variable = vel_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_slipweakening_z
      coupled = disp_z
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacment_x]
      type = CompVar
      variable = disp_slipweakening_x
      coupled = disp_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_y]
      type = CompVar
      variable = disp_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_z]
      type = CompVar
      variable = disp_slipweakening_z
      coupled = disp_z
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_x]
      type = CompVar
      variable = resid_slipweakening_x
      coupled = resid_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
      type = CompVar
      variable = resid_slipweakening_y
      coupled = resid_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_z]
      type = CompVar
      variable = resid_slipweakening_z
      coupled = resid_z
      execute_on = 'TIMESTEP_BEGIN'
    []
  []
  
  [Kernels]
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
    [./inertia_z]
      type = InertialForce
      use_displaced_mesh = false
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
        type = SlipWeakeningFriction3dv3
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        disp_slipweakening_z     = disp_slipweakening_z
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        reaction_slipweakening_z = resid_slipweakening_z
        nodal_area = nodal_area
        boundary = 'Block0_Block1'
    [../]
  []
  
  [UserObjects]
    [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = 'interface secondary_interface'
      execute_on = 'initial TIMESTEP_BEGIN'
    [../]
  []
  
  [Executioner]
    type = Transient
    dt = 0.005
    end_time = 3.0
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
  []
  
  [Outputs]
    exodus = true
    interval = 10
  []

  [MultiApps]
    [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'tpv2053D_v3_sub.i'
      execute_on = 'TIMESTEP_BEGIN'
    [../]
  []

  [Transfers]
    [pull_resid]
      type = MultiAppCopyTransfer
      from_multi_app = sub_app
      source_variable = 'resid_sub_x resid_sub_y resid_sub_z'
      variable = 'resid_x resid_y resid_z'
      execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'disp_x disp_y disp_z'
      variable = 'disp_sub_x disp_sub_y disp_sub_z'
      execute_on = 'TIMESTEP_BEGIN'
    []
  []