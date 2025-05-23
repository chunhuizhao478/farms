# Verification of Benchmark Problem TPV205-3D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #

# [Version 3] #
# This file serves to: #
# Generalize the computation of sticking traction using consistent displacement jump and nodal reaction forces #

[Mesh] 
    [./msh]
      type = GeneratedMeshGenerator
      dim = 2
      nx = 150
      ny = 150
      xmin = -15000
      xmax = 15000
      ymin = -15000
      ymax = 15000
      #elem_type = QUAD9
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
    displacements = 'disp_x disp_y' 
    fluid_vel = 'fluid_vel_x fluid_vel_y'
    porepressure = 'p'
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
    #area = 100
    mu_d = 0.525
  []
  
  [Variables]
    [./disp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./disp_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./fluid_vel_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./fluid_vel_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./p]
        order = FIRST
        family = LAGRANGE
    [../]
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
    [./disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./disp_slipweakening_y]
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
    [./fluid_disp_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./fluid_disp_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./fluid_vel_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./fluid_vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./interface_pressure]
        order = FIRST
        family = LAGRANGE
    []
  []
  
  [Actions/PoroCohesiveZoneAction]
      [./czm_ik]
          boundary = 'Block0_Block1'
          strain = SMALL
          p = p   
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
    [Fluid_Vel_x]
      type = CompVar
      variable = fluid_vel_slipweakening_x
      coupled = fluid_vel_x
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Fluid_Vel_y]
      type = CompVar
      variable = fluid_vel_slipweakening_y
      coupled = fluid_vel_y
      execute_on = 'TIMESTEP_BEGIN'
    []
    [Fluid_disp_x]
      type = VariableTimeIntegrationAux
      variable_to_integrate = fluid_vel_x
      variable = fluid_disp_slipweakening_x
      order = 2
      execute_on = 'TIMESTEP_BEGIN'
    [../]
    [Fluid_disp_y]
      type = VariableTimeIntegrationAux
      variable_to_integrate = fluid_vel_y
      variable = fluid_disp_slipweakening_y
      order = 2
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
    [pressure_global_plus]
        type = ProjectionAux
        variable = interface_pressure
        v = p
        boundary = 'Block0_Block1'
    []
  []
  
[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_y
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_x]
        type = CoupledFluidInertialForce
        variable = disp_x
        fluid_vel = fluid_vel_x
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_y]
        type = CoupledFluidInertialForce
        variable = disp_y
        fluid_vel = fluid_vel_y
        use_displaced_mesh = false
    [../]
    [./darcyflow_x]
        type = DynamicDarcyFlow2
        variable = fluid_vel_x
        skeleton_acceleration = disp_x
    [../]
    [./darcyflow_y]
        type = DynamicDarcyFlow2
        variable = fluid_vel_y
        skeleton_acceleration = disp_y
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_x
        porepressure = p
        component = 0
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_x
        porepressure = p
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_y
        porepressure = p
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p
        displacements = 'disp_x disp_y'
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p
        u = fluid_vel_x
        v = fluid_vel_y
        pressure = p
    [../]
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
    [./Reaction_x_porocoupling]
        type = PoroPropDampingCoupling
        variable = 'disp_x'
        component = '0'
    []
    [./Reaction_y_porocoupling]
        type = PoroPropDampingCoupling
        variable = 'disp_y'
        component = '1'
    []
[]

  
  [Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 35.7e9
        shear_modulus = 32e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [Strain]
        type = ComputeSmallStrain
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2320
    []
    [./rhof]
        type = GenericConstantMaterial
        prop_names = rhof
        prop_values = 1000
    [../]
    [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 2.236067977
    [../]
    [./porosity]
        type = GenericConstantMaterial
        prop_names = porosity
        prop_values = 0.2
    [../]
    [./hydconductivity]
        type = GenericConstantMaterial
        prop_names = hydconductivity
        prop_values = 1.128052989e-9
    [../]
    [./hydconductivity_layer]
        type = GenericConstantMaterial
        prop_names = hydconductivity_layer
        prop_values = 1.128052989e-9
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.567
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 10.084081087e+9
    [../]
    [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
    [../]
    [./czm_mat]
        type = PermeableSlipWeakening2dv3
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        fluid_disp_slipweakening_x  = fluid_disp_slipweakening_x
        fluid_disp_slipweakening_y  = fluid_disp_slipweakening_y
        fluid_vel_slipweakening_x  = fluid_vel_slipweakening_x
        fluid_vel_slipweakening_y  = fluid_vel_slipweakening_y
        pressure_interface = interface_pressure
        nodal_area = nodal_area
        boundary = 'Block0_Block1'
        p = p
    [../]
  []
  
  [BCs]
    [./fault_p]
       type = FunctionDirichletBC
       variable = fluid_vel_y
       boundary = 'Block0_Block1'
       function = 0.0
    [../]
    [./fault_n]
      type = FunctionDirichletBC
      variable = fluid_vel_y
      boundary = 'Block0_Block1'
      function = 0.0
    [../]
      [./flux_top]
      type = FunctionDirichletBC
      variable = fluid_vel_y
      boundary = top
      function = 0.0
    [../]
    [./flux_bot]
      type = FunctionDirichletBC
      variable = fluid_vel_y
      boundary = bottom
      function = 0.0
    [../]
    [./flux_left]
      type = FunctionDirichletBC
      variable = fluid_vel_x
      boundary = left
      function = 0.0
    [../]
    [./flux_right]
      type = FunctionDirichletBC
      variable = fluid_vel_x
      boundary = right
      function = 0.0
    [../]
  []

  [UserObjects]
    [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = 'Block0_Block1'
      execute_on = 'initial TIMESTEP_BEGIN'
    [../]
  []
  
  [Executioner]
    type = Transient
    dt = 0.0025
    end_time = 4
    [TimeIntegrator]
      type = CentralDifference
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
      input_files = 'poro_tpv2052D_v3_sub.i'
      execute_on = 'TIMESTEP_BEGIN'
    [../]
  []

  [Transfers]
    [pull_resid]
      type = MultiAppCopyTransfer
      from_multi_app = sub_app
      source_variable = 'resid_sub_x resid_sub_y'
      variable = 'resid_x resid_y'
      execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'disp_x disp_y'
      variable = 'disp_sub_x disp_sub_y'
      execute_on = 'TIMESTEP_BEGIN'
    []
    [push_fluid_vel]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'fluid_vel_x fluid_vel_y'
      variable = 'fluid_vel_sub_x fluid_vel_sub_y'
      execute_on = 'TIMESTEP_BEGIN'
    []
    [push_pressure]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'p'
      variable = 'p_sub'
      execute_on = 'TIMESTEP_BEGIN'
    []
  []