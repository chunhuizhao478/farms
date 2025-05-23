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

[Variables]
    [disp_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_vel_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [fluid_vel_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    [p_sub]
        order = FIRST
        family = LAGRANGE
    []
[]

[AuxVariables]
    [./resid_sub_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_sub_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y'
[]

[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_sub_x
        component = 0
        displacements = 'disp_sub_x disp_sub_y'
        use_displaced_mesh = false
        save_in = 'resid_sub_x'
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_sub_y
        component = 1
        displacements = 'disp_sub_x disp_sub_y'
        use_displaced_mesh = false
        save_in = 'resid_sub_y'
    [../]
    [./skeletoninertia_x]
        type = InertialForce
        variable = disp_sub_x
        use_displaced_mesh = false
    [../]
    [./skeletoninertia_y]
        type = InertialForce
        variable = disp_sub_y
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_x]
        type = CoupledFluidInertialForce
        variable = disp_sub_x
        fluid_vel = fluid_vel_sub_x
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_y]
        type = CoupledFluidInertialForce
        variable = disp_sub_y
        fluid_vel = fluid_vel_sub_y
        use_displaced_mesh = false
    [../]
    [./darcyflow_x]
        type = DynamicDarcyFlow2
        variable = fluid_vel_sub_x
        skeleton_acceleration = disp_sub_x
    [../]
    [./darcyflow_y]
        type = DynamicDarcyFlow2
        variable = fluid_vel_sub_y
        skeleton_acceleration = disp_sub_y
    [../]
    [./poromechskeletoncoupling_x]
        type = PoroMechanicsCoupling
        variable = disp_sub_x
        porepressure = p_sub
        component = 0
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_sub_y
        porepressure = p_sub
        component = 1
    [../]
    [./poromechfluidcoupling_x]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_sub_x
        porepressure = p_sub
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_vel_sub_y
        porepressure = p_sub
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p_sub
        displacements = 'disp_sub_x disp_sub_y'
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p_sub
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p_sub
        u = fluid_vel_sub_x
        v = fluid_vel_sub_y
        pressure = p_sub
    [../]
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
[]

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        #solve_type = lumped
    []
[]