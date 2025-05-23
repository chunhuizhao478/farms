# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 200
        ny = 200
        xmin = -15000
        xmax = 15000
        ymin = -15000
        ymax = 15000
        elem_type = QUAD9
        
    []
    [./new_block]
        type = ParsedSubdomainMeshGenerator
        input = msh
        combinatorial_geometry = 'y>0'
        block_id = 1
    []
    [./split]
        type = BreakMeshByBlockGenerator
        input = new_block
        split_interface = true
        add_interface_on_two_sides = true
    []   
[]

[GlobalParams]
    displacements = 'disp_sub_x disp_sub_y' 
    Dc = 0.4
    elem_size = 150
    T2_o = 120e6
    mu_d = 0.525
[]

[Variables]
    [./disp_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_sub_vel_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_sub_vel_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p_sub]
        order = FIRST
        family = LAGRANGE
    [../]
[]


[AuxVariables]
    [./disp_czm_x]
        order = SECOND
        family = LAGRANGE
    []
    [./disp_czm_y]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_sub_x]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_sub_y]
        order = SECOND
        family = LAGRANGE
    []
    [./nodal_area]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damp_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damp_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p_sub_main]
        order = FIRST
        family = LAGRANGE
    [../]
    [./traction_sub_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./traction_sub_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./traction_sub_strike]
        order = SECOND
        family = MONOMIAL
    [../]
    [./traction_sub_normal]
        order = SECOND
        family = MONOMIAL
    [../]

[]

[AuxKernels]
    [mat_traction_strike]  # Changed name to avoid conflict
        type = MaterialRealAux
        property = traction_x
        variable = traction_sub_strike
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [mat_traction_normal]  # Changed name to avoid conflict
        type = MaterialRealAux
        property = traction_y
        variable = traction_sub_normal
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [mat_traction_strike_s]  # Changed name to avoid conflict
        type = MaterialRealAux
        property = traction_x
        variable = traction_sub_strike
        boundary = 'Block1_Block0'
        execute_on = 'TIMESTEP_END'
    []
    [mat_traction_normal_s]  # Changed name to avoid conflict
        type = MaterialRealAux
        property = traction_y
        variable = traction_sub_normal
        boundary = 'Block1_Block0'
        execute_on = 'TIMESTEP_END'
    []  
    [strike_tract]
        type = ProjectionAux
        v = traction_sub_strike
        variable = traction_sub_x
    []
    [normal_tract]
        type = ProjectionAux
        v = traction_sub_normal
        variable = traction_sub_y
    []
[]


[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y traction_global_x traction_global_y'
    [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik_2]
        boundary = 'Block1_Block0'
        strain = SMALL
        generate_output='traction_x traction_y traction_global_x traction_global_y'
    [../]
[]


[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_sub_x
        component = 0
        use_displaced_mesh = false   
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_sub_y
        component = 1
        use_displaced_mesh = false
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
        fluid_vel = fluid_sub_vel_x
        use_displaced_mesh = false
    [../]
    [./porefluidIFcoupling_y]
        type = CoupledFluidInertialForce
        variable = disp_sub_y
        fluid_vel = fluid_sub_vel_y
        use_displaced_mesh = false
    [../]
    [./darcyflow_x]
        type = DynamicDarcyFlow2
        variable = fluid_sub_vel_x
        skeleton_acceleration = disp_sub_x
    [../]
    [./darcyflow_y]
        type = DynamicDarcyFlow2
        variable = fluid_sub_vel_y
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
        variable = fluid_sub_vel_x
        porepressure = p_sub
        component = 0
    [../]
    [./poromechfluidcoupling_y]
        type = PoroMechanicsCoupling2
        variable = fluid_sub_vel_y
        porepressure = p_sub
        component = 1
    [../]
    [./massconservationskeleton]
        type = INSmassSolid
        variable = p_sub
    [../]
    [./massconservationpressure]
        type = FluidStorage
        variable = p_sub
    [../]
    [./massconservationfluid]
        type = INSmassFluid
        variable = p_sub
        u = fluid_sub_vel_x
        v = fluid_sub_vel_y
        pressure = p_sub
    [../]
[]


[Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 15.46e9
        shear_modulus = 13.86e9
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
        prop_values = 2.24
    [../]
    [./porosity]
        type = GenericConstantMaterial
        prop_names = porosity
        prop_values = 0.2
    [../]
    [./hydconductivity]
        type = GenericConstantMaterial
        prop_names = hydconductivity
        prop_values = 1.1280533319e-9
    [../]
    [./hydconductivity_layer]
        type = GenericConstantMaterial
        prop_names = hydconductivity_layer
        prop_values = 1.1280533319e-9
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.567
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 1.0084e10
    [../]
    [./constants]
        type = GenericConstantMaterial
        prop_names = 'rho mu'
        prop_values = '1  1'
    [../]
    [./czm_mat]
        type = PoroSlipWeakening2d_sub
        boundary = 'Block0_Block1'
        interface_disp_x = disp_czm_x
        interface_disp_y = disp_czm_y
        pressure_plus = p_sub_main
        pressure_minus = p_sub_main
        react_x = resid_sub_x
        react_y = resid_sub_y
        react_pressure_x = resid_pressure_sub_x
        react_pressure_y = resid_pressure_sub_y 
        react_damp_x = resid_damp_sub_x
        react_damp_y = resid_damp_sub_y
        nodal_area = nodal_area
        fluid_vel_x = fluid_disp_sub_x
        fluid_vel_y = fluid_disp_sub_y
        fluid_disp_x = fluid_disp_sub_x
        fluid_disp_y = fluid_disp_sub_y
        permeability_type = 'impermeable'
    [../]
    [./czm_mat_sec]
        type = PoroSlipWeakening2d_sub
        boundary = 'Block1_Block0'
        interface_disp_x = disp_czm_x
        interface_disp_y = disp_czm_y
        pressure_plus = p_sub_main
        pressure_minus = p_sub_main
        react_x = resid_sub_x
        react_y = resid_sub_y
        react_pressure_x = resid_pressure_sub_x
        react_pressure_y = resid_pressure_sub_y 
        react_damp_x = resid_damp_sub_x
        react_damp_y = resid_damp_sub_y
        nodal_area = nodal_area
        fluid_vel_x = fluid_disp_sub_x
        fluid_vel_y = fluid_disp_sub_y
        fluid_disp_x = fluid_disp_sub_x
        fluid_disp_y = fluid_disp_sub_y
        permeability_type = 'impermeable'
    [../]
    
[]


[UserObjects]
 #   compute element side volume (using CONTACT modulus)
     [element_side_volume]
         type = NodalArea
         variable = nodal_area
         boundary = 'Block0_Block1 Block1_Block0'
         execute_on = 'initial'
    []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
    type = Transient
    [TimeIntegrator]
        type = CentralDifference
        #solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]
