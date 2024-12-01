# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  './Planar_fault_unstructured.msh'
    []
    [subdomain1]
        input = msh
        type = SubdomainBoundingBoxGenerator
        bottom_left = '-15000 -15000 0'
        top_right = '15000 15000 0'
        block_id = 0
    []
    [./new_block_1]
        type = ParsedSubdomainMeshGenerator
        input = subdomain1
        combinatorial_geometry = 'y > 0'
        block_id = 1
    []
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_1
        split_interface = true
        add_interface_on_two_sides = true
        block_pairs = '0 1'
    []     
[]

[GlobalParams]
    displacements = 'disp_x disp_y' 
    fluid_vel = 'fluid_vel_x fluid_vel_y'
    porepressure = 'p'
    q = 0.05
    Dc = 0.4
    elem_size = 12.5
    T2_o = 120e6
    mu_d = 0.525
[]

[Variables]
    [./disp_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_vel_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./fluid_vel_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./p]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxVariables]
    [./vel_x]
        order = SECOND
        family = LAGRANGE
    []
    [./accel_x]
    []
    [./vel_y]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_x]
        order = SECOND
        family = LAGRANGE
    []
    [./fluid_disp_y]
        order = SECOND
        family = LAGRANGE
    []
    [./accel_y]
    []
    [./nodal_area]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_primary_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_primary_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_primary_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_primary_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damping_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_damping_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_damping_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_damping_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./resid_pressure_y]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_pressure_x]
        order = SECOND
        family = LAGRANGE
    [../]
    [./jacob_pressure_y]
        order = SECOND
        family = LAGRANGE
    [../]
[]

[AuxKernels]
    [velocity_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
    []
    [velocity_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
    [] 
[]

[Actions/PoroCohesiveZoneAction]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y jump_vel_x jump_vel_y normal_traction tangent_traction normal_jump tangent_jump'
    [../]
[]


[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false   
        save_in = 'resid_primary_x' 
        save_in_diag = 'jacob_primary_x' 
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
        save_in = 'resid_primary_y' 
        save_in_diag = 'jacob_primary_y' 
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
        save_in = 'resid_pressure_x'
        save_in_diag = 'jacob_pressure_x' 
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
        save_in = 'resid_pressure_y'
        save_in_diag = 'jacob_pressure_y' 
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
        save_in = 'resid_damping_x'
        save_in_diag = 'jacob_damping_x' 
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = 'disp_y'
        component = '1'
        save_in = 'resid_damping_y'
        save_in_diag = 'jacob_damping_y' 
    []
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
    [./czm_stress_derivative]
        type = StressDerivative2
        boundary = 'Block0_Block1'
    [../]
    [./czm_mat]
        type = PoroSlipWeakening2d
        boundary = 'Block0_Block1'
        pressure_plus = p
        pressure_minus = p
        react_x = resid_primary_x
        react_y = resid_primary_y
        jacob_x = jacob_primary_x
        jacob_y = jacob_primary_y
        react_pressure_x = resid_pressure_x
        react_pressure_y = resid_pressure_y 
        jacob_pressure_x = jacob_pressure_x
        jacob_pressure_y = jacob_pressure_y 
        react_damp_x = resid_damping_x
        react_damp_y = resid_damping_y
        jacob_damp_x = jacob_damping_x
        jacob_damp_y = jacob_damping_y
        nodal_area = nodal_area
        fluid_disp_x = fluid_disp_x
        fluid_disp_y = fluid_disp_y
        permeability_type = 'impermeable'
    [../]
    
[]

[BCs]
    [./fault_p]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = Block0_Block1
        value = 0.0
    [../]
    [./fault_n]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = Block1_Block0
        value = 0.0
    [../]
    [./left_vf]
        type = DirichletBC
        variable = fluid_vel_x
        boundary = left
        value = 0.0
    [../]
    [./right_vf]
        type = DirichletBC
        variable = fluid_vel_x
        boundary = right
        value = 0.0
    [../]
    [./top_vf]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = top
        value = 0.0
    [../]
    [./bot_vf]
        type = DirichletBC
        variable = fluid_vel_y
        boundary = bottom
        value = 0.0
    [../]
  ##non-reflecting bc
    [./dashpot_top_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = top
    []
    [./dashpot_top_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = top
    []
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = bottom
    []
    [./dashpot_left_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = left
    []
    [./dashpot_right_x]
        type = NonReflectDashpotBC
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        p_wave_speed = 4003
        shear_wave_speed = 2444
        boundary = right
    []
[]

[UserObjects]
    [./nodal_area]
      type = NodalArea
      variable = nodal_area
      boundary = Block0_Block1
      execute_on = 'initial TIMESTEP_BEGIN'
    [../]
[]


[Preconditioning]
    [./smp]
        type = SMP
        full = true
       # petsc_options_iname = '-snes_test_jacobian  -pc_factor_mat_solver_type -pc_factor_shift_type'
        #petsc_options_value = 'lu umfpack NONZERO'
    [../]
[]

[Executioner]
    type = Transient
    dt = 0.0002
    end_time = 3.6
    automatic_scaling = true
    [TimeIntegrator]
         type = CentralDifference
    []
    
[]

[Outputs]
    exodus = true
    time_step_interval = 50
[]

