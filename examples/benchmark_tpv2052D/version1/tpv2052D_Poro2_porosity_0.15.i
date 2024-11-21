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
        #elem_type = QUAD9
        
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
    displacements = 'disp_x disp_y' 
    fluid_vel = 'fluid_vel_x fluid_vel_y'
    porepressure = 'p'
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
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
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_x]
    []
    [./vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_y]
    []
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damping_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damping_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_pressure_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_pressure_y]
        order = FIRST
        family = LAGRANGE
    [../]
[]

[AuxKernels]
    [velocity_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
    []
   # [velocity_y]
   ##     type = CompVarRate
   #     variable = vel_y
   ##     coupled = disp_y
    #[]
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
        save_in = 'resid_x'
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
        save_in = 'resid_y'
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
    [../]
    [./poromechskeletoncoupling_y]
        type = PoroMechanicsCoupling
        variable = disp_y
        porepressure = p
        component = 1
        save_in = 'resid_pressure_y'
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
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = 'disp_y'
        component = '1'
        save_in = 'resid_damping_y'
    []
[]

[Materials]
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 20.12e9
        shear_modulus = 18.03e9
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
        prop_values = 2402
    []
    [./rhof]
        type = GenericConstantMaterial
        prop_names = rhof
        prop_values = 1000
    [../]
    [./turtuosity]
        type = GenericConstantMaterial
        prop_names = taut
        prop_values = 2.58
    [../]
    [./porosity]
        type = GenericConstantMaterial
        prop_names = porosity
        prop_values = 0.15
    [../]
    [./hydconductivity]
        type = GenericConstantMaterial
        prop_names = hydconductivity
        prop_values = 3.177893026e-10
    [../]
    [./hydconductivity_layer]
        type = GenericConstantMaterial
        prop_names = hydconductivity_layer
        prop_values = 3.177893026e-10
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.4363
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 1.338876574e10
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
        type = PoroSlipWeakening2dImpermeable
        boundary = 'Block0_Block1'
        pressure_plus = p
        pressure_minus = p
        react_x = resid_x
        react_y = resid_y
        react_pressure_x = resid_pressure_x
        react_pressure_y = resid_pressure_y 
        react_damp_x = resid_damping_x
        react_damp_y = resid_damping_y
        nodal_area = nodal_area
    [../]
[]

[BCs]
  [./fault_p]
    type = FunctionDirichletBC
    variable = fluid_vel_y
    boundary = Block0_Block1
    function = 0.0
  [../]
  [./fault_n]
    type = FunctionDirichletBC
    variable = fluid_vel_y
    boundary = Block1_Block0
    function = 0.0
  [../]
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
   # solve_type = 'PJFNK'
    dt = 0.002
    end_time = 2.4
    #verbose = true
    automatic_scaling = true
    [TimeIntegrator]
         type = CentralDifference
         #type = NewmarkBeta
        #solve_type = lumped
    []
    
[]

[Outputs]
    exodus = true
    time_step_interval = 20
     #print_linear_residuals = true
[]


#[Debug]
#  show_material_props = true
#[]
#[Debug]
#    show_actions = true
#[]

