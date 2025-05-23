# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

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
    [./vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_x]
    []
    [./vel_y]
    []
    [./accel_y]
    []
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_primary_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_primary_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_secondary_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_secondary_y]
        order = FIRST
        family = LAGRANGE
    [../]
    [./interface_p_plus]
        order = FIRST
        family = LAGRANGE
    []
    [./interface_p_minus]
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
[]

[AuxKernels]
    [velocity_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
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
    [pressure_global_plus]
        type = ProjectionAux
        variable = interface_p_plus
        v = p
        boundary = 'Block1_Block0'
    []
    [pressure_global_minus]
        type = ProjectionAux
        variable = interface_p_minus
        v = p
        boundary = 'Block0_Block1'
    []   
[]

[Actions/PoroCohesiveZoneAction]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y jump_vel_x jump_vel_y normal_traction tangent_traction normal_jump tangent_jump'
        save_in_master = 'resid_primary_x resid_primary_y'
        save_in_slave  = 'resid_secondary_x  resid_secondary_y'      
    [../]
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
    [./czm_stress_derivative]
        type = StressDerivative2
    #    boundary = 'Block0_Block1'
    [../]
    [./czm_mat]
        type = PermeableSlipWeakening2dv2
        boundary = 'Block0_Block1'
        pressure_plus = interface_p_plus
        pressure_minus = interface_p_minus
        react_x = resid_primary_x
        react_y = resid_primary_y
        react_neighbor_x = resid_secondary_x
        react_neighbor_y = resid_secondary_y
        fluid_disp_x = fluid_disp_slipweakening_x
        fluid_disp_y = fluid_disp_slipweakening_y
        nodal_area = nodal_area
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


[BCs]
  [./fault_p]
    type = FunctionDirichletBC
    variable = p
    boundary = Block0_Block1
    function = 0.0
  [../]
  [./fault_n]
    type = FunctionDirichletBC
    variable = p
    boundary = Block1_Block0
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
    dt = 0.0025
    end_time = 4
    #verbose = true
    automatic_scaling = true
    [TimeIntegrator]
        type = CentralDifference
        #type = NewmarkBeta
       # solve_type = lumped
    []
    
[]

[Outputs]
    exodus = true
    time_step_interval = 10
[]

#[Debug]
#  show_material_props = true
#[]
#[Debug]
#    show_actions = true
#[]

