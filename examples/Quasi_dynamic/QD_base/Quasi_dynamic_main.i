[Mesh]
    [./msh]
        type = GeneratedMeshGenerator
        dim = 2
        nx = 200
        ny = 200
        xmin = -4
        xmax = 4
        ymin = -4
        ymax = 4
    [../]
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
        add_interface_on_two_sides = true
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
    PorousFlowDictator = dictator
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
    [lambda_x]
        # Lagrange multiplier for tangential traction
        family = LAGRANGE
        order = FIRST
    []
    [lambda_y]
        # Lagrange multiplier for normal traction
        family = LAGRANGE
        order = FIRST
    []
    [p]
        family = LAGRANGE
        order = FIRST
    []
[]

[AuxVariables]
    # Velocity calculations
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    
    # Variables to receive slip from subapp
    [slip_x]
        family = MONOMIAL
        order = FIRST
    []
    [slip_y]
        family = MONOMIAL
        order = FIRST
    []
    [sliprate]
        family = MONOMIAL
        order = FIRST
    []
    [min_elem_size]
         order = CONSTANT
        family = MONOMIAL
    []
[]

[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'p disp_x disp_y'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[]

[AuxKernels]
    # Calculate velocity
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [min_elem_size_calculator]
        type = ElementLengthAux
        variable = min_elem_size
        method = min
        execute_on = initial
    []
[]

[Modules]
    [./TensorMechanics]
        [./Master]
        [./all]
            strain = SMALL
            displacements = 'disp_x disp_y'
            planar_formulation = PLANE_STRAIN
            generate_output = 'stress_xx stress_yy stress_xy'
        [../]
        [../]
    [../]
[]

[Kernels]
    # Poroelastic coupling
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        variable = disp_x
        component = 0
        biot_coefficient = 0.8145
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        variable = disp_y
        component = 1
        biot_coefficient = 0.6
    []
    
    # Mass conservation for poroelasticity
    [mass]
        type = PorousFlowFullySaturatedMassTimeDerivative
        variable = p
        biot_coefficient = 0.6
    []
    
    # Darcy flow
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p
        gravity = '0 0 0'
    []
    # Kernels for Lagrange multiplier variables
    [lambda_x_kernel]
        type = NullKernel
        variable = lambda_x
    []
    
    [lambda_y_kernel]
        type = NullKernel
        variable = lambda_y
    []
[]


[InterfaceKernels]
     # Normal slip in y-direction
     [y_slip]
         type = FaultSlipInterfaceKernel
         variable = disp_y
         neighbor_var = disp_y
         slip = slip_y
         boundary = 'Block1_Block0'
     []
     
     # Tangential slip in x-direction
     [x_slip]
        type = FaultSlipInterfaceKernel
         variable = disp_x
         neighbor_var = disp_x
         slip = slip_x
         boundary = 'Block1_Block0'
    []
    
    # Normal traction continuity with Lagrange multiplier
    [y_traction]
        type = FaultNormalTractionLagrangeMultiplier
        variable = disp_y
        neighbor_var = disp_y
        lambda_y = lambda_y
        boundary = 'Block1_Block0'
    []
    
    # Shear traction continuity with Lagrange multiplier
    [x_traction]
        type = FaultShearTractionLagrangeMultiplier
        variable = disp_x
        neighbor_var = disp_x
        lambda_x = lambda_x
        boundary = 'Block1_Block0'
    []
[]

[FluidProperties]
    [simple_fluid]
        type = SimpleFluidProperties
        bulk_modulus = 2.2e9
        density0 = 1000
        viscosity = 1e-3
    []
[]

[Materials]
    # Elasticity and mechanics materials
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 20e9
        shear_modulus = 30e9 
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    # Poroelasticity materials
    [porosity]
        type = PorousFlowPorosityConst
        porosity =  0.0274
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.8145
        solid_bulk_compliance = 2.5e-11
        fluid_bulk_modulus = 2e9
    []
    [permeability_lower]
        type = PorousFlowPermeabilityConst
        permeability = '0.025e-9 0 0 0 0.025e-9 0 0 0 0.025e-9'
        block = 1
    []
    [temperature]
        type = PorousFlowTemperature
    []
    [eff_fluid_pressure]
        type = PorousFlowEffectiveFluidPressure
    []
    [vol_strain]
        type = PorousFlowVolumetricStrain
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = p
    []
    [massfrac]
        type = PorousFlowMassFraction
    []
    [simple_fluid_qp]
        type = PorousFlowSingleComponentFluid
        fp = simple_fluid
        phase = 0
    []
[]


[BCs]
    # Fixed boundaries
    [no_x_left]
        type = DirichletBC
        variable = disp_x
        boundary = 'left'
        value = 0.0
    []
    [no_x_right]
        type = DirichletBC
        variable = disp_x
        boundary = 'right'
        value = 0.0
    []
    [no_y_bottom]
        type = DirichletBC
        variable = disp_y
        boundary = 'bottom'
        value = 0.0
    []
    
    # Pressure boundary conditions
    [p_left]
        type = DirichletBC
        variable = p
        boundary = 'left'
        value = 0.0
    []
    [p_right]
        type = DirichletBC
        variable = p
        boundary = 'right'
        value = 0.0
    []
    
    # Applied loading on top boundary
    [top_load]
        type = FunctionDirichletBC
        variable = disp_y
        boundary = 'top'
        function = '-0.001*t'  # Gradual loading
    []
[]


[Postprocessors]
  [max_sliprate_value]
    type = ElementExtremeValue
    variable = sliprate  
    execute_on = 'TIMESTEP_END'
  []
  [min_elem_size_pp]
    type = ElementExtremeValue
    variable = min_elem_size
    value_type = min
    execute_on = 'INITIAL'
  []
[]


[Preconditioning]
    [./smp]
        type = SMP
        full = true
        petsc_options = '-snes_ksp_ew'
        petsc_options_iname = '-ksp_gmres_restart -pc_type'
        petsc_options_value = '100 asm'
    [../]
[]

[Executioner]
    type = Transient
    [TimeStepper]
        type = AdaptiveTimeStepCalculator
        cp = 5500                # P-wave velocity
        a_prem = 0.1             # Proportionality index for stable time step
        shear_modulus = 30e9     # Shear modulus
        permeability = 0.025e9   # Permeability value (k_py)
        biot_modulus = 57.63e9   # Inverse of Biot modulus (c_o)
        a_o = 0.008              # Direct effect parameter
        b_o = 0.012              # State evolution parameter
        L = 0.02                 # Characteristic slip distance
        normal_stress = 50e6    # Normal stress on the fault
        zeta_max = 0.5           # Maximum allowed value for zeta
        dx_min = min_elem_size_pp              # Minimum element size
        max_slip_rate = max_sliprate_value  # Reference to the postprocessor that computes max slip rate
    []
    end_time = 8.0
[]

[Outputs]
    exodus = true
    interval = 10
[]

[MultiApps]
    # SubApp to solve rate-and-state friction
    [./rate_state_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'Quasi_dynamic_sub.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    # Send Lagrange multipliers and pressure to subapp
    [send_tractions]
        type = MultiAppCopyTransfer
        to_multi_app = rate_state_app
        source_variable = 'lambda_x lambda_y p'
        variable = 'traction_tangential traction_normal pore_pressure'
        execute_on = 'TIMESTEP_BEGIN'
    []
    
    # Get slip from subapp
    [get_slip]
        type = MultiAppCopyTransfer
        from_multi_app = rate_state_app
        source_variable = 'slip_x slip_y sliprate_sub'
        variable = 'slip_x slip_y sliprate'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]