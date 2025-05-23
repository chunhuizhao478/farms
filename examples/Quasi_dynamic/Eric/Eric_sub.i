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
    displacements = 'disp_sub_x disp_sub_y'
    PorousFlowDictator = dictator
    
    ##radiation damping parameters
    shear_modulus = 30e9
    shear_wave_velocity = 3400  # m/s
    
    ##rate-and-state coefficients
    V_o = 1e-6  # Reference slip rate
    Vini = 6.7e-7  # Initial slip rate
    f_o = 0.6   # Initial friction coefficient
    a = 0.015    # Direct effect parameter
    b = 0.02    # State variable evolution parameter
    L = 1e-5    # State variable characteristic distance
    statevarini = 0.606 # Initial state variable 

    ##Enhanced weakening parameters
    enhanced_weakening = false
    f_w = 0.4    # Weakened friction coefficient
    V_w = 1e-3   # Weakening slip rate

    ##initial normal traction (Pa)
    background_normal_stress = 50e6

    ##initial shear traction (Pa)
    background_tangential_stress = 30e6
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
    [p_sub]
        family = LAGRANGE
        order = FIRST
    []
[]

[AuxVariables]
    # Velocity calculations
    [vel_sub_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_sub_y]
        order = FIRST
        family = LAGRANGE
    []
    
    # Traction and input variables
    [traction_tangential]
        order = FIRST
        family = LAGRANGE
    []
    [traction_normal]
        order = FIRST
        family = LAGRANGE
    []
    [p_fault]
        order = FIRST
        family = LAGRANGE
    []
    
    #sliprate
    [sliprate_sub]
        order = FIRST
        family = MONOMIAL
    []
    #slip
    [slip_x]
        order = FIRST
        family = MONOMIAL
    []
    [slip_y]
        order = FIRST
        family = MONOMIAL
    []
    #statevar
    [statevar_sub]
        order = FIRST
        family = MONOMIAL
    []
    #statevar_dot
    [statevar_dot_sub]
        order = FIRST
        family = MONOMIAL
    []
    [min_elem_size]
         order = CONSTANT
        family = MONOMIAL
    []
[]

[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'p_sub disp_sub_x disp_sub_y'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[]

[AuxKernels]
    # Calculate velocity
    [Vel_sub_x]
        type = CompVarRate
        variable = vel_sub_x
        coupled = disp_sub_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_sub_y]
        type = CompVarRate
        variable = vel_sub_y
        coupled = disp_sub_y
        execute_on = 'TIMESTEP_END'
    []

    # Output material properties
    [output_sliprate]
        type = MaterialRealAux
        property = sliprate
        variable = sliprate_sub
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_slip]
        type = MaterialRealAux
        property = slip
        variable = slip_x
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_statevar_dot]
        type = MaterialRealAux
        property = statevar_dot
        variable = statevar_dot_sub
        boundary = 'Block0_Block1'
        execute_on = 'TIMESTEP_END'
    []
    [output_statevar]
        type = MaterialRealAux
        property = statevar
        variable = statevar_sub
        boundary = 'Block0_Block1'
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
            displacements = 'disp_sub_x disp_sub_y'
            planar_formulation = PLANE_STRAIN
            generate_output = 'stress_xx stress_yy stress_xy'
        [../]
        [../]
    [../]
[]

[Kernels]
    # Poroelastic coupling
    [poro_sub_x]
        type = PorousFlowEffectiveStressCoupling
        variable = disp_sub_x
        component = 0
        biot_coefficient = 0.8145
    []
    [poro_sub_y]
        type = PorousFlowEffectiveStressCoupling
        variable = disp_sub_y
        component = 1
        biot_coefficient = 0.8145
    []
    
    # Mass conservation for poroelasticity
    [mass_sub]
        type = PorousFlowFullySaturatedMassTimeDerivative
        variable = p_sub
        biot_coefficient = 0.8145
    []
    
    # Darcy flow
    [flux_sub]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p_sub
        gravity = '0 0 0'
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
        porosity = 0.1
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.8145
        solid_bulk_compliance = 2.5e-11
        fluid_bulk_modulus = 2e9
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '0.025e-9 0 0 0 0.025e-9 0 0 0 0.025e-9'
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
        porepressure = p_sub
    []
    [massfrac]
        type = PorousFlowMassFraction
    []
    [simple_fluid_qp]
        type = PorousFlowSingleComponentFluid
        fp = simple_fluid
        phase = 0
    [] 
    # RSF Material
    [./czm_mat]
        type = QuasiDynamicRSF
        normal_traction = traction_normal
        tangential_traction = traction_tangential
        pore_pressure = p_fault
        boundary = 'Block0_Block1'
        output_properties = 'sliprate slip statevar statevar_dot'
    [../]
[]


[Postprocessors]
  [max_sliprate_value]
    type = ElementExtremeValue
    variable = sliprate_sub  
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
    [../]
[]

[Executioner]
    type = Transient
    [TimeStepper]
        type = AdaptiveTimeStepCalculator
        cp = 5500                # P-wave velocity
        a_prem = 0.1             # Proportionality index for stable time step
        shear_modulus = 30e9     # Shear modulus
        permeability = 1e-9   # Permeability value (k_py)
        biot_modulus = 1   # Inverse of Biot modulus (c_o)
        a_o = 0.015              # Direct effect parameter
        b_o = 0.02               # State evolution parameter
        L = 1e-5                # Characteristic slip distance
        normal_stress = 50e6    # Normal stress on the fault
        zeta_max = 0.5           # Maximum allowed value for zeta
        dx_min = min_elem_size_pp              # Minimum element size
        max_slip_rate = max_sliprate_value  # Reference to the postprocessor that computes max slip rate
    []
    end_time = 8.0
[]

[Outputs]
    exodus = true
    interval = 2
[]
