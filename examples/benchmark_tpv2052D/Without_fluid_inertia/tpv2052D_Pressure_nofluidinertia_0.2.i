# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  './New_Planar_fault.msh'
    []
    [subdomain1]
        input = msh
        type = SubdomainBoundingBoxGenerator
        bottom_left = '-20000 -20000 0'
        top_right = '10000 10000 0'
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
    PorousFlowDictator = dictator
    q = 0.5
    Dc = 0.4
    T2_o = 120e6
    elem_size = 50
    mu_d = 0.525
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'p'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
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
    [./resid_primary_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_primary_y]
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


[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0.0
    bulk_modulus = 2.25e9
    viscosity = 0.001
    density0 = 1000
  []
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

[Modules/TensorMechanics/CohesiveZoneMaster]
    [./czm_ik]
        boundary = 'Block0_Block1'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y normal_traction tangent_traction normal_jump tangent_jump'
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
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        displacements = 'disp_x disp_y'
        use_displaced_mesh = false
        save_in = 'resid_primary_y' 
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
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.4092
        variable = disp_x
        component = 0
        save_in = 'resid_pressure_x'
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.4092
        variable = disp_y
        component = 1
        save_in = 'resid_pressure_y'
    []
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        biot_coefficient = 0.4092
        coupling_type = HydroMechanical
        variable = p
        multiply_by_density = false
    []
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p
        gravity = '0 0 0'
        multiply_by_density = false
    []
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
        lambda = 6.22219e9
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
        prop_values = 1.1430653319e-9
    [../]
    [./hydconductivity_layer]
        type = GenericConstantMaterial
        prop_names = hydconductivity_layer
        prop_values = 1.1430653319e-9
    [../]
    [./biotcoeff]
        type = GenericConstantMaterial
        prop_names = biot_coefficient
        prop_values = 0.5669
    [../]
    [./biotmodulus]
        type = GenericConstantMaterial
        prop_names = biot_modulus
        prop_values = 1.00841e10
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

[Materials]
    [temperature]
        type = PorousFlowTemperature
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 21.09e9
        shear_modulus = 18.9e9
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
        prop_values = 2419
    []
    [eff_fluid_pressure_qp]
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
    [porosity]
        type = PorousFlowPorosityConst # only the initial value of this is ever used
        porosity = 0.14
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.4092
        solid_bulk_compliance = 4.7412329e-11
        fluid_bulk_modulus = 2.25e9
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '2.3e-13 0 0   0 2.3e-13 0   0 0 2.3e-13'
    []
    [./czm_stress_derivative]
        type = StressDerivative2
        boundary = 'Block0_Block1'
    [../]
    [./czm_mat]
        type = PoroSlipWeakeningFriction2dNoInertia
        boundary = 'Block0_Block1'
        pressure_plus = p
        pressure_minus = p
        react_x = resid_primary_x
        react_y = resid_primary_y
        react_pressure_x = resid_pressure_x
        react_pressure_y = resid_pressure_y 
        react_damp_x = resid_damping_x
        react_damp_y = resid_damping_y
        nodal_area = nodal_area
    [../]
[]

[BCs]
  [./fault_p]
    type = FunctionNeumannBC
    variable = p
    boundary = Block0_Block1
    function = 0.0
  [../]
  [./fault_n]
    type = FunctionNeumannBC
    variable = p
    boundary = Block1_Block0
    function = 0.0
  [../]
  [./flux_top]
    type = FunctionNeumannBC
    variable = p
    boundary = top
    function = 0.0
  [../]
  [./flux_bot]
    type = FunctionNeumannBC
    variable = p
    boundary = bottom
    function = 0.0
  [../]
  [./flux_left]
    type = FunctionNeumannBC
    variable = p
    boundary = left
    function = 0.0
  [../]
  [./flux_right]
    type = FunctionNeumannBC
    variable = p
    boundary = right
    function = 0.0
  []
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
    #solve_type = 'PJFNK'
    dt = 0.001
    end_time = 4.8
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
[]

#[Debug]
#  show_material_props = true
#[]
#[Debug]
#    show_actions = true
#[]

