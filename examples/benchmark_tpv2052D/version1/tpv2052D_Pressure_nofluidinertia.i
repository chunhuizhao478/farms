# Verification of Benchmark Problem TPV205-2D from the SCEC Dynamic Rupture Validation exercises
# Reference:
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126.

[Mesh]

    #second_order = true
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
    PorousFlowDictator = dictator
    q = 0.1
    Dc = 0.4
    T2_o = 120e6
    elem_size = 100
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
        order = SECOND
        family = LAGRANGE
    [../]
    [./disp_y]
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
    [./interface_p_plus]
        order = FIRST
        family = LAGRANGE
    []
    [./interface_p_minus]
        order = FIRST
        family = LAGRANGE
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
    [pressure_global_plus]
        type = ProjectionAux
        variable = interface_p_plus
        v = porepressure
        boundary = 'Block1_Block0'
    []
    [pressure_global_minus]
        type = ProjectionAux
        variable = interface_p_minus
        v = porepressure
        boundary = 'Block0_Block1'
    []   
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
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.567
        variable = disp_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.567
        variable = disp_y
        component = 1
    []
    [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        biot_coefficient = 0.567
        coupling_type = HydroMechanical
        variable = porepressure
        multiply_by_density = false
    []
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = porepressure
        gravity = '0 0 0'
        multiply_by_density = false
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
[]

[Materials]
    [temperature]
        type = PorousFlowTemperature
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        bulk_modulus = 15.46e9
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
    [eff_fluid_pressure_qp]
        type = PorousFlowEffectiveFluidPressure
    []
    [vol_strain]
        type = PorousFlowVolumetricStrain
    []
    [ppss]
        type = PorousFlow1PhaseFullySaturated
        porepressure = porepressure
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
        porosity = 0.2
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 0.567
        solid_bulk_compliance = 6.4683053e-11
        fluid_bulk_modulus = 2.25e9
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '1.128052989e-12 0 0   0 1.128052989e-12 0   0 0 1.128052989e-12'
    []
    [./czm_stress_derivative]
        type = StressDerivative2
        boundary = 'Block0_Block1'
    [../]
    [./czm_mat]
        type = PoroSlipWeakeningFriction2dNoInertia
        boundary = 'Block0_Block1'
        pressure_plus = interface_p_plus
        pressure_minus = interface_p_minus
        nodal_area = nodal_area
    [../]
[]

[BCs]
  [./fault_p]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = Block0_Block1
    function = 0.0
  [../]
  [./fault_n]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = Block1_Block0
    function = 0.0
  [../]
  [./flux_top]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = top
    function = 0.0
  [../]
  [./flux_bot]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = bottom
    function = 0.0
  [../]
  [./flux_left]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = left
    function = 0.0
  [../]
  [./flux_right]
    type = FunctionNeumannBC
    variable = porepressure
    boundary = right
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
    #solve_type = 'PJFNK'
    dt = 0.002
    end_time = 4
    #verbose = true
    #automatic_scaling = true
    [TimeIntegrator]
        type = CentralDifference
        #type = NewmarkBeta
        #solve_type = lumped
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

