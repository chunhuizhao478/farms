[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2000 
  xmax = 1000
  elem_type = EDGE3
[]

[GlobalParams]
  displacements = 'disp_x'
  use_displaced_mesh = false
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]


[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'p'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[] 

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0.0
    bulk_modulus = 2.5e9
    viscosity = 0.001
    density0 = 1000
  []
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
  []
  [p]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Solid mechanics kernels
  [solid_inertia]
        type = InertialForce
        variable = disp_x
        use_displaced_mesh = false
  []
  [stress_divergence]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        use_displaced_mesh = false
  []
  [viscoelastic_x]
        type = ViscoelasticStressKernel
        variable = disp_x
        component = 0
        kelvin_voigt_viscosity = 10 
        shear_modulus = 8e9      
  []
  [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 1
        variable = disp_x
        component = 0
  []
  [mass0]
        type = PorousFlowFullySaturatedMassTimeDerivative
        biot_coefficient = 1
        coupling_type = HydroMechanical
        variable = p
        multiply_by_density = false
   []
   [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p
        multiply_by_density = false
   []
[]

[Functions]
    [./bc_func]
        type = ParsedFunction
        expression = 'if(t>0,-1e6,0)'
    [../]
[]

[BCs]
  [left_fixed]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [right_pressure]
    type = DirichletBC
    variable = p
    boundary = right
    value = 0
  []
  [left_pressure]
    type = NeumannBC
    variable = p
    boundary = left
    value = 0
  []
  [step_load]
    type = FunctionNeumannBC
    variable = disp_x
    boundary = right
    function = bc_func
  []
[]

[Materials]
[temperature]
        type = PorousFlowTemperature
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        youngs_modulus = 20e9
        poissons_ratio = 0.25
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
        prop_values = 2350
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
        porosity = 0.1
    []
    [biot_modulus]
        type = PorousFlowConstantBiotModulus
        biot_coefficient = 1
        solid_bulk_compliance = 0  # 1/Ks
        fluid_bulk_modulus = 2.5e9        # Kf
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '10e-14 0 0   0 10e-14 0   0 0 10e-14'
    []
[]

[Executioner]
    type = Transient
    solve_type = PJFNK
    start_time = 0.0
    end_time = 3.0
    dt = 1e-4  
    [TimeIntegrator]
        type = NewmarkBeta  
    []
[]

[Outputs]
  exodus = true
  interval = 100
[]
