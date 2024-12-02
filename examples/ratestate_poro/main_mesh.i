[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = './Inclined_fault_with_injection.msh'
    []
    [subdomain1]
        input = msh
        type = SubdomainBoundingBoxGenerator
        bottom_left = '-5000 -5000 0'
        top_right = '5000 5000 0'
        block_id = 0
    []
    [./inner_block]
        type = ParsedSubdomainMeshGenerator
        input = subdomain1
        combinatorial_geometry = 'x > -1732.050808 & x < 1732.050808'
        block_id = 1
    []
    [./fault_block_upper]
        type = ParsedSubdomainMeshGenerator
        input = inner_block
        combinatorial_geometry = 'y > (-0.5773502692 * x) & x > -1732.050808 & x < 1732.050808'
        block_id = 2
    []
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = fault_block_upper
        split_interface = true
        add_interface_on_two_sides = true
        block_pairs = '1 2'
        show_info = true
    [../]
[]

[GlobalParams]
    displacements = 'disp_x disp_y' 
    PorousFlowDictator = dictator
    gravity = '0 0 0'
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
  [stress_xx]
    order = SECOND
    family = MONOMIAL
  []
  [stress_yy]
    order = SECOND
    family = MONOMIAL
  []
  [shear]
    order = SECOND
    family = MONOMIAL
  []
[]

[AuxKernels]
  [s_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = 'initial timestep_end'  
[]
  [s_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = shear
    index_i = 1
    index_j = 0
    execute_on = 'initial timestep_end' 
  []
  [s_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = 'initial timestep_end' 
  []
[]

[Kernels]
    [./stressdiv_x]
        type = StressDivergenceTensors
        variable = disp_x
        component = 0
        use_displaced_mesh = false    
    [../]
    [./stressdiv_y]
        type = StressDivergenceTensors
        variable = disp_y
        component = 1
        use_displaced_mesh = false    
    [../]
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.95
        variable = disp_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.95
        variable = disp_y
        component = 1
    []
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        variable = p
        multiply_by_density = false
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
        biot_coefficient = 0.95
        solid_bulk_compliance = 3.75e-11  # 1/Ks
        fluid_bulk_modulus = 2.5e9        # Kf
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '10e-14 0 0   0 10e-14 0   0 0 10e-14'
    []
[]

[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'p'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[]

[BCs]
    #assign displacement boundary condition
    [./matchval_primary_x]
        type = DirichletBC
        variable = disp_x
        boundary = 'domain_Block2'
        value = 0
    []
    [./matchval_secondary_x]
        type = DirichletBC
        variable = disp_x
        boundary = 'Block2_domain'
        value = 0
    []
    [./matchval_primary_y]
        type = DirichletBC
        variable = disp_y
        boundary = 'domain_Block2'
        value = 0
    []
    [./matchval_secondary_y]
        type = DirichletBC
        variable = disp_y
        boundary = 'Block2_domain'
        value = 0
    []
    [./pressure_primary_y]
        type = NeumannBC
        variable = p
        boundary = 'domain_Block2'
        value = 0
    []
    [./pressure_secondary_y]
        type = NeumannBC
        variable = p
        boundary = 'Block2_domain'
        value = 0
    []
    # Displacement and tractions boundary conditions
    [./disp_bottom]
        type = DirichletBC
        variable = disp_y
        boundary = '1'
        value = 0
    [../]
    [./disp_left]
        type = DirichletBC
        variable = disp_x
        boundary = 'left'
        value = 0
    [../]
    [./traction_right]
        type = NeumannBC
        variable = disp_x
        boundary = 'right'
        value = -20e6
    [../]
    [./traction_top]
        type = NeumannBC
        variable = disp_y
        boundary = 'top'
        value = -10e6
    [../]
    [./pressure_top]
        type = DirichletBC
        variable = p
        boundary = 'top'
        value = 0.0
    [../]
    [./pressure_right]
        type = DirichletBC
        variable = p
        boundary = 'right'
        value = 0.0
    [../]
    [./flux_bot]
        type = NeumannBC
        variable = p
        boundary = '1'
        value = 0.0
    [../]
    [./flux_left]
        type = NeumannBC
        variable = p
        boundary = 'left'
        value = 0.0
    [../]
[]



[Executioner]
    type = Steady
    automatic_scaling = true

    # Add line search to help convergence
    line_search = 'basic'  # Enable basic line search
    
    # Improve linear solver behavior
    petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type'
    petsc_options_value = '200 hypre boomeramg'
[]


[Outputs]
  exodus = true
[]