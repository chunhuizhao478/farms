[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../mesh/2dborehole.msh'
    []
    displacements = 'disp_x disp_y'
[]
    
[GlobalParams]
    displacements = 'disp_x disp_y'
    PorousFlowDictator = dictator #All porous modules must contain 
[]
  
[Variables]
    #displacement components
    [disp_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
    []
    #pore pressure
    [pp]
        order = FIRST
        family = LAGRANGE  
    []
[]
  
[AuxVariables]
    [vel_x]
        order = FIRST
        family = LAGRANGE
    []
    [accel_x]
        order = FIRST
        family = LAGRANGE
    []
    [vel_y]
        order = FIRST
        family = LAGRANGE
    []
    [accel_y]
        order = FIRST
        family = LAGRANGE
    []
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []  
    #
    [./biot_modulus_aux]
        order = CONSTANT
        family = MONOMIAL      
    []   
    [./effective_perm00_aux]
        order = CONSTANT
        family = MONOMIAL      
    [] 
    [./effective_perm11_aux]
        order = CONSTANT
        family = MONOMIAL      
    []
    [./effective_perm22_aux]
        order = CONSTANT
        family = MONOMIAL      
    [] 
    #
    [./xi]
        order = CONSTANT
        family = MONOMIAL         
    []      
[]
  
[Kernels]
    #effective stress tensor
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_x
        component = 0
        use_displaced_mesh = false
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y'
        variable = disp_y
        component = 1
        use_displaced_mesh = false
    []
    #pressure coupling on stress tensor
    [poro_x]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_x
        component = 0
    []
    [poro_y]
        type = PorousFlowEffectiveStressCoupling
        biot_coefficient = 0.5
        variable = disp_y
        component = 1
    []
    #flux * grad(test)
    [flux]
        type = PorousFlowFullySaturatedDarcyBase
        multiply_by_density = false
        variable = pp
        gravity = '0 0 0'
    []
[]
    
[BCs]
    #fix bottom boundary
    [fix_bottom_x]
        type = DirichletBC
        variable = disp_x
        boundary = 7
        value = 0
    []
    [fix_bottom_y]
        type = DirichletBC
        variable = disp_y
        boundary = 7
        value = 0
    []
    [out_right]
        type = PorousFlowPiecewiseLinearSink
        variable = pp
        boundary = 'right'
        pt_vals = '-1e9 1e9' # x coordinates defining g
        multipliers = '-1e9 1e9' # y coordinates defining g
        PT_shift = 3.4e6   # BC pressure
        flux_function = 1E-6 # Variable C
        fluid_phase = 0
        save_in = fluxes_in
    []
    #applied load on top boundary
    [applied_top_z_dispload]
        type = FunctionDirichletBC
        variable = disp_y
        boundary = 6
        function = applied_load_top
    []
    #applied confining pressure on the outer boundary
    [./Pressure]
        [./outer_boundary]
          boundary = 4
          factor = 20.6e6
          displacements = 'disp_x disp_y disp_z'
        [../]
        [./inner_boundary]
          boundary = 5
          factor = 3.4e6
          displacements = 'disp_x disp_y disp_z'
        [../]
    []
[]

[Functions]
    [applied_load_top]
        type = ParsedFunction
        expression = '-2.6477e-5'
    []
[]

[AuxKernels]
    [get_xi]
        type = CompXi3D
        variable = xi
    []
[]
    
[Materials]
    [./elasticity_tensor]
      type = ComputeIsotropicElasticityTensor
      youngs_modulus = 48.5e9
      poissons_ratio = 0.22
    [../]
    [./elastic_stress]
      type = ComputeLinearElasticStress
      outputs = exodus
    []
    [strain]
      type = ComputeSmallStrain
      displacements = 'disp_x disp_y'
      outputs = exodus
    []
    #
    [temperature]
      type = PorousFlowTemperature
    []
    [eff_fluid_pressure_qp]
      type = PorousFlowEffectiveFluidPressure
    []
    #compute volumetric strain and its rate
    [vol_strain]
      type = PorousFlowVolumetricStrain
      outputs = exodus
    []
    #This Material is used for the fully saturated single-phase situation "
    #"where porepressure is the primary variable", saturation = 1.0
    [ppss]
      type = PorousFlow1PhaseFullySaturated
      porepressure = pp
    []
    #List of variables that represent the mass fractions.
    #If no "variables are provided then num_phases=1=num_components."
    [massfrac]
      type = PorousFlowMassFraction
    []
    #compute porosity
    [porosity]
      type = PorousFlowPorosityConst # only the initial value of this is ever used
      porosity = 0.008 #slide
    []
    #comopute permeability
    [permeability]
      type = PorousFlowPermeabilityConst
      permeability = '1E-20 0 0 0 1E-20 0 0 0 1E-20' #slide
    []
    #compute biot modulus
    [biot_modulus]
      type = PorousFlowConstantBiotModulus
      biot_coefficient = 0.5 #paper
      solid_bulk_compliance = 3.46e-11 #calculated
      fluid_bulk_modulus = 2E+09
    []
    #Compute density and viscosity
    [simple_fluid_qp]
      type = PorousFlowSingleComponentFluid
      fp = the_simple_fluid
      phase = 0
    []
[]
  
#provide fluid properties for porous flow 
[FluidProperties]
    [the_simple_fluid]
        type = SimpleFluidProperties
        bulk_modulus = 2.2e+9
        density0 = 1000
        thermal_expansion = 0
        viscosity = 1e-3
    []
[]
  
[Preconditioning]
    [smp]
        type = SMP
        full = true
    []
[]
  
#this user object must contain for porous flow
[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pp'
        number_fluid_phases = 1
        number_fluid_components = 1
    []
[]

[Executioner]
    type = Steady
    solve_type = Newton
    petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    petsc_options_value = '101                asm      lu'
    automatic_scaling = true
    line_search = 'none'
    # num_steps = 1
    l_max_its = 100
    nl_max_its = 10
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-8
    l_tol = 1e-5
[]
  
[Outputs]
    exodus = true
[]