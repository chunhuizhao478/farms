#solid properties
#----------------------------------------------------#
E = 40e9 # Young's modulus
nu = 0.2 # Poisson's ratio
Gc_const = 100  # critical energy release rate, N * m
ft = 5e6 # tensile strength, Pa
solid_density = 2000 # kg/m^3 
dx_min = 0.25 # minimum mesh size, m
K = '${fparse E/3.0/(1.0-2.0*nu)}'
G = '${fparse E/2.0/(1.0+nu)}'
l =  ${fparse 5 * dx_min} # length scale, m
#'${fparse 3.0/8.0 * E*Gc_const/(ft*ft)}' # AT1 model, N * h, N: number of elements, h: element size -> l = 1.64e-3 m -> this only works for CZM model
Cs = '${fparse sqrt(G/solid_density)}'
Cp = '${fparse sqrt((K + 4.0/3.0 * G)/solid_density)}'
#----------------------------------------------------#

#hydraulic properties
#----------------------------------------------------#
initial_pore_pressure = 0
fluid_density = 1000
biot_coefficient = 1.0
fluid_bulk_modulus = 2.24e+9
viscosity = 1e-6
porosity = 0.137 #(determined from given biot modulus = 10 GPa)
solid_bulk_modulus_compliance = 4.5e-11
intrinsic_permeability = 1e-15 # m^2

##exponential permeability model
# coeff_b = 10 # coefficient for the exponential function in the effective permeability

##darcy-poiseuille permeability model: ultimate crack opening width
wc = ${fparse 2 * Gc_const / ft } # m
perm_exponent = 50 # exponent for the Darcy-Poiseuille model for the effective permeability
#----------------------------------------------------#

#finite element properties
#----------------------------------------------------#
newmark_beta = 0.25
newmark_gamma = 0.5
hht_alpha = 0
#----------------------------------------------------#

#fieldscale small: dx = 1e-3 < l = 1.64e-3, 3x adaptivity levels

# [Adaptivity]
#   max_h_level = 5
#   marker = 'combo'
#   cycles_per_step = 1
#   [Markers]
#       [./combo]
#         type = FarmsComboMarker
#         markers = 'damage_marker strain_energy_marker'
#         meshsize_marker = 'meshsize_marker'
#       [../]
#       [damage_marker]
#         type = ValueThresholdMarker
#         variable = d
#         refine = 0.01
#       []
#       [strain_energy_marker]
#         type = ValueThresholdMarker
#         variable = psie_active
#         refine = '${fparse 1.0*3/8*Gc_const/l}'
#       []   
#       # if mesh_size > dxmin, refine
#       # if mesh_size < dxmin/100, coarsen (which never happens)
#       # otherwise, do nothing
#       [meshsize_marker]
#         type = ValueThresholdMarker
#         variable = mesh_size
#         refine = '${dx_min}'
#         coarsen = '${fparse dx_min/100}'
#         third_state = DO_NOTHING
#       [] 
#   []
# []

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = fracture.i
    cli_args = 'Gc_const=${Gc_const};l=${l};dx_min=${dx_min}'
    execute_on = 'TIMESTEP_END'
    clone_parent_mesh = true
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    from_multi_app = 'fracture'
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    to_multi_app = 'fracture'
    variable = 'psie_active mesh_size'
    source_variable = 'psie_active mesh_size'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator #All porous modules must contain
[]

[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 160
    ny = 160
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 40
    elem_type = QUAD4
  []
  [./damage_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '0 20 0'
    top_right = '4.0 20.25 0'
  [../]
  [./sidesets]
    input = damage_block
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0'
    new_boundary = 'left right bottom top'
  [] 
  displacements = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE  
    scaling = 1e-6
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE  
    scaling = 1e-6
  []
  [pp]
    order = FIRST
    family = LAGRANGE  
  []
[]

[AuxVariables]
  [fy]
  []
  [d]
    family = LAGRANGE
    order = FIRST
  []
  #err measurement of active strain energy
  [eng_err]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [accel_x]
    family = LAGRANGE
    order = FIRST
  []
  [accel_y]
    family = LAGRANGE
    order = FIRST
  []
  #
  [pulse_load_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  #
  [mesh_size]
    family = MONOMIAL
    order = CONSTANT
  []
  #
  [effective_perm00_aux]
    family = MONOMIAL
    order = FIRST
  []
  [effective_perm11_aux]
    family = MONOMIAL
    order = FIRST
  []
  [effective_perm01_aux]
    family = MONOMIAL
    order = FIRST
  []
[]

[AuxKernels]
  #
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  #
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'TIMESTEP_END'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'TIMESTEP_END'
  []
  # #get pulse load aux
  # [get_pulse_load_aux]
  #   type = FunctionAux 
  #   variable = pulse_load_aux
  #   function = func_tri_pulse
  #   execute_on = timestep_end
  # []
  #mesh size aux
  [./max]
    type = ElementLengthAux
    variable = mesh_size
    method = max
    execute_on = TIMESTEP_BEGIN
  [../]
  ### PorousFlow Aux ###
  #effective permeability
  [effective_permeability_00]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 0
    variable = effective_perm00_aux
  []
  [effective_permeability_11]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 1
    column = 1
    variable = effective_perm11_aux
  []
  [effective_permeability_01]
    type = MaterialRealTensorValueAux
    property = effective_perm
    row = 0
    column = 1
    variable = effective_perm01_aux
  []
[]

[Functions]
[]

[Kernels]
  #solid
  [inertia_x]
      type = InertialForce
      variable = disp_x
      acceleration = accel_x
      velocity = vel_x
      beta = 0.25
      gamma = 0.5
      eta = 0
  []
  [inertia_y]
      type = InertialForce
      variable = disp_y
      acceleration = accel_y
      velocity = vel_y
      beta = 0.25
      gamma = 0.5
      eta = 0
  []
  [dispkernel_x]
      type = StressDivergenceTensors
      variable = disp_x
      component = 0
  []
  [dispkernel_y]
      type = StressDivergenceTensors
      variable = disp_y
      component = 1
  []
  #pressure coupling on stress tensor
  [poro_x]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_x
      component = 0
  []
  [poro_y]
      type = PorousFlowEffectiveStressCoupling
      biot_coefficient = ${biot_coefficient}
      variable = disp_y
      component = 1
  []
  #alpha * volumetric strain rate * test + 1 / biot modulus * pressure rate * test
  [mass0]
      type = PorousFlowFullySaturatedMassTimeDerivative
      biot_coefficient = ${biot_coefficient}
      coupling_type = HydroMechanical
      multiply_by_density = false
      variable = pp
  []
  #flux * grad(test)
  [flux]
      type = PorousFlowFullySaturatedDarcyBase
      variable = pp
      multiply_by_density = false
      gravity = '0 0 0'
  []  
[]

[BCs]
  #fix top displacements
  [fix_top_x]
    type = DirichletBC
    variable = disp_x
    boundary = top
    value = 0
  []
  [fix_top_y]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0
  []
  #fix bottom displacements
  [fix_bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [fix_bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  #fix left displacement x
  [fix_left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  #fix right displacements
  [fix_right_x]
    type = DirichletBC
    variable = disp_y
    boundary = right
    value = 0
  []
  [fix_right_y]
    type = DirichletBC
    variable = disp_y
    boundary = right
    value = 0
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${E}
    poissons_ratio = ${nu}
  [../]
  [strain]
    type = ComputeSmallStrain
  []
  [bulk]
    type = GenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [elasticity]
    type = NDSmallDeformationIsotropicElasticity
    # material property names
    ##---------------------------------------------##
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    strain_energy_density = psie
    strain_energy_density_active = psie_active
    strain_energy_density_derivative = dpsie_dd
    degradation_function = g
    degradation_function_derivative = dg_dd
    degradation_function_second_derivative = d2g_dd2
    ##---------------------------------------------##
    # decomposition type
    ##---------------------------------------------##
    decomposition = SPECTRAL
    ##---------------------------------------------##
    # model type
    ##---------------------------------------------##
    model_type = AT1
    ##---------------------------------------------##
    # constants
    ##---------------------------------------------##
    eta = 1e-6
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
    ##---------------------------------------------##
    # porous flow coupling
    ##---------------------------------------------##
    porous_flow_coupling = true
    ##-----darcy_poiseuille_permeability_model-----##
    darcy_poiseuille_permeability_model = true
    intrinsic_permeability = ${intrinsic_permeability}
    wc = ${wc}
    perm_exponent = ${perm_exponent}
    ##---------------------------------------------##
  []
  [stress]
    type = NDComputeSmallDeformationStress ###
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  #solid properties
  ##-------------------------------------------------------------------------##
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${solid_density}
  []
  #define initial bulk modulus material property
  #check with youngs_modulus = 50e9, poissons_ratio = 0.373
  [solid_bulk_modulus_compliance]
    type = GenericConstantMaterial
    prop_names = solid_bulk_modulus_compliance
    prop_values = ${solid_bulk_modulus_compliance}
  []
  ##-------------------------------------------------------------------------##
  #porous flow related properties
  ##-------------------------------------------------------------------------##
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
    porosity = ${porosity}
  []
  #comopute permeability
  [permeability] #take effective_perm
    type = ElkPorousFlowPermeabilityDamaged
  []
  # #compute biot modulus #include damaged solid compliance
  # [biot_modulus]
  #   type = ElkPorousFlowDamagedBiotModulus
  #   biot_coefficient = ${biot_coefficient}
  #   solid_bulk_compliance = ${solid_bulk_modulus_compliance}
  #   fluid_bulk_modulus = ${fluid_bulk_modulus}
  # []
  ##----------------------------------------------------------##
  #compute permeability
  # [permeability_constant]
  #     type = PorousFlowPermeabilityConst
  #     permeability = ${permeability}
  # []
  #compute biot modulus
  [biot_modulus_constant]
      type = PorousFlowConstantBiotModulus
      biot_coefficient = ${biot_coefficient}
      solid_bulk_compliance = ${solid_bulk_modulus_compliance}
      fluid_bulk_modulus = ${fluid_bulk_modulus}
  []  
  ##----------------------------------------------------------##
  #Compute density and viscosity
  [simple_fluid_qp]
    type = PorousFlowSingleComponentFluid
    fp = the_simple_fluid
    phase = 0
  []
  #define relative permeability as 1 (used in PorousFlowDarcyVelocityComponent)
  [relperm]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
    kr = 1
  []
[]

#provide fluid properties for porous flow 
[FluidProperties]
  [the_simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = ${fluid_bulk_modulus}
    density0 = ${fluid_density}
    thermal_expansion = 0
    viscosity = ${viscosity}
  []
[]

#this user object must contain for porous flow
[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp disp_x disp_y'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
  type = Transient

  solve_type = NEWTON

  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist                 '

  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -ksp_gmres_restart'
  petsc_options_value = ' lu       mumps       100'

  # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -ksp_initial_guess_nonzero'
  # petsc_options_value = 'gmres     hypre  boomeramg True'

  # automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 30

  dt = 50e-6
  end_time = 100

  fixed_point_max_its = 10
  accept_on_max_fixed_point_iteration = false
  fixed_point_rel_tol = 1e-8
  fixed_point_abs_tol = 1e-10

  [./TimeIntegrator]
    type = NewmarkBeta
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
  [../]
[]

[DiracKernels]
  [sink1]
    type = PorousFlowSquarePulsePointSource
    start_time = 0
    end_time = 100
    point = '0 20.125 0'
    mass_flux = 2 # kg/s
    variable = pp
  []
[]

[Outputs]
  exodus = true
  time_step_interval = 1
  print_linear_residuals = false
  csv = true
  [checkpoint]
      type = Checkpoint
      time_step_interval = 20
      num_files = 2
  []
[]