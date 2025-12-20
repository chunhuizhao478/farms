#parameters

##mesh parameters
###need to check gmsh file for changing the parameters here
bottom_nodes_coord =' -60000 -60000 -60000;
                      60000 -60000 -60000;
                      60000 60000  -60000;
                     -60000 60000  -60000'

##element size
elem_size = 200 #!!! element size near the fault, need to be consistent with the mesh file

##mesh domain
nonlocal_eqstrain_blocks = '10 12 100 200'

#here we avoid the cross-fault averaging by defining separate averaging blocks
nonlocal_eqstrain_blocks_1 = '10 12'
nonlocal_eqstrain_blocks_2 = '100'
nonlocal_eqstrain_blocks_3 = '200'

#local_eqstrain_blocks = '11'

##main fault parameters
xmin_fault = -22500 #xmin of fault
xmax_fault = 22500 #xmax of fault
zmin_fault = -20000 #zmin of fault
# zmax_fault = 0 #zmax of fault

#nonlocal length applied region along ydir
ymin_fault = -1500
ymax_fault = 1500
nonlocal_averaging_length_scale = 250 #your length scale must be resolved by multiple elements, or xi_out_of_range error will occur, for 400m testing, use 800m
nonlocal_averaging_radius = 500 #for 800m testing, use 1600m

##-------------------------##
##material properties##
density = 2670 #density
lambda_o = 2.1547e10 #first lame constant
shear_modulus_o = 2.975e10 #second lame constant
# Cs = '${fparse shear_modulus_o / density }' #shear wave speed
# Cp = '${fparse (lambda_o + 2 * shear_modulus_o) / density }' #pressure wave speed
##-------------------------##

##Slip weakening parameters##
Dc = 0.8 #characteristic length (m)
q = 0.4 #damping ratio
mu_s = 0.8 #static friction coefficient
mu_d = 0.6 #dynamic friction coefficient
##-------------------------##

##Cohesion parameters##
cohesion_depth = 6500 #cohesion depth (m)
cohesion_slope = 0.01 #cohesion slope (MPa/m)
cohesion_min = 0.4 #minimum cohesion value (MPa)
##---------------------------------------------##

##CDB model parameters##
xi_0 = -0.8 #strain invariants ratio: onset of damage evolution
xi_d = -0.8 #strain invariants ratio: onset of breakage healing

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
use_strain_rate_dependent_Cd = false #use strain rate dependent Cd
m_exponent = 0.8 #strain rate dependent parameters
strain_rate_hat = 1e-4 #strain rate dependent parameters
cd_hat = 10 #strain rate dependent parameters
###

CdCb_multiplier = 100 #multiplier between Cd and Cb
CBH_constant = 0 #coefficient of healing for breakage evolution
C_1 = 0 #coefficient of healing for damage evolution
C_2 = 0.05 #coefficient of healing for damage evolution
beta_width = 0.05 #coefficient gives width of transitional region
C_g = 1e-12 #material parameter: compliance or fluidity of the fine grain granular material
m1 = 10 #coefficient of power law indexes
m2 = 1 #coefficient of power law indexes
chi = 0.8 #energy ratio
##-------------------------##

### poroelastic properties

fluid_bulk_modulus = 2.2e9   # Water bulk modulus (2.2 GPa)
permeability_solid_o = 1e-19 # Initial permeability (1 milli-darcy)
porosity_solid_o = 0.0026 # Initial porosity (8%)
solid_bulk_modulus_g = 57e9 # Granular bulk modulus (50.3 GPa)
solid_bulk_modulus_s = 57e9 # Solid grains bulk modulus (36 GPa - typical for quartz)
permeability_evolution_with_damage = 3
initial_grain_size = 1
ultimate_grain_size = 1
initial_viscosity_fluid = 1e-3
anand_param_go_mat = 0.25
anand_param_eta_cv_mat = 0.01
anand_param_p_mat = 1

##------------------------------------------------------------------##

##initial stress parameters##
#background stress
#fluid_density = 1000
#gravity = 9.8
#bxx = 0.926793
#byy = 1.073206
#bxy = -0.8
##------------------------------------------------------------------##

##tapering parameters##
#use_tapering = true #use tapering to reduce deviatoric stress components at shallow depth
#tapering_depth_A = 15000 #depth at which tapering starts to be applied (m)
#tapering_depth_B = 20000 #depth at which tapering stops to be applied (m)
##------------------------------------------------------------------##

#nucleation parameters
nucl_center_x = -21000 #nucleation center x coordinate
nucl_center_y = 0 #nucleation center y coordinate
nucl_center_z = -10000 #nucleation center y coordinate
r_crit = 3000 #critical distance to hypocenter (m)
Vs = 3340 #shear wave speed (m/s)
t0 = 0.1 #nucleation time (s)
##------------------------------------------------------------------##

##model parameters##
dt = 0.0025 #time step size

end_time = 10.0 #end time for simulation

# num_steps = 40 #end_time or num_steps only one of them is needed
exodus_time_step_interval = 10 #time step interval for output
sample_snapshots_time_step_interval = 400 #time step interval for sample snapshots output
csv_time_step_interval = 100 #time step interval for csv output
checkpoint_time_step_interval = 40 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##------------------------------------------------------------------------##

[Mesh]
  [msh]
    type = FileMeshGenerator
      file = '../../mesh/tpv26_200m_unifrom.msh'
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y > 0 & y < ${ymax_fault}'
    block_id = 100
  []
  [./new_block_2]
    type = ParsedSubdomainMeshGenerator
    input = new_block_1
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y < 0 & y > ${ymin_fault}'
    block_id = 200
  []

  [./split_1]
    type = BreakMeshByBlockGenerator
    input = new_block_2
    split_interface = true
    block_pairs = '100 200'
    show_info = true
  []

  [sidesets]
    type = SideSetsFromNormalsGenerator
    input = split_1
    normals = '-1  0  0
                1  0  0
                0 -1  0
                0  1  0
                0  0 -1
                0  0  1'
    fixed_normal = true
    new_boundary = 'left right bottom top back front'
  []

  [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = ${bottom_nodes_coord}
      new_boundary = corner_ptr
      input = sidesets
  []
[]

[GlobalParams]

  ##------------slip weakening------------##
  displacements = 'disp_x disp_y disp_z'
  porepressure = 'porepressure'

  #damping ratio
  q = ${q}

  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = ${lambda_o}

  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = ${shear_modulus_o}

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = ${xi_0}

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = ${xi_d}

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1000000000000

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1000000000000

  #if option 2, use Cd_constant
  Cd_constant = ${Cd_constant}

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = ${CdCb_multiplier}

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  # CBCBH_multiplier = 0.0
  CBH_constant = ${CBH_constant}

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = ${C_1}

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = ${C_2}

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = ${beta_width}

  #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  C_g = ${C_g}

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  m1 = ${m1}

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  m2 = ${m2}

  # energy ratio
  chi = ${chi}

  # Water bulk modulus (2.2 GPa)
  fluid_bulk_modulus = ${fluid_bulk_modulus}

  # Initial permeability (1 milli-darcy)
  permeability_solid_o = ${permeability_solid_o}

  # Initial porosity (15%)
  porosity_solid_o = ${porosity_solid_o}

  # Granular bulk modulus (15 GPa)
  solid_bulk_modulus_g = ${solid_bulk_modulus_g}

  # Solid grains bulk modulus (36 GPa - typical for quartz)
  solid_bulk_modulus_s = ${solid_bulk_modulus_s}

  permeability_evolution_with_damage = ${permeability_evolution_with_damage}
  initial_grain_size = ${initial_grain_size}
  ultimate_grain_size = ${ultimate_grain_size}
  initial_viscosity_fluid = ${initial_viscosity_fluid}

  anand_param_go_mat = ${anand_param_go_mat}
  anand_param_eta_cv_mat = ${anand_param_eta_cv_mat}
  anand_param_p_mat = ${anand_param_p_mat}
[]


[Variables]
    [disp_x]
        order = FIRST
        family = LAGRANGE
    []
    [disp_y]
        order = FIRST
        family = LAGRANGE
    []
    [disp_z]
        order = FIRST
        family = LAGRANGE
    []
    [porepressure]
        order = FIRST
        family = LAGRANGE
        scaling = 1E9
    []
[]


[AuxVariables]
  ###
  #slip weakening friction parameters
  [./resid_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./resid_y]
      order = FIRST
      family = LAGRANGE
  []
  [./resid_z]
    order = FIRST
    family = LAGRANGE
  []
  [./resid_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  [../]
  [./resid_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  [../]
  [./resid_slipweakening_z]
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
  [./resid_pressure_z]
      order = FIRST
      family = LAGRANGE
  [../]
  [./disp_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  []
  [./disp_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./disp_slipweakening_z]
    order = FIRST
    family = LAGRANGE
  []
  [./vel_slipweakening_x]
    order = FIRST
    family = LAGRANGE
  []
  [./vel_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./vel_slipweakening_z]
    order = FIRST
    family = LAGRANGE
  []
  ###
  #output initial shear stress
  [cohesion_aux]
    order = FIRST
    family = LAGRANGE
  []
  [forced_rupture_aux]
    order = FIRST
    family = LAGRANGE
  []
  [fluid_pressure_aux]
    order = FIRST
    family = LAGRANGE
  []
  ###
  #output jump, jump rate, traction quantities
  [jump_strike_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_strike_rate_aux]
    order = FIRST
    family = MONOMIAL
  []
  [traction_strike_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_normal_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_normal_rate_aux]
    order = FIRST
    family = MONOMIAL
  []
  [traction_normal_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_dip_aux]
    order = FIRST
    family = MONOMIAL
  []
  [jump_dip_rate_aux]
    order = FIRST
    family = MONOMIAL
  []
  [traction_dip_aux]
    order = FIRST
    family = MONOMIAL
  []
  ###
  #output CDB model properties
  [alpha_damagedvar_aux]
      order = FIRST
      family = MONOMIAL
  []
  [B_aux]
      order = FIRST
      family = MONOMIAL
  []
  [xi_aux]
      order = FIRST
      family = MONOMIAL
  []
  ###
  [deviatoric_strain_rate_aux]
    order = FIRST
    family = MONOMIAL
  []
  ##
  [eqstrain_nonlocal_aux]
    order = FIRST
    family = MONOMIAL
  []
  ##
  [eqstrain_nonlocal_initial_aux]
    order = FIRST
    family = MONOMIAL
  []
  ##
  [stress_xx]
    order = FIRST
    family = MONOMIAL
  []
  [stress_yy]
    order = FIRST
    family = MONOMIAL
  []
  [stress_xy]
    order = FIRST
    family = MONOMIAL
  []
  # sts_total components (from material property)
  [sts_total_xx]
    order = FIRST
    family = MONOMIAL
  []
  [sts_total_yy]
    order = FIRST
    family = MONOMIAL
  []
  [sts_total_zz]
    order = FIRST
    family = MONOMIAL
  []
  [sts_total_xy]
    order = FIRST
    family = MONOMIAL
  []
  [sts_total_xz]
    order = FIRST
    family = MONOMIAL
  []
  [sts_total_yz]
    order = FIRST
    family = MONOMIAL
  []
  [I1_aux]
    order = FIRST
    family = MONOMIAL
  []

[]

[Physics/SolidMechanics/CohesiveZone]
  [./czm_ik]
    boundary = 'Block100_Block200'
    strain = SMALL
    generate_output='traction_x traction_y traction_z jump_x jump_y jump_z'
  [../]
[]

[Problem]
  extra_tag_vectors = 'restore_tag'
[]

[AuxKernels]
  [Displacment_x]
    type = ProjectionAux
    variable = disp_slipweakening_x
    v = disp_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_y]
    type = ProjectionAux
    variable = disp_slipweakening_y
    v = disp_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_z]
    type = ProjectionAux
    variable = disp_slipweakening_z
    v = disp_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_x]
    type = CompVarRate
    variable = vel_slipweakening_x
    coupled = disp_x
    execute_on = 'TIMESTEP_END'
  []
  [Vel_y]
    type = CompVarRate
    variable = vel_slipweakening_y
    coupled = disp_y
    execute_on = 'TIMESTEP_END'
  []
  [Vel_z]
    type = CompVarRate
    variable = vel_slipweakening_z
    coupled = disp_z
    execute_on = 'TIMESTEP_END'
  []
  [Residual_x]
    type = ProjectionAux
    variable = resid_slipweakening_x
    v = resid_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_y]
    type = ProjectionAux
    variable = resid_slipweakening_y
    v = resid_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_z]
    type = ProjectionAux
    variable = resid_slipweakening_z
    v = resid_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [restore_x]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_x'
    variable = 'resid_x'
  []
  [restore_y]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_y'
    variable = 'resid_y'
  []
  [restore_z]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_z'
    variable = 'resid_z'
  []
  ### slip weakening cohesion
  [get_cohesion_aux]
    type = FunctionAux
    variable = cohesion_aux
    function = func_cohesion
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  ### slip weakening forced rupture
  [get_forced_rupture_aux]
    type = FunctionAux
    variable = forced_rupture_aux
    function = func_forced_rupture
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  ### slip weakening strike direction
  [get_jump_strike_aux]
    type = MaterialRealAux
    property = displacement_jump_strike
    variable = jump_strike_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  [get_jump_strike_rate_aux]
    type = MaterialRealAux
    property = displacement_jump_rate_strike
    variable = jump_strike_rate_aux
    execute_on = 'TIMESTEP_END'
    boundary = 'Block100_Block200'
  []
  [get_traction_strike_aux]
    type = MaterialRealAux
    property = traction_strike
    variable = traction_strike_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  ### slip weakening normal direction
  [get_jump_normal_aux]
    type = MaterialRealAux
    property = displacement_jump_normal
    variable = jump_normal_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  [get_jump_normal_rate_aux]
    type = MaterialRealAux
    property = displacement_jump_rate_normal
    variable = jump_normal_rate_aux
    execute_on = 'TIMESTEP_END'
    boundary = 'Block100_Block200'
  []
  [get_traction_normal_aux]
    type = MaterialRealAux
    property = traction_normal
    variable = traction_normal_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  ### slip weakening dip direction
  [get_jump_dip_aux]
    type = MaterialRealAux
    property = displacement_jump_dip
    variable = jump_dip_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  [get_jump_dip_rate_aux]
    type = MaterialRealAux
    property = displacement_jump_rate_dip
    variable = jump_dip_rate_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  [get_traction_dip_aux]
    type = MaterialRealAux
    property = traction_dip
    variable = traction_dip_aux
    boundary = 'Block100_Block200'
    execute_on = 'TIMESTEP_END'
  []
  ### get CDB model properties
  [get_alpha_damagedvar]
    type = MaterialRealAux
    variable = alpha_damagedvar_aux
    property = alpha_damagedvar
    execute_on = 'TIMESTEP_END'
  []
  [get_B]
    type = MaterialRealAux
    variable = B_aux
    property = B
    execute_on = 'TIMESTEP_END'
  []
  [get_xi]
    type = MaterialRealAux
    variable = xi_aux
    property = xi
    execute_on = 'TIMESTEP_END'
  []
  [get_I1]
    type = MaterialRealAux
    variable = I1_aux
    property = I1
    execute_on = 'TIMESTEP_END'
  []
  ###
  [get_deviatoric_strain_rate]
    type = MaterialRealAux
    variable = deviatoric_strain_rate_aux
    property = deviatoric_strain_rate
    execute_on = 'TIMESTEP_END'
  []
  ###
  [get_eqstrain_nonlocal_aux]
    type = MaterialRealAux
    variable = eqstrain_nonlocal_aux
    property = eqstrain_nonlocal
    execute_on = 'TIMESTEP_END'
   block = ${nonlocal_eqstrain_blocks}
  []
  ###
  [get_eqstrain_nonlocal_initial]
    type = SolutionAux
    variable = eqstrain_nonlocal_initial_aux
    solution = init_sol_components
    from_variable = xi
    execute_on = 'INITIAL'
  []
  ###
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  []

  # Extract sts_total components (material property)
  [sts_total_xx]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_xx
    index_i = 0
    index_j = 0
    execute_on = 'TIMESTEP_END'
  []
  [sts_total_yy]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_yy
    index_i = 1
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [sts_total_zz]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_zz
    index_i = 2
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []
  [sts_total_xy]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_xy
    index_i = 0
    index_j = 1
    execute_on = 'TIMESTEP_END'
  []
  [sts_total_xz]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_xz
    index_i = 0
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []
  [sts_total_yz]
    type = RankTwoAux
    rank_two_tensor = sts_total
    variable = sts_total_yz
    index_i = 1
    index_j = 2
    execute_on = 'TIMESTEP_END'
  []


[]

[Kernels]
  [dispkernel_x]
      type = TotalStressDivergenceTensor
      displacements = 'disp_x disp_y disp_z'
      variable = disp_x
      component = 0
      extra_vector_tags = 'restore_tag'
  []
  [dispkernel_y]
      type = TotalStressDivergenceTensor
      displacements = 'disp_x disp_y disp_z'
      variable = disp_y
      component = 1
      extra_vector_tags = 'restore_tag'
  []
  [dispkernel_z]
      type = TotalStressDivergenceTensor
      displacements = 'disp_x disp_y disp_z'
      variable = disp_z
      component = 2
      extra_vector_tags = 'restore_tag'
  []
  [./mass1]
      type = SmallStrainFluidSolidCouplingExplicit
      variable = porepressure
  [../]
  [./mass2]
      type = SmallStrainPorePressureTimeDerivativeExplicit
      variable = porepressure
  [../]
  [./darcy_flow]
      type = SmallStrainFluidDiffusion
      variable = porepressure
      fluid_density = 1000
      gravity_vector = '0 0 0'
  []
  [./darcy_flow_granular]
      type = SmallStrainFluidDiffusionGranular
      variable = porepressure
  []
  [./plastic_volumetric]
      type = SmallStrainPlasticVolumetricStrainCoupling
      variable = porepressure
  []
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_y
  []
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_z
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
  [./Reactionz]
    type = StiffPropDamping
    variable = 'disp_z'
    component = '2'
  []
  [pspg_stabilization]
    type = PoroStabilization
    variable = porepressure
    displacements = 'disp_x disp_y disp_z'
  []

  
[]

[Materials]
  [strain3]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
  []
  [stabilization_param]
    type = StabilizationMaterial
    stabilization_coeff = 1
  []
  #damage breakage model
  [stress_medium]
      type = ComputePoroDamageBreakageStress3DSlipWeakeningnonlocal
      output_properties = 'B alpha_damagedvar xi I1 I2 deviatoric_strain_rate'
      use_strain_rate_dependent_Cd = ${use_strain_rate_dependent_Cd}
      m_exponent = ${m_exponent}
      strain_rate_hat = ${strain_rate_hat}
      cd_hat = ${cd_hat}
      use_dilatancy = false
      use_nonlocal_eqstrain = true
      nonlocal_eqstrain_blocks = ${nonlocal_eqstrain_blocks}
      outputs = exodus
  []
  [dummy_material]
      type = GenericConstantMaterial
      prop_names = 'initial_damage initial_breakage damage_perturbation density'
      prop_values = '0 0 0 ${density}'
  []
  [eqstrain_nonlocal_initial_xi]
      type = FarmsCoupledVariableValueMaterial
      coupled_variable = eqstrain_nonlocal_initial_aux
      prop_name = eqstrain_nonlocal_initial
      output_properties = 'eqstrain_nonlocal_initial'
      outputs = exodus
  []
  [./czm_mat]
      type = PoroSlipWeakeningFrictionczm3dCDBM
      mu_s = ${mu_s}
      mu_d = ${mu_d}
      Dc = ${Dc}
      len = ${elem_size}
      disp_slipweakening_x     = disp_slipweakening_x
      disp_slipweakening_y     = disp_slipweakening_y
      disp_slipweakening_z     = disp_slipweakening_z
      vel_slipweakening_x      = vel_slipweakening_x
      vel_slipweakening_y      = vel_slipweakening_y
      vel_slipweakening_z      = vel_slipweakening_z
      reaction_slipweakening_x = resid_slipweakening_x
      reaction_slipweakening_y = resid_slipweakening_y
      reaction_slipweakening_z = resid_slipweakening_z
      reaction_pressure_x = resid_pressure_x
      reaction_pressure_y = resid_pressure_y
      reaction_pressure_z = resid_pressure_z
      fault_pressure = porepressure
      #---------------------------------------------#
      use_forced_rupture = true
      t0 = ${t0}
      cohesion_aux = cohesion_aux
      forced_rupture_aux = forced_rupture_aux
      #---------------------------------------------#
      boundary = 'Block100_Block200'
  [../]
  [./static_initial_strain_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_strain_tensor
      tensor_functions = 'func_initial_strain_xx   func_initial_strain_xy      func_initial_strain_xz
                          func_initial_strain_xy   func_initial_strain_yy      func_initial_strain_yz
                          func_initial_strain_xz   func_initial_strain_yz      func_initial_strain_zz'
      output_properties = 'static_initial_strain_tensor'
      outputs = exodus
  [../]
  [./static_initial_stress_tensor] #this is used in the SlipWeakeningFrictionczm3dCDBM
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz
                          func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                          func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
      output_properties = 'static_initial_stress_tensor'
      outputs = exodus
  [../]
  # Define initial pore pressure from function
  [./initial_porepressure]
    type = GenericFunctionMaterial
    prop_names = 'initial_porepressure'
    prop_values = 'func_initial_porepressure'
    output_properties = 'initial_porepressure'
    outputs = exodus
  [../]
  #nonlocal eqstrain #set initial value to be eqstrain_nonlocal_initial for the first step
  #the ComputeDamageBreakageStress3DSlipWeakeningNonlocal takes old value for updating damage/breakage
  [nonlocal_eqstrain_block1]
    type = ElkNonlocalEqstrainUpdated
    average_UO = eqstrain_averaging_block1
    block = ${nonlocal_eqstrain_blocks_1}
  []
  [nonlocal_eqstrain_block2]
    type = ElkNonlocalEqstrainUpdated
    average_UO = eqstrain_averaging_block2
    block = ${nonlocal_eqstrain_blocks_2}
  []
  [nonlocal_eqstrain_block3]
    type = ElkNonlocalEqstrainUpdated
    average_UO = eqstrain_averaging_block3
    block = ${nonlocal_eqstrain_blocks_3}
  []
  #for the block outside the region, nonlocal strain is equal to the local strain
  #[nonlocal_eqstrain_block]
  #  type = ParsedMaterial
  #  property_name = eqstrain_nonlocal
  #  coupled_variables = 'xi_aux'
  #  expression = 'xi_aux'
  #  block = ${local_eqstrain_blocks}
  #[]
[]

[Functions]
  ###strain field###
  [./func_initial_strain_xx]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_00'
  []
  [./func_initial_strain_xy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_01'
  []
  [./func_initial_strain_xz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_02'
  []
  [./func_initial_strain_yy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_11'
  []
  [./func_initial_strain_yz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_12'
  []
  [./func_initial_strain_zz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_22'
  []
  ###stress field###
  [./func_initial_stress_xx]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_00'
  []
  [./func_initial_stress_xy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_01'
  []
  [./func_initial_stress_xz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_02'
  []
  [./func_initial_stress_yy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_11'
  []
  [./func_initial_stress_yz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_12'
  []
  [./func_initial_stress_zz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_22'
  []
  [./func_initial_porepressure]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'porepressure'
  []
  ###cohesion###
  [./func_cohesion]
    type = InitialCohesionCDBMv2
    depth = ${cohesion_depth}
    slope = ${cohesion_slope}
    min_cohesion = ${cohesion_min}
  []
  ###forcedrupture###
  [./func_forced_rupture]
    type = ForcedRuptureTimeCDBMv2
    loc_x = ${nucl_center_x}
    loc_y = ${nucl_center_y}
    loc_z = ${nucl_center_z}
    r_crit = ${r_crit}
    Vs = ${Vs}
  []
[]

[UserObjects]
  [recompute_residual_tag]
      type = ResidualEvaluationUserObject
      vector_tag = 'restore_tag'
      force_preaux = true
      execute_on = 'TIMESTEP_END'
  []
  [./init_sol_components]
    type = SolutionUserObject
    mesh = '../static_solve/static_solve_out.e'
    system_variables = 'elastic_strain_00 elastic_strain_01 elastic_strain_02
                        elastic_strain_11 elastic_strain_12 elastic_strain_22
                        stress_00 stress_01 stress_02 stress_11 stress_12 stress_22 xi porepressure'
    timestep = LATEST
    force_preaux = true
    execute_on = 'INITIAL'
  [../]
  #here we avoid the cross-fault averaging by defining separate averaging blocks
  [eqstrain_averaging_block1]
    type = ElkRadialAverageUpdated
    length_scale = ${nonlocal_averaging_length_scale}
    prop_name = xi
    radius = ${nonlocal_averaging_radius}
    weights = BAZANT3D
    execute_on = TIMESTEP_END
    block = ${nonlocal_eqstrain_blocks_1}
  []
  [eqstrain_averaging_block2]
    type = ElkRadialAverageUpdated
    length_scale = ${nonlocal_averaging_length_scale}
    prop_name = xi
    radius = ${nonlocal_averaging_radius}
    weights = BAZANT3D
    execute_on = TIMESTEP_END
    block = ${nonlocal_eqstrain_blocks_2}
  []
  [eqstrain_averaging_block3]
    type = ElkRadialAverageUpdated
    length_scale = ${nonlocal_averaging_length_scale}
    prop_name = xi
    radius = ${nonlocal_averaging_radius}
    weights = BAZANT3D
    execute_on = TIMESTEP_END
    block = ${nonlocal_eqstrain_blocks_3}
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
  dt = ${dt}
  end_time = ${end_time}
  # num_steps = ${num_steps}
  [TimeIntegrator]
    type = CentralDifference
  []
[]

[Outputs]
  [exodus]
    type = Exodus
    execute_on = 'timestep_end'
    show = 'sts_total_xx sts_total_xy sts_total_zz sts_total_yy
            traction_strike_aux traction_normal_aux traction_dip_aux
            vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z
            stress_xx stress_yy stress_xy
            porepressure'
    time_step_interval = ${exodus_time_step_interval}
  []
  [csv]
    type = CSV
    execute_on = 'timestep_end'
    time_step_interval = ${csv_time_step_interval}
  []
  [out]
    type = Checkpoint
    time_step_interval = ${checkpoint_time_step_interval}
    num_files = ${checkpoint_num_files}
  []
  [sample_snapshots]
    type = Exodus
    execute_on = 'timestep_end'
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z alpha_damagedvar_aux B_aux xi_aux stress_xx stress_yy stress_xy deviatoric_strain_rate_aux eqstrain_nonlocal_aux porepressure'
    time_step_interval = ${sample_snapshots_time_step_interval}
  []
[]

[VectorPostprocessors]
  [main_fault]
    type = SideValueSampler
    variable = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z jump_strike_aux jump_dip_aux jump_normal_aux jump_strike_rate_aux jump_dip_rate_aux jump_normal_rate_aux traction_strike_aux traction_dip_aux traction_normal_aux alpha_damagedvar_aux B_aux xi_aux porepressure'
    boundary = 'Block100_Block200'
    sort_by = x
  []
  [off_fault]
    type = PositionsFunctorValueSampler
    functors = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z porepressure'
    positions = 'pos'
    sort_by = x
    execute_on = TIMESTEP_END
    discontinuous = false
  []
[]

#should use negative y as in tensile side
[Positions]
  [pos]
    type = InputPositions
    positions = '-24000 -1000 0
                 -20000 -1000 0
                 -16000 -1000 0
                 -12000 -1000 0
                 -8000 -1000 0
                 -4000 -1000 0
                 0 -1000 0
                 4000 -1000 0
                 8000 -1000 0
                 12000 -1000 0
                 16000 -1000 0
                 20000 -1000 0
                 24000 -1000 0
                 -24000 -2000 0
                 -20000 -2000 0
                 -16000 -2000 0
                 -12000 -2000 0
                 -8000 -2000 0
                 -4000 -2000 0
                 0 -2000 0
                 4000 -2000 0
                 8000 -2000 0
                 12000 -2000 0
                 16000 -2000 0
                 20000 -2000 0
                 24000 -2000 0
                 -24000 -3000 0
                 -20000 -3000 0
                 -16000 -3000 0
                 -12000 -3000 0
                 -8000 -3000 0
                 -4000 -3000 0
                 0 -3000 0
                 4000 -3000 0
                 8000 -3000 0
                 12000 -3000 0
                 16000 -3000 0
                 20000 -3000 0
                 24000 -3000 0
                 -24000 -4000 0
                 -20000 -4000 0
                 -16000 -4000 0
                 -12000 -4000 0
                 -8000 -4000 0
                 -4000 -4000 0
                 0 -4000 0
                 4000 -4000 0
                 8000 -4000 0
                 12000 -4000 0
                 16000 -4000 0
                 20000 -4000 0
                 24000 -4000 0'
  []
[]
