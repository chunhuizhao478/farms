#parameters

##mesh parameters
###need to check gmsh file for changing the parameters here
bottom_nodes_coord =' -60000 -60000 -60000;
                      60000 -60000 -60000;
                      60000 60000  -60000;
                     -60000 60000  -60000'

xmin_fault = -20000 #xmin of fault
xmax_fault = 20000 #xmax of fault
zmin_fault = -20000 #zmin of fault
# zmax_fault = 0 #zmax of fault
elem_size = 200 #!!! element size near the fault, need to be consistent with the mesh file
##-------------------------##

##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
# Cs = '${fparse shear_modulus_o / density }' #shear wave speed
# Cp = '${fparse (lambda_o + 2 * shear_modulus_o) / density }' #pressure wave speed
##-------------------------##

##Slip weakening parameters##
Dc = 0.8 #0.4 #characteristic length (m)
q = 0.4 #damping ratio
mu_s = 0.8 #static friction coefficient
mu_d = 0.6 #dynamic friction coefficient
##-------------------------##

##Cohesion parameters##
cohesion_depth = 5000 #cohesion depth (m)
cohesion_slope = 0.00072 #cohesion slope (MPa/m)
cohesion_min = 0.4 #minimum cohesion value (MPa)
##---------------------------------------------##

##CDB model parameters##
xi_0 = -0.8 #strain invariants ratio: onset of damage evolution
xi_d = -0.8 #strain invariants ratio: onset of breakage healing

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
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

##initial stress parameters##

#background stress
fluid_density = 1000
gravity = 9.8
bxx = 0.926793
byy = 1.073206
bxy = -0.8
##------------------------------------------------------------------##

##tapering parameters##
use_tapering = true #use tapering to reduce deviatoric stress components at shallow depth
tapering_depth_A = 15000 #depth at which tapering starts to be applied (m)
tapering_depth_B = 20000 #depth at which tapering stops to be applied (m)
##------------------------------------------------------------------##

#nucleation parameters
nucl_center_x = -16000 #nucleation center x coordinate
nucl_center_y = 0 #nucleation center y coordinate
nucl_center_z = -10000 #nucleation center y coordinate
r_crit = 4000 #critical distance to hypocenter (m)
Vs = 3464 #shear wave speed (m/s)
t0 = 0.5 #nucleation time (s)
##------------------------------------------------------------------##

##initial damage parameters
sigma = 5e2
peak_val = 0
len_of_fault_strike = 40000
len_of_fault_dip = 20000
fault_center = '0 0 -10000'
##-------------------------##

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../../mesh/tpv26_100m_nonlocal_occ_40kmfault_uniform200m.msh'
    #file = '../../mesh/tpv26_400m_nonlocal_occ.msh'
  []
  [./sidesets]
    input = msh
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0
                0 0 -1
                0 0 1'
    new_boundary = 'left right back front bottom top'
  []
  [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = ${bottom_nodes_coord}
      new_boundary = corner_ptr
      input = sidesets
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

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
  xi_max = 1.8

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

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
[]

[AuxVariables]
[]

[AuxKernels]
[]

[Kernels]
  [SolidMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
  [gravity_z]
    type = EffectiveBodyForceTPV26VaryingDensity
    variable = disp_z
    fluid_density = ${fluid_density}
    gravity = ${gravity}
    use_overpressure = false
  []
[]

[Materials]
  [./depth_dependent_props]
    type = TPV32DepthSeismicProperties
    depth_axis = z           # z is vertical in this model
    flip_sign = true         # domain uses negative z for depth; flip to get positive depth
    clamp_nonpositive_depth = true
    vs_property_name = "Vs"
    vp_property_name = "Vp"
    density_property_name = "density"
    shear_modulus_property_name = "shear_modulus_input"
    lambda_property_name = "lambda_input"
    output_properties = 'lambda_input shear_modulus_input density'
    outputs = exodus
  [../]
  [./elasticity_tensor]
    type = FarmsComputeIsotropicElasticityTensor
    shear_modulus = ${shear_modulus_o}   # fallback value; overridden by *_input below
    lambda = ${lambda_o}                 # fallback value; overridden by *_input below
    shear_modulus_input = shear_modulus_input #material property from depth_dependent_props
    lambda_input = lambda_input #material property from depth_dependent_props
    xi_o = ${xi_0}
    initial_damage_property = "initial_damage"
  [../]
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = ini_stress
    outputs = exodus
  []
  [stress_medium]
    type = ComputeDamageBreakageStress3DStatic
    shear_modulus_input = shear_modulus_input #material property from depth_dependent_props
    lambda_input = lambda_input #material property from depth_dependent_props
    output_properties = 'B alpha_damagedvar xi I1 I2 stress elastic_strain'
    outputs = exodus
  []
  [initial_damage_surround]
    type = InitialDamageCycleSim3DPlane
    sigma = ${sigma}
    peak_val = ${peak_val}
    len_of_fault_strike = ${len_of_fault_strike}
    len_of_fault_dip = ${len_of_fault_dip}
    nucl_center = ${fault_center}
    output_properties = 'initial_damage'
    outputs = exodus
  []
  [dummy_material]
      type = GenericConstantMaterial
      prop_names = 'initial_breakage damage_perturbation'
      prop_values = '0 0'
  []
  [./static_initial_strain_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz
                      func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                      func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    eigenstrain_name = ini_stress
  [../]
  [./static_initial_stress_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening, SlipWeakeningFrictionczm3dCDBM
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz
                          func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                          func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
      output_properties = 'static_initial_stress_tensor'
      outputs = exodus
  [../]
  [./comp_xi]
    type = ComputeXi
    output_properties = 'strain_invariant_ratio'
    outputs = exodus
  []
[]

[BCs]
  [fix_bottom_z]
      type = DirichletBC
      variable = disp_z
      boundary = bottom
      value = 0
  []
  [static_pressure_left]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = left
      function = func_pos_xx_stress
      displacements = 'disp_x disp_y disp_z'
  []
  [static_pressure_right]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = right
      function = func_neg_xx_stress
      displacements = 'disp_x disp_y disp_z'
  []
  #
  [static_pressure_front]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = front
      function = func_neg_yy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  [static_pressure_back]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = back
      function = func_pos_yy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  #
  [static_pressure_front_shear]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = front
      function = func_pos_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  [static_pressure_back_shear]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = back
      function = func_neg_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  [static_pressure_left_shear]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = left
      function = func_neg_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  [static_pressure_right_shear]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = right
      function = func_pos_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []
  #
  [fix_node_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'corner_ptr'
    value = 0
  [../]
  [fix_node_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'corner_ptr'
    value = 0
  [../]
  [fix_node_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'corner_ptr'
    value = 0
  [../]
[]

[Functions]
  ###stress field###
  [./func_initial_stress_xx]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 1
    j = 1
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    # rock density is taken from TPV32 profile inside the function
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  [./func_pos_xx_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xx'
    scale_factor = '-1'
  [../]
  [./func_neg_xx_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xx'
    scale_factor = '1'
  [../]
  ##
  [./func_initial_stress_xy]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 1
    j = 2
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  [./func_pos_xy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xy'
    scale_factor = '1'
  [../]
  [./func_neg_xy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xy'
    scale_factor = '-1'
  [../]
  ##
  [./func_initial_stress_xz]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 1
    j = 3
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  ##
  [./func_initial_stress_yy]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 2
    j = 2
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  [./func_pos_yy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_yy'
    scale_factor = '-1'
  [../]
  [./func_neg_yy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_yy'
    scale_factor = '1'
  [../]
  ##
  [./func_initial_stress_yz]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 2
    j = 3
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  [./func_initial_stress_zz]
    type = InitialStressStrainTPV26VaryingDensity
    get_initial_stress = true
    i = 3
    j = 3
    lambda_o = ${lambda_o}
    shear_modulus_o = ${shear_modulus_o}
    fluid_density = ${fluid_density}
    flip_sign = true
    gravity = ${gravity}
    bxx = ${bxx}
    byy = ${byy}
    bxy = ${bxy}
    use_tapering = ${use_tapering}
    tapering_depth_A = ${tapering_depth_A}
    tapering_depth_B = ${tapering_depth_B}
  []
  [./func_pos_zz_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_zz'
    scale_factor = '1'
  [../]
  [./func_neg_zz_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_zz'
    scale_factor = '-1'
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  solve_type = NEWTON
  type = Steady

  nl_abs_tol = 1E-10
  nl_rel_tol = 1E-12
  l_tol = 1E-7
  l_max_its = 200
  nl_max_its = 400
  line_search  = 'basic'
  automatic_scaling = true
  verbose = true
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -ksp_type -ksp_gmres_restart'
  petsc_options_value = ' asm      2              hypre             gmres     200'
[]

[Outputs]
  exodus = true
[]
