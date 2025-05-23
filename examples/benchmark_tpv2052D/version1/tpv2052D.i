# Method of Manufactured Solutions for 2D Dynamic Rupture with Slip-Weakening
# Based on the MATLAB script provided
[Mesh]
 [./msh]
 type = GeneratedMeshGenerator
 dim = 2
 nx = 51
 ny = 51
 xmin = -4000
 xmax = 4000
 ymin = -4000
 ymax = 4000
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
[]
[]
[GlobalParams]
 displacements = 'disp_x disp_y'
 q = 0.1
 Dc = 0.75
 T2_o = 50e6
 mu_d = 0.7
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
[]
[AuxVariables]
 [./vel_x]
 order = FIRST
 family = LAGRANGE
[]
 [./accel_x]
[]
 [./vel_y]
[]
 [./accel_y]
[]
 [./nodal_area]
 order = FIRST
 family = LAGRANGE
 [../]
 [./mms_exact_disp_x]
 order = FIRST
 family = LAGRANGE
[]
[]
[Functions]
 [./mms_displacement_x]
 type = ParsedFunction
 expression = '(delta/2) * ((1 + erf((t - tbar)/(sqrt(2)*tw))) + (Vmin * t)) * R_m / sqrt(x*x + y*y + R_m*R_m)'
 vars = 'delta R_m tbar tw Vmin'
 vals = '1.0 1000.0 0.25 0.1 1e-9'
[]
 [./initial_velocity_x]
 type = ParsedFunction
 expression = '(delta/(sqrt(2*pi)*tw)) * exp(-(0-tbar)^2/(2*tw^2)) * R_m / sqrt(x*x + y*y + R_m*R_m) + Vmin'
 vars = 'delta R_m tbar tw Vmin pi'
 vals = '1.0 1000.0 0.25 0.1 1e-9 3.14159265359'
[]
[]
[Modules/TensorMechanics/CohesiveZoneMaster]
 [./czm_ik]
 boundary = 'Block0_Block1'
 strain = SMALL
 generate_output='traction_x traction_y jump_x jump_y normal_traction tangent_traction normal_jump tangent_jump'
 [../]
[]
[Modules]
 [./TensorMechanics]
 [./Master]
 [./all]
 strain = SMALL
 add_variables = true
 planar_formulation = PLANE_STRAIN
 generate_output = 'stress_xx stress_yy stress_xy'
 [../]
 [../]
 [../]
[]
[AuxKernels]
 [velocity_x]
 type = CompVarRate
 variable = vel_x
 coupled = disp_x
[]
 [./mms_solution]
 type = FunctionAux
 function = mms_displacement_x
 variable = mms_exact_disp_x
 execute_on = 'INITIAL TIMESTEP_END'
[]
[]
[Kernels]
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
 [./mms_source_x]
 type = MMSSourceTermX
 variable = disp_x
 delta = 1.0
 R_m = 1000.0
 tbar = 0.25
 tw = 0.1
 Vmin = 1e-9
[]
 [./mms_source_y]
 type = MMSSourceTermY
 variable = disp_y
 delta = 1.0
 R_m = 1000.0
 tbar = 0.25
 tw = 0.1
 Vmin = 1e-9
[]
[]
[ICs]
 [./disp_x_ic]
 type = FunctionIC
 variable = disp_x
 function = mms_displacement_x
[]
 [./vel_x_ic]
 type = FunctionIC
 variable = vel_x
 function = initial_velocity_x
[]
[]
[Materials]
 [elasticity]
 type = ComputeIsotropicElasticityTensor
 lambda = 20e9
 shear_modulus = 30e9
 use_displaced_mesh = false
[]
 [stress]
 type = ComputeLinearElasticStress
[]
 [density]
 type = GenericConstantMaterial
 prop_names = density
 prop_values = 2850
[]
 [./czm_stress_derivative]
 type = StressDerivative2
 boundary = 'Block0_Block1'
 args = 'disp_x disp_y'
 [../]
 [./czm_mat]
 type = SlipWeakeningFriction2dxx
 boundary = 'Block0_Block1'
 nodal_area = nodal_area
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
[Postprocessors]
 [./l2_error]
 type = ElementL2Error
 function = mms_displacement_x
 variable = disp_x
[]
 [./max_vel_x]
 type = ElementExtremeValue
 variable = vel_x
[]
[]
[Executioner]
 type = Transient
 dt = 0.005
 end_time = 0.75
 automatic_scaling = true
 [TimeIntegrator]
 type = CentralDifference
[]
[]
[Outputs]
 exodus = true
 interval = 5
 [./console]
 type = Console
 output_on = 'timestep_end'
 execute_postprocessors_on = 'TIMESTEP_END'
[]
[]