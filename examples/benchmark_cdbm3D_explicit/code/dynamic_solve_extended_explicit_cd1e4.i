#continuum damage-breakage model dynamics

##########################################################################################################################################
#Mesh section
#FileMeshGenerator: read mesh file
#SideSetsFromNormalsGenerator: generate side sets from normals
#ExtraNodesetGenerator: generate extra nodeset - here we use it to define corner points associated with the bottom boundary
##########################################################################################################################################
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_large_extended.msh'
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
        new_boundary = 'left right front back bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = ' -40000 -40000 -40000;
                   40000 -40000 -40000;
                   40000 40000  -40000;
                  -40000 40000  -40000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.8
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 1e4

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    # energy ratio
    chi = 0.7

    #
    D = 0
    
[]

##############################################
#Variables section
#disp_x, disp_y, disp_z: displacement field
##############################################
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

############################################################################################
#AuxVariables section
#gradient of damage variable (not used): alpha_grad_x, alpha_grad_y, alpha_grad_z
#initial_shear_stress_aux: initial shear stress
############################################################################################
[AuxVariables]
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [vel_x]
    []  
    [vel_y]
    []
    [vel_z]
    []
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [alpha_damagedvar_aux]
        order = FIRST
        family = MONOMIAL
    []
    [xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    [B_aux]
        order = FIRST
        family = MONOMIAL
    []
[]

#############################################################################################################
#AuxKernels section
#CompVarRate: compute the rate of a variable: vel_x, vel_y, vel_z
#SolutionAux: get a solution from static solve and define in an auxiliary variable: initial_shear_stress_aux
#############################################################################################################
[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
        type = CompVarRate
        variable = vel_z
        coupled = disp_z
        execute_on = 'TIMESTEP_END'
    []
    [get_initial_damage]
        type = MaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
    []
    [get_damage]
        type = MaterialRealAux
        variable = alpha_damagedvar_aux
        property = alpha_damagedvar
        block = '1 3'
    []
    [get_strain_invariant_ratio]
        type = MaterialRealAux
        variable = xi_aux
        property = xi
        block = '1 3'
    []
    [get_B]
        type = MaterialRealAux
        variable = B_aux
        property = B
        block = '1 3'
    []
[]

##############################################
#Kernel section
#StressDivergenceTensors: compute the divergence of stress tensor, in all three directions
#InertialForce: compute the inertial force, in all three directions
##############################################
[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
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
    [gravity]
        type = Gravity
        variable = disp_z
        value = -9.81
    []       
[]

###################################################################
#Materials section
#ComputeSmallStrain: compute the small strain
#GenericConstantMaterial: define the density
#ComputeDamageBreakageStress3D: compute the stress field
#ComputeLinearElasticStress: compute the elastic stress field
#ComputeIsotropicElasticityTensor: compute the elasticity tensor
#InitialDamageCycleSim3DPlane: define the initial damage field
#InitialBreakageCycleSim3DPlane: define the initial breakage field
#PerturbationRadial: define the perturbation field
#ParsedMaterial: define the initial shear stress field
###################################################################
#Block 1, 3: use continuum damage breakage model
#Block 2: use linear elastic model
###################################################################
[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        # outputs = exodus
    [] 
    [density]
        type = GenericConstantMaterial
        prop_names = 'density nonADdensity'
        prop_values = '2700 2700'
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBM
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'eps_p eps_e I1 I2 stress xi strain_rate cd_ratedependent alpha_damagedvar B xi'
        block = '1 3'
        outputs = exodus
    [] 
    [stress_elastic]
        type = ComputeLinearElasticStress
        block = '2'
        output_properties = 'elastic_strain stress'
        outputs = exodus
    []
    [elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 30e9
        shear_modulus = 30e9
    []
    ################################################################################
    #initial damage field
    #sigma = 5e2: sigma value
    #peak_val = 0.7: peak value of the initial damage
    #len_of_fault_strike = 8000: length of the fault in the x-direction
    #len_of_fault_dip = 3000: length of the fault in the z-direction
    #nucl_center = '0 0 -7500': nucleation center
    ################################################################################
    #within the damage zone plane: len_of_fault_strike by len_of_fault_dip
    #the initial damage value = 0.7
    #outside the damage zone plane: initial damage experience exponential decay
    #Real alpha_o = _peak_val * std::exp(-1.0 * (r * r) / (_sigma * _sigma));
    #r = std::sqrt(dx * dx + dy * dy + dz * dz);
    ################################################################################
    [initial_damage_surround]
        type = InitialDamageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.7
        len_of_fault_strike = 14000
        len_of_fault_dip = 10000
        nucl_center = '0 0 -10000'
        use_damage_perturb = true
        output_properties = 'initial_damage'      
        outputs = exodus
    []
    ################################################################################
    #initial breakage field
    #sigma = 5e2: sigma value
    #peak_val = 0.1: peak value of the initial breakage
    #len_of_fault_strike = 8000: length of the fault in the x-direction
    #len_of_fault_dip = 3000: length of the fault in the z-direction
    #nucl_center = '0 0 -7500': nucleation center
    ################################################################################
    #within the breakage zone plane: len_of_fault_strike by len_of_fault_dip
    #the initial breakaege value = 0.7
    #outside the breakage zone plane: initial damage experience exponential decay
    #Real alpha_o = _peak_val * std::exp(-1.0 * (r * r) / (_sigma * _sigma));
    #r = std::sqrt(dx * dx + dy * dy + dz * dz);
    ################################################################################
    [initial_breakage_surround]
        type = InitialBreakageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.0
        len_of_fault_strike = 14000
        len_of_fault_dip = 10000
        nucl_center = '0 0 -10000'
        use_breakage_perturb = true
        output_properties = 'initial_breakage'      
        outputs = exodus
    []
    ################################################################################
    #perturbation field
    #perturbation field is defined as a radial shear stress perturbation
    #peak_value = 10e6: peak value of the perturbation
    #thickness = 200: thickness of the perturbation (along normal y direction)
    #length = 1000: length of the perturbation (along normal x,z direction)
    #duration = 1.0: duration of the perturbation
    #perturbation_type = 'shear_stress': perturbation type
    #sigma_divisor = 2.0: divisor of the sigma
    #see the code snippet below for more details
    # // Get the current point coordinates in the mesh
    # const Real xcoord = _q_point[_qp](0); // strike direction
    # const Real ycoord = _q_point[_qp](1); // normal direction
    # const Real zcoord = _q_point[_qp](2); // dip direction
    #
    # // We define a 2D Gaussian in the XZ plane, ignoring y in the exponent
    # // The characteristic "sigma" is length / sigma_divisor
    # const Real sigma_x = _length / _sigma_divisor;
    # const Real sigma_z = _length / _sigma_divisor;
    # const Real gaussian_factor = _peak_value; // maximum amplitude of the Gaussian
    #
    # // Distance in XZ from the nucleation center
    # // If your center is (0, 0, -7500), then _nucl_center might be [0, 0, -7500].
    # const Real dx = xcoord - _nucl_center[0];
    # const Real dz = zcoord - _nucl_center[2];
    #
    # // 2D Gaussian distribution in the XZ plane:
    # //    G(x,z) = peak_value * exp( - [dx^2 / (2*sigma_x^2) + dz^2 / (2*sigma_z^2)] )
    # const Real gaussian_value = gaussian_factor *
    #                            std::exp(-((dx * dx) / (2.0 * sigma_x * sigma_x) +
    #                                       (dz * dz) / (2.0 * sigma_z * sigma_z)));
    ################################################################################
    [shear_stress_perturbation]
        type = PerturbationRadial
        nucl_center = '0 0 -7500'
        peak_value = 0.3
        thickness = 200
        length = 1000
        duration = 0.01
        perturbation_type = 'damage'
        sigma_divisor = 2.0
        output_properties = 'breakage_perturbation shear_stress_perturbation damage_perturbation'
        outputs = exodus
    []
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'initial_shear_stress'
        prop_values = '0'
    []
[]  

[Functions]
[]

[Postprocessors]
    [./maxvelx]
        type = NodalExtremeValue
        variable = vel_x
    [../]
    [./maxvely]
        type = NodalExtremeValue
        variable = vel_y
    [../]
    [./maxvelz]
        type = NodalExtremeValue
        variable = vel_z
    [../]
[../]

#########################################################################################################
#Executioner Section
#type = Transient: Specify the type of executioner, in this case, a transient (time-dependent) simulation
#dt = 1e-4: Specify the time step size
#end_time = 10.0: Specify the end time of the simulation
#num_steps = 10: Optionally, you can specify the number of steps instead of end_time
##TimeIntegrator: Specify the time integrator to use
##type = CentralDifference: Specify the type of time integrator, in this case, a central difference scheme
##solve_type = consistent: Specify the solve type, in this case, a consistent solve (as opposed to lumped)
##use_constant_mass = true: Optionally, you can specify to use constant mass matrix
#CFL condition needs to be satisfied: dt < factor * dx / pressure_wave_speed
#########################################################################################################
[Executioner]
    type = Transient
    dt = 1e-3
    end_time = 100.0
    # num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = consistent
        # use_constant_mass = true
    []
[]

#########################################################################################################
#Outputs Section
#exodus = true: Specify that we want to output the solution to an Exodus file
#time_step_interval = 1: Optionally, you can specify the time step interval at which to output the solution
#show = : Optionally, you can specify which variables to output to the Exodus file
##[./csv]: Optionally, you can specify to output the solution to a CSV file
##type = CSV: Specify the type of output, in this case, a CSV file
##time_step_interval = 1: Optionally, you can specify the time step interval at which to output the solution
##show = : Optionally, you can specify which variables to output to the CSV file
#########################################################################################################
[Outputs] 
    ### save the solution to a exodus file every [time_step_interval] time steps]
    exodus = true
    time_step_interval = 100
    #############################################
    ##disp_x, disp_y, disp_z: displacement field
    ##vel_x, vel_y, vel_z: velocity field
    ##alpha_damagedvar: damage variable
    ##B: breakage variable
    ##initial_damage: initial damage field
    ##initial_breakage: initial breakage field
    ##stress_00, stress_01, stress_02, stress_11, stress_12, stress_22: stress field
    ##eps_e_00, eps_e_01, eps_e_02, eps_e_11, eps_e_12, eps_e_22: elastic strain field
    ##eps_p_00, eps_p_01, eps_p_02, eps_p_11, eps_p_12, eps_p_22: plastic strain field
    ##I1, I2: strain invariants
    ##xi: strain invariants ratio
    ##shear_stress_perturbation: perturbation field
    #############################################
    show = 'alpha_damagedvar B xi vel_x vel_y vel_z alpha_damagedvar_aux B_aux xi_aux'
    # [./csv]
    #     type = CSV
    #     time_step_interval = 1
    #     show = 'maxvelx maxvely maxvelz'
    # [../]
[]

#We assume the simulation is loaded with compressive pressure and shear stress
#############################################################################################################################################################
#BCs Section
#Note: use neuamnnBC gives minimum waves than pressureBC
#Coordinate system: x: strike, y: normal , z: dip
#The following BCs are applied to the top, bottom, left, right, front, and back boundaries
#Note depends on the coordinate system, the normal direction is different
#Confinement pressure: static_pressure_top, static_pressure_bottom, static_pressure_left, static_pressure_right, static_pressure_front, static_pressure_back
#Shear stress: static_pressure_front_shear, static_pressure_back_shear
#Constraints on corner_ptr: fix_cptr1_x, fix_cptr1_y, fix_cptr1_z
#############################################################################################################################################################
[BCs]
    # [static_pressure_top]
    #     type = NeumannBC
    #     variable = disp_z
    #     boundary = top
    #     value = 0
    #     displacements = 'disp_x disp_y disp_z'
    # []
    # [static_pressure_bottom]
    #     type = NeumannBC
    #     variable = disp_z
    #     boundary = bottom
    #     value = 5.2974e8
    #     displacements = 'disp_x disp_y disp_z'
    # []
    [fix_bottom_z]
        type = DirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []
    #Note: use neuamnnBC gives minimum waves than pressureBC  
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
        function = func_pos_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = FunctionNeumannBC
        variable = disp_y
        boundary = back
        function = func_neg_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []
    #
    [static_pressure_front_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = front
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = FunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_pos_xy_stress
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
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = DirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]

[Functions]
    [func_pos_yy_stress]
        type = ParsedFunction      
        expression = 'if(-z<15600, -1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_yy_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (1.073206 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'  
    []
    [func_pos_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600, -1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), -1 * (-2700 * 9.81 * (-z)))'
    []
    [func_neg_xx_stress]
        type = ParsedFunction
        expression = 'if(-z<15600,  1 * (0.926793 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) ) - (1000 * 9.81 * (-z))), 1 * (-2700 * 9.81 * (-z)))'
    []
    [func_pos_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, -1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, -1 * (-0.75 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
    [func_neg_xy_stress]
        type = ParsedFunction
        # expression = 'if(-z<15600, 1 * (-0.169029 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
        expression = 'if(-z<15600, 1 * (-0.75 * ( (-2700 * 9.81 * (-z)) + (1000 * 9.81 * (-z)) )), 0.0)'
    []
[]

#####################################
#BCs Section
#Absorbing boundary conditions
#####################################
# [BCs]
#     ##non-reflecting bc
#     #
#     [./dashpot_top_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = top
#     []
#     [./dashpot_top_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = top
#     []
#     [./dashpot_top_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = top
#     []
#     #
#     [./dashpot_bottom_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = bottom
#     []
#     [./dashpot_bottom_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = bottom
#     []
#     [./dashpot_bottom_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = bottom
#     []
#     #
#     [./dashpot_left_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = left
#     []
#     [./dashpot_left_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = left
#     []
#     [./dashpot_left_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = left
#     []
#     #
#     [./dashpot_right_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = right
#     []
#     [./dashpot_right_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = right
#     []
#     [./dashpot_right_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = right
#     []
#     #
#     [./dashpot_front_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = front
#     []
#     [./dashpot_front_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = front
#     []
#     [./dashpot_front_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = front
#     []
#     #
#     [./dashpot_back_x]
#         type = NonReflectDashpotBC3d
#         component = 0
#         variable = disp_x
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = back
#     []
#     [./dashpot_back_y]
#         type = NonReflectDashpotBC3d
#         component = 1
#         variable = disp_y
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = back
#     []
#     [./dashpot_back_z]
#         type = NonReflectDashpotBC3d
#         component = 2
#         variable = disp_z
#         disp_x = disp_x
#         disp_y = disp_y
#         disp_z = disp_z
#         p_wave_speed = 5773.5
#         shear_wave_speed = 3333.3
#         boundary = back
#     []
# [] 

#############################################################################################################################################################
#UserObjects Section
#init_sol_components: SolutionUserObject
##mesh: The mesh file to use, in this case, the static solution mesh
##system_variables: The system variables to use, in this case, disp_x, disp_y, disp_z
##timestep: The time step to use, in this case, the latest time step (LATEST), static solve has two time steps, initial and latest
##force_preaux: Whether to force the preaux to be used, in this case, true, then SolutioAux will be used 
#############################################################################################################################################################
[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_extended_explicit_out.e'
      system_variables = 'disp_x disp_y disp_z'
      timestep = LATEST
      force_preaux = true
    [../]
[]

#############################################################################################################################################################
#ICs Section
#Initial conditions are applied to the solution variables from the "init_sol_components" in [UserObjects] section
#for all three direction: disp_x_ic, disp_y_ic, disp_z_ic
#############################################################################################################################################################
[ICs]
    [disp_x_ic]
      type = SolutionIC
      variable = disp_x
      solution_uo = init_sol_components
      from_variable = disp_x
    []
    [disp_y_ic]
      type = SolutionIC
      variable = disp_y
      solution_uo = init_sol_components
      from_variable = disp_y
    []
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]
