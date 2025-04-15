#continuum damage-breakage model dynamics

##########################################################################################################################################
#User Parameters section
dt = 1e-3 #time step size
end_time = 20.0 #end time of the simulation
time_step_interval = 100 #output interval

##########################################################################################################################################
#Mesh section
#FileMeshGenerator: read mesh file
#SideSetsFromNormalsGenerator: generate side sets from normals
#ExtraNodesetGenerator: generate extra nodeset - here we use it to define corner points associated with the bottom boundary
##########################################################################################################################################
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh.msh'
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
        coord = ' -15000 -15000 -15000;
                   15000 -15000 -15000;
                   15000 15000  -15000;
                  -15000 15000  -15000'
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
    xi_d = -0.9
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 0

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
    C_g = 1e-8
    
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
#velocity field: vel_x, vel_y, vel_z
#acceleration field: accel_x, accel_y, accel_z
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
    [accel_x]
    []
    [accel_y]
    []
    [accel_z]
    []
    [initial_shear_stress_aux]
        order = CONSTANT
        family = MONOMIAL
    []
[]

#############################################################################################################
#AuxKernels section
#NewmarkVelAux: compute the rate of a variable: vel_x, vel_y, vel_z using Nemwark time integration
#NewmarkAccelAux: compute the acceleration of a variable: accel_x, accel_y, accel_z using Nemwark time integration
#SolutionAux: get a solution from static solve and define in an auxiliary variable: initial_shear_stress_aux
#############################################################################################################
[AuxKernels]
    [accel_x]
        type = NewmarkAccelAux
        variable = accel_x
        displacement = disp_x
        velocity = vel_x
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_x]
        type = NewmarkVelAux
        variable = vel_x
        acceleration = accel_x
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_y]
        type = NewmarkAccelAux
        variable = accel_y
        displacement = disp_y
        velocity = vel_y
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_y]
        type = NewmarkVelAux
        variable = vel_y
        acceleration = accel_y
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [accel_z]
        type = NewmarkAccelAux
        variable = accel_z
        displacement = disp_z
        velocity = vel_z
        beta = 0.25
        execute_on = 'TIMESTEP_END'
    []
    [vel_z]
        type = NewmarkVelAux
        variable = vel_z
        acceleration = accel_z
        gamma = 0.5
        execute_on = 'TIMESTEP_END'
    []
    [initial_shear_stress_aux]
        type = SolutionAux
        variable = initial_shear_stress_aux
        solution = init_sol_components
        from_variable = stress_01
    []
[]

##############################################
#Kernel section
#StressDivergenceTensors: compute the divergence of stress tensor, in all three directions
#InertialForce: compute the inertial force, in all three directions, using Newmark time integration
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
        acceleration = accel_x
        velocity = vel_x
        beta = 0.25
        gamma = 0.5
        eta = 0
    []
    [./inertia_y]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_y
        acceleration = accel_y
        velocity = vel_y
        beta = 0.25
        gamma = 0.5
        eta = 0
    [] 
    [./inertia_z]
        type = InertialForce
        use_displaced_mesh = false
        variable = disp_z
        acceleration = accel_z
        velocity = vel_z
        beta = 0.25
        gamma = 0.5
        eta = 0
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
        prop_names = 'density'
        prop_values = '2700'
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3DDynamicCDBM
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 xi stress'
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
    #nucl_center = '0 0 -7500': plane center
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
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
        output_properties = 'initial_damage'      
        outputs = exodus
    []
    ################################################################################
    #initial breakage field
    #sigma = 5e2: sigma value
    #peak_val = 0.1: peak value of the initial breakage
    #len_of_fault_strike = 8000: length of the fault in the x-direction
    #len_of_fault_dip = 3000: length of the fault in the z-direction
    #nucl_center = '0 0 -7500': plane center
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
        peak_val = 0.1
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
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
        peak_value = 10e6
        thickness = 200
        length = 1000
        duration = 1.0
        perturbation_type = 'shear_stress'
        sigma_divisor = 2.0
        output_properties = 'shear_stress_perturbation damage_perturbation'
        outputs = exodus
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
[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
    type = Transient
    solve_type = 'NEWTON'
    # solve_type = 'PJFNK'
    start_time = -1e-12
    end_time = ${end_time}
    # num_steps = 10
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-6
    nl_max_its = 8
    nl_abs_tol = 1e-8
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    dt = ${dt}
    [./TimeIntegrator]
        type = NewmarkBeta
        beta = 0.25
        gamma = 0.5
    [../]
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
    time_step_interval = ${time_step_interval}
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
    show = 'disp_x disp_y disp_z vel_x vel_y vel_z alpha_damagedvar B initial_damage initial_breakage stress_00 stress_01 stress_02 stress_11 stress_12 stress_22 eps_e_00 eps_e_01 eps_e_02 eps_e_11 eps_e_12 eps_e_22 eps_p_00 eps_p_01 eps_p_02 eps_p_11 eps_p_12 eps_p_22 xi shear_stress_perturbation'
    [./csv]
        type = CSV
        time_step_interval = ${time_step_interval}
        show = 'maxvelx maxvely maxvelz'
    [../]
[]

#############################################################################################################
#Controls Section
#This is used for first steady state solve, we close all inertial terms, absorbing boundary conditions
#############################################################################################################
[Controls] # turns off inertial terms for the FIRST time step
  [./period0]
    type = TimePeriod
    disable_objects = '*/vel_x */vel_y */vel_z */accel_x */accel_y */accel_z */inertia_x */inertia_y */inertia_z */dashpot_front_x */dashpot_front_y */dashpot_front_z */dashpot_back_x */dashpot_back_y */dashpot_back_z */dashpot_left_x */dashpot_left_y */dashpot_left_z */dashpot_right_x */dashpot_right_y */dashpot_right_z */dashpot_top_x */dashpot_top_y */dashpot_top_z */dashpot_bottom_x */dashpot_bottom_y */dashpot_bottom_z'
    start_time = -1e-12
    end_time = ${dt} # dt used in the simulation
  []
[../]

#parameters for the initial stress field
################################################
bxx = 0.926793
byy = 1.073206
bxy = -0.169029
linear_variation_cutoff_distance = 15600
################################################
[Functions]
    [func_pos_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_yy_stress]
        type = InitialDepthDependentStress
        i = 2
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xx_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 1
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_pos_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_xy_stress]
        type = InitialDepthDependentStress
        i = 1
        j = 2
        pos_sign = false
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
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
    #     type = ADNeumannBC
    #     variable = disp_z
    #     boundary = top
    #     value = 0
    #     displacements = 'disp_x disp_y disp_z'
    # []
    # [static_pressure_bottom]
    #     type = ADNeumannBC
    #     variable = disp_z
    #     boundary = bottom
    #     value = 5.2974e8
    #     displacements = 'disp_x disp_y disp_z'
    # []
    [fix_bottom_z]
        type = ADDirichletBC
        variable = disp_z
        boundary = bottom
        value = 0
    []
    #Note: use neuamnnBC gives minimum waves than pressureBC  
    [static_pressure_left]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = left
        function = func_pos_xx_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = right
        function = func_neg_xx_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    #
    [static_pressure_front]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = front
        function = func_pos_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = back
        function = func_neg_yy_stress
        displacements = 'disp_x disp_y disp_z'
    []
    #
    [static_pressure_front_shear]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = front
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_left_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = left
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = right
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []   
    # fix ptr
    [./fix_cptr1_x]
        type = ADDirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = ADDirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = ADDirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]

#####################################
#BCs Section
#Absorbing boundary conditions
#####################################
[BCs]
    ##non-reflecting bc
    #
    [./dashpot_front_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_front_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_front_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = front
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    #
    [./dashpot_back_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_back_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_back_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = back
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    #
    [./dashpot_top_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_top_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_top_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = top
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    #
    [./dashpot_bottom_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_bottom_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_bottom_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = bottom
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    #
    [./dashpot_left_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_left_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_left_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = left
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    #
    [./dashpot_right_x]
        type = FarmsNonReflectDashpotBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 0
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_right_y]
        type = FarmsNonReflectDashpotBC
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 1
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
    [./dashpot_right_z]
        type = FarmsNonReflectDashpotBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        velocities = 'vel_x vel_y vel_z'
        accelerations = 'accel_x accel_y accel_z'
        component = 2
        boundary = right
        beta = 0.25
        gamma = 0.5
        shear_wave_speed = 3333.33
        p_wave_speed = 5773.5
        density = 2700
    []
[]    

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
      mesh = '../static_solve/static_solve_out.e'
      system_variables = 'disp_x disp_y disp_z stress_01'
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