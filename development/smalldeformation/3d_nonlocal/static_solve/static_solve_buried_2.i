#continuum damage-breakage model dynamics

#material properties
lambda_o = 32.04e9
shear_modulus_o = 32.04e9
xi_o = -0.8
xi_d = -0.85
chi = 0.8
fluid_density = 1000   
solid_density = 2700
gravity_pos = 9.81
gravity_neg = -9.81

sigma = 5e2
peak_val = 0.7
len_of_fault_strike = 14000
len_of_fault_dip = 10000
nucl_center = '0 0 -10000'

##########################################################################################################################################
#Mesh section
#FileMeshGenerator: read mesh file
#SideSetsFromNormalsGenerator: generate side sets from normals
#ExtraNodesetGenerator: generate extra nodeset - here we use it to define corner points associated with the bottom boundary
##########################################################################################################################################
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../mesh/mesh_large_buried_2.msh'
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
        coord = ' -60000 -60000 -60000;
                   60000 -60000 -60000;
                   60000 60000  -60000;
                  -60000 60000  -60000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
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

##############################################
#AuxVariables section
#gradient of damage variable (not used): alpha_grad_x, alpha_grad_y, alpha_grad_z
#initial_shear_stress_aux: initial shear stress
##############################################
[AuxVariables]
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_breakage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_I2_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_cd_aux]
        order = FIRST
        family = MONOMIAL
    []
    [initial_xi_aux]
        order = FIRST
        family = MONOMIAL
    []
    #stress function checks
    [stress_xx_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_yy_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_xy_initial]
        order = CONSTANT
        family = MONOMIAL
    []
    [stress_zz_initial]
        order = CONSTANT
        family = MONOMIAL
    []
[]

#######################################################################
#AuxKernels section
#get initial damage material property and save as an auxiliary variable
#######################################################################
[AuxKernels]
    [get_initial_stress_xx]
        type = FunctionAux 
        variable = stress_xx_initial
        function = func_neg_xx_stress
    []
    [get_initial_stress_yy]
        type = FunctionAux 
        variable = stress_yy_initial
        function = func_neg_yy_stress
    []
    [get_initial_stress_xy]
        type = FunctionAux 
        variable = stress_xy_initial
        function = func_pos_xy_stress
    []
    [get_initial_stress_zz]
        type = FunctionAux 
        variable = stress_zz_initial
        function = func_neg_zz_stress
    []
    [get_initial_damage]
        type = ADMaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
    []
    [get_initial_breakage]
        type = ADMaterialRealAux
        variable = initial_breakage_aux
        property = initial_breakage
    []
    [get_initial_I2]
        type = ADMaterialRealAux
        variable = initial_I2_aux
        property = I2_initial
    []
    [get_initial_xi]
        type = ADMaterialRealAux
        variable = initial_xi_aux
        property = xi_initial
    []
[]

###############################################################################################
#Kernels section
#ADStressDivergenceTensors: compute the divergence of stress tensor, in all three directions
###############################################################################################
[Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
    [gravity]
        type = ADGravity
        variable = disp_z
        value = ${gravity_neg}
    []
[]

###############################################################################################
#Materials section
#ADComputeSmallStrain: compute the small strain
#ADComputeDamageStressStaticDistribution: compute the stress field
#ADComputeLinearElasticStress: compute the elastic stress field
#ADComputeIsotropicElasticityTensor: compute the elasticity tensor
#ADComputeXi: compute the strain invariants ratio
#ADInitialDamageCycleSim3DPlane: define the initial damage field
#ADInitialBreakageCycleSim3DPlane: define the initial breakage field
###############################################################################################
#Block 1, 3: use continuum damage breakage model
#Block 2: use linear elastic model
###############################################################################################
[Materials]
    [density]
        type = ADGenericConstantMaterial
        prop_names = 'density'
        prop_values = ${solid_density}
    []
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        eigenstrain_names = ini_stress_to_strain
        outputs = exodus
    [] 
    ###################################################################
    #lambda_o = 30e9: lambda value
    #shear_modulus_o = 30e9: shear modulus value
    #xi_o = -0.8: strain invariants ratio: onset of damage evolution
    #chi = 0.7: ratio of solid energy and granular energy
    #xi_d = -0.9: strain invariants ratio: onset of breakage healing
    #this only applies to block 1, 3
    ###################################################################
    [stress]
        type = ADComputeDamageStressStaticDistributionDynamicCDBM
        lambda_o = ${lambda_o}
        shear_modulus_o = ${shear_modulus_o}
        xi_o = ${xi_o}
        chi = ${chi}
        xi_d = ${xi_d}
        outputs = exodus
        block = '1 3'
    [] 
    ###################################################################
    #lambda = 30e9: lambda value
    #shear_modulus = 30e9: shear modulus value
    #this only applies to block 2
    ###################################################################
    [stress_elastic]
        type = ADComputeLinearElasticStress
        output_properties = 'stress'
        outputs = exodus
        block = '2'
    []
    [elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = ${lambda_o}
        shear_modulus = ${shear_modulus_o}
    []
    ###################################################################
    [getxi]
        type = ADComputeXi
        outputs = exodus
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
        type = ADInitialDamageCycleSim3DPlane
        sigma = ${sigma}
        peak_val = ${peak_val}
        len_of_fault_strike = ${len_of_fault_strike}
        len_of_fault_dip = ${len_of_fault_dip}
        nucl_center = '${nucl_center}'
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
        type = ADInitialBreakageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0
        len_of_fault_strike = ${len_of_fault_strike}
        len_of_fault_dip = ${len_of_fault_dip}
        nucl_center = '${nucl_center}'
        output_properties = 'initial_breakage'      
        outputs = exodus
    []
    [./strain_from_initial_stress]
        type = ADComputeDamageBreakageEigenstrainFromInitialStress
        initial_stress = 'func_neg_xx_stress func_pos_xy_stress 0  
                          func_pos_xy_stress func_neg_yy_stress 0  
                          0 0 func_neg_zz_stress'
        eigenstrain_name = ini_stress_to_strain
        lambda_o = ${lambda_o}
        shear_modulus_o = ${shear_modulus_o}
        xi_o = ${xi_o}
    [../]
[]  

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]

[Executioner]
    type = Steady
    solve_type = 'NEWTON'
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-10
    nl_max_its = 20
    nl_abs_tol = 1e-12
    # this is very robust, use as default
    petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
[]  

[Outputs]
    exodus = true    
[]

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
    #
    [func_pos_zz_stress]
        type = InitialDepthDependentStress
        i = 3
        j = 3
        pos_sign = true
        fluid_density = ${fluid_density}
        rock_density = ${solid_density}
        gravity = ${gravity_pos}
        bxx = ${bxx}
        byy = ${byy}
        bxy = ${bxy}
        linear_variation_cutoff_distance = ${linear_variation_cutoff_distance}
    []
    [func_neg_zz_stress]
        type = InitialDepthDependentStress
        i = 3
        j = 3
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

[BCs]
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
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = ADFunctionNeumannBC
        variable = disp_x
        boundary = back
        function = func_neg_xy_stress
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_left_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = left
        function = func_pos_xy_stress
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right_shear]
        type = ADFunctionNeumannBC
        variable = disp_y
        boundary = right
        function = func_neg_xy_stress
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