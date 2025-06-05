# G (GPa): 2.170 
# rho (kg/m^3): 1180.000 
# nu (-): 0.370 
# mus (-): 0.700 
# mud (-): 0.100 
# Dc (mu m): 32.000 
# sigmaxx (MPa): -5.550 
# sigmayy (MPa): -15.000 
# sigmaxy (MPa): 0.000 
# use dx (mm): 0.900 
# dx (mm): 1.132 
# use dx < dx_max (mm): true 
# -------------------------main fault--------------------------
# tau_S_main (MPa): -4.007 
# sigma_N_main (MPa): -12.779 
# tau_S_strength_main (MPa): 8.945 
# Lfric_main (mm): 9.057 
# ------------------------branch fault-------------------------
# tau_S_branch (MPa): 4.007 
# sigma_N_branch (MPa): -7.771 
# tau_S_strength_branch (MPa): 5.440 
# Lfric_branch (mm): 9.057 
# -------------------------------------------------------------
# numofelem (-): 11.000 use 15 elems

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../../mesh/me.msh'
    []
    #element_ids = ids - 1
    [./subdomain_id] 
        input = msh
        type = SubdomainPerElementGenerator 
        element_ids = '
        863 2104 3508 4233 5387 6008 6158 6658 7405 7856 8020 8147 8512 8714 9330 9913 9988 11014 11057 12446 12727 15434 17787 19356 19793 19968 20889 21070 21797 22951 22976 23010 23761 26038 26490 26579 27706 27942 28413 29073 29102 29984 30302 30659 30902 34520 34723 34792 35043 35948 36761 37337 37654 38623 38650 39222 39355 39776 39786 40865 41028 41925 42877 43713 44068 44098 44965 46822 46875 48320 48727 49074 49740 49794 49829 49857 50196 51174 51303 51723 52544 53251 54560 55106 56152 56866 57260 59062 59213 59574 59771 59946 60323 62552 63058 63601 64614 64711 65254 65389 65425 65513 65622 65851 65896 67729 68110 68612 68772 69426 69706 70000 70262 70467 70703 70756 71205 71364 71860 71957 72581 72636 73007 73551 74333 74394 74773 75506 75733 76132 76730 77012 77064 77469 77801 78421 78630 79060 79126 79149 79331 79785 79855 80093 80352 80631 81135 81170 81492 81928 82343 82469 82821 83180 83285 83642 83727 84077 84135 84833 85294 85578 85872 85874 86244 86814 87005 87352 87850 87898 87926 88153 88288 88737 88743 89081 89261 89628 90670 90690 90966 91194 91382 91843 92100 92171 92236 92744 93028 93543 94102 94408 94448 94927 96043 96259 96759 97037 97283 97300 97366 97593 98022 98416 99068 99553 99699 99775 99919 100456 101246 101613 101787 101789 101869 102063 102551 103253 103678 104226 104228 104496 104845 105528 106145 107032 107138 107288 107373 107483 107655 109594 110288 110855 110965 111149 111469 111604 111759 111982 112208 112476 112523 113570 114269 114578 114923 116639 116847 117265 117442 118311 118408 119437 120742 121407 121544 122359 122710 122796 123130 124106 125115 126359 126710 126796 126890 127383 127389 128610 128765 130357 130914 131257 132266 132441 132953 133327 133449 133465 133715 134016 134423 134596 135176 135217 135245 135891 136214 136452 136773 137041 137057 137207 137284 137326 137919 137955 138039 138148 138187 138242 138506 139085 139304 139384 139499 139554 139599 139664 139715 139901 139927 140465 140808 143241 143426 145689 147923 148273 148470 151224 151863 156605 157229 158405 162490 163269 164509 164541 164659 166609 167672 168493 1098 1211 1230 3058 3212 3219 4531 4739 5029 5408 5559 6344 7652 8306 8751 8985 9919 10952 11904 11947 12014 12388 12486 12924 14278 15821 16124 16375 16747 16812 16872 17674 18712 18981 19679 20537 20872 21047 21272 21950 22166 22766 22915 23240 25183 25939 25983 26031 26045 26491 26712 27044 28154 28495 29277 29405 30845 31017 31964 32162 32803 32918 33706 33751 34207 34886 36990 37139 37530 38728 39006 39448 39679 40095 40370 40468 40796 42873 43893 44026 44095 44802 46466 47038 47537 47821 48205 49564 49780 49953 50662 50858 50882 51160 51625 52170 52776 53527 54410 54666 54685 54688 55721 56821 56891 57133 57196 57567 57730 58299 58860 58999 59231 59298 60088 60390 60498 61148 61777 62556 64243 64253 65099 65296 65833 65971 66037 66069 66424 68552 70061 70528 70889 71003 71208 72105 73136 73881 74362 74936 75048 75425 76938 77377 77399 77766 78027 78943 79897 80026 80388 80446 80775 81038 81163 81230 82377 84120 84298 84339 84354 84933 85201 85423 85478 85499 85611 86493 87040 87349 87944 88127 88246 88708 88724 88869 90747 90988 91197 91384 92195 92559 92665 92799 92951 93102 93145 93355 93818 93821 93897 94498 94823 94979 95507 96171 96298 96992 97807 97973 98285 98731 98808 99023 100667 100961 101170 101175 101533 102366 102661 102666 102863 103083 103525 103842 103937 104035 104056 104589 104600 104811 104951 105267 105328 105655 106276 106373 106528 106927 107101 107112 107120 107833 108291 108646 108757 109582 109758 109768 109787 109836 110311 110556 110718 111245 111299 111661 112356 112496 112746 114184 114397 115011 115019 115118 115874 116529 116831 116948 117441 117669 117781 118650 118861 119057 119302 119594 119630 120249 121197 121983 123075 123254 123393 123510 124484 125484 126316 126828 127035 127332 127729 127788 128233 128283 129625 129651 130524 131010 131576 131776 131983 132257 132315 132454 132962 133031 133149 133154 133452 134745 135161 135729 135765 135968 136332 136333 136672 137290 137325 137505 137696 137890 138150 138475 138614 138943 139345 139859 139868 139918 140473 140818 140994 144651 146537 148852 151157 156834 159056 159444 161308 162591 166988 2952 5381 9019 10023 10079 16181 16883 17539 18497 21552 22105 22191 22200 23724 26466 31847 32252 32849 33272 34002 34125 35928 40323 40782 41251 42355 42806 44065 44992 45361 47179 47207 48197 48474 48999 54225 54889 56739 56953 59555 59687 61915 65955 66080 67966 68273 68425 68911 69211 70424 71366 73126 73308 75163 77793 77798 79057 79391 79924 82585 83089 84871 85428 85776 86264 86512 86981 87335 87992 88039 88175 88332 88852 88976 89297 90731 90803 93598 93841 93857 93903 93905 94163 95327 95930 97067 97077 97845 97865 97926 98754 98944 99343 99539 100199 100232 101742 103223 105018 105208 105569 105712 105871 105903 107985 107995 108335 108377 108631 108707 108884 108957 109336 110548 110803 111292 111458 112156 113132 113693 113741 115754 116328 116650 117099 117247 117476 117904 118674 120529 121610 123538 123928 126154 126203 127465 128650 131444 131952 132463 133896 134638 134768 135978 136105 137102 137148 137249 137291 137632 138142 139022 139372 139783 139828 142339 143942 144701 145027 145268 146118 146718 146944 147009 148639 150108 152864 153109 154592 155221 157839 164966 166915 4920 5448 5472 5896 8156 8171 13375 14243 16207 17558 17643 17692 18265 18974 20443 21057 21972 22201 23841 23879 23995 24351 24562 26355 26917 27513 28115 28598 30784 32167 32507 33174 34346 35201 35220 36161 37939 38222 38250 38839 41280 41910 43701 45101 46472 48822 51064 52972 53023 53704 55408 57455 57846 59800 60171 60275 60559 60881 61319 61974 63022 65317 65370 66944 67247 67760 68239 68283 68973 69226 69972 70645 71415 72312 73299 73961 74116 74246 74434 76401 76577 77369 78067 78573 79194 79237 79558 80582 80716 82954 84332 84373 84615 88274 90912 93300 93733 94068 94253 94384 96211 97139 97328 98972 99917 100250 101590 101769 101773 102567 102880 102907 105057 105088 105094 106027 106774 108579 108760 108974 112183 112239 112885 114589 114658 116399 118191 118858 118914 119193 119395 119699 121427 121792 123360 123467 127467 128268 128395 129111 130937 132020 132450 133324 134721 135283 135846 136500 137875 138322 138626 138758 138864 139377 139451 139705 146591 147665 149859 152719 157827 160069 164761 166530 166601 1687 5423 20827 29402 30474 31296 36862 37551 58336 61035 83450 84197 97825 103582 103647 105556 107134 117269 117379 119173 121805 127442 140156 140642 146082 158995 159597 160981 3078 9437 10322 13633 24165 34543 37564 38867 46270 46410 47798 49344 58360 73102 73534 74303 86916 102970 105866 120986 121905 127628 129892 131503 139096 142840 151946 154603 159527 
        '
        subdomain_ids = '
        100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 
        '
    []
    [./split]
        input = subdomain_id
        type = BreakMeshByBlockGenerator
        block_pairs = '100 200; 300 200; 300 100'
        split_interface = true
    []
    [./sidesets]
        input = split
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0'
        new_boundary = 'left right bottom top'
    []
[]

[GlobalParams]
    displacements = 'disp_x disp_y'
    q = 0.2
    # Dc = 2.4e-5 #mm -> m
    # area = 8e-4
    out_of_plane_strain = strain_zz
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
    [./strain_zz]
    []
[]

[AuxVariables]
    [./resid_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_y]
        order = FIRST
        family = LAGRANGE
    []
    [./nodal_area]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_slipweakening_y]
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
    [./vel_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./vel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    #import from csv
    [./ini_shear_sts]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_normal_sts]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_d]
        order = CONSTANT
        family = MONOMIAL
    []
    [./D_c]
        order = CONSTANT
        family = MONOMIAL
    []   
    #output special case parameters
    [./flag_track_opening]
        order = CONSTANT
        family = MONOMIAL
    []
    [./flag_track_activecase]
        order = CONSTANT
        family = MONOMIAL
    []
    [./jump_track_opening]
        order = CONSTANT
        family = MONOMIAL
    []
    [./jump_track_reversal]
        order = CONSTANT
        family = MONOMIAL
    []
    [./jump_effective]
        order = CONSTANT
        family = MONOMIAL
    []
    #output csv files
    [./jump_x_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [./jump_y_rate_aux]
        order = FIRST
        family = MONOMIAL
    []
    [./jump_x_aux]
        order = FIRST
        family = MONOMIAL
    []
    [./jump_y_aux]
        order = FIRST
        family = MONOMIAL
    [] 
    [./T1_aux]
        order = FIRST
        family = MONOMIAL
    []
    [./T2_aux]
        order = FIRST
        family = MONOMIAL
    []
[]

[Physics/SolidMechanics/CohesiveZone]
    [./czm_ik]
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        strain = SMALL
        generate_output='traction_x traction_y jump_x jump_y'
    [../]
[]

[Physics/SolidMechanics/QuasiStatic]
    [plane_stress]
      planar_formulation = WEAK_PLANE_STRESS
      strain = SMALL
      generate_output = 'stress_xx stress_yy stress_xy'
    []
[]

[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_slipweakening_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_slipweakening_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    #jump rate (element)
    [TJump_rate]
        type = FDCompVarRate
        variable = jump_x_rate_aux
        coupled = jump_x_aux
        execute_on = 'TIMESTEP_BEGIN'
    []
    [NJump_rate]
        type = FDCompVarRate
        variable = jump_y_rate_aux
        coupled = jump_y_aux
        execute_on = 'TIMESTEP_BEGIN'
    []
    #
    [Displacment_x]
        type = CompVar
        variable = disp_slipweakening_x
        coupled = disp_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Displacement_y]
        type = CompVar
        variable = disp_slipweakening_y
        coupled = disp_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_x]
        type = CompVar
        variable = resid_slipweakening_x
        coupled = resid_x
        execute_on = 'TIMESTEP_BEGIN'
    []
    [Residual_y]
        type = CompVar
        variable = resid_slipweakening_y
        coupled = resid_y
        execute_on = 'TIMESTEP_BEGIN'
    []
    #func->aux
    [ShearStress]
        type = FunctionAux
        variable = ini_shear_sts
        function = func_ini_shear_sts
        execute_on = 'initial TIMESTEP_BEGIN'
    []
    [NormalStress]
        type = FunctionAux
        variable = ini_normal_sts
        function = func_ini_normal_sts
        execute_on = 'initial TIMESTEP_BEGIN'
    []
    [StaticFricCoeff]
        type = FunctionAux
        variable = mu_s
        function = func_mus
        execute_on = 'initial TIMESTEP_BEGIN'
    []
    [DynamicFricCoeff]
        type = FunctionAux
        variable = mu_d
        function = func_mud
        execute_on = 'initial TIMESTEP_BEGIN'
    []
    [CharacteristicLength]
        type = FunctionAux
        variable = D_c
        function = func_dc
        execute_on = 'initial TIMESTEP_BEGIN'
    []
    #fault length
    [fault_len]
        type = ConstantAux
        variable = nodal_area
        value = 9e-4
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #output material properties
    [flag_track_opening]
        type = MaterialRealAux
        property = flag_track_opening
        variable = flag_track_opening
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [flag_track_activecase]
        type = MaterialRealAux
        property = flag_track_activecase
        variable = flag_track_activecase
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [jump_track_opening]
        type = MaterialRealAux
        property = jump_track_opening
        variable = jump_track_opening
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [jump_track_reversal]
        type = MaterialRealAux
        property = jump_track_reversal
        variable = jump_track_reversal
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [jump_effective]
        type = MaterialRealAux
        property = jump_effective
        variable = jump_effective
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    #output properties in csv files
    [jump_x]
        type = MaterialRealAux
        property = jump_x
        variable = jump_x_aux
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [jump_y]
        type = MaterialRealAux
        property = jump_y
        variable = jump_y_aux
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [T1]
        type = MaterialRealAux
        property = T1
        variable = T1_aux
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [T2]
        type = MaterialRealAux
        property = T2
        variable = T2_aux
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
        execute_on = 'INITIAL TIMESTEP_BEGIN'
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
[]

[Materials]
    #take PMMA: poisson's ratio 0.37, shear modulus 2.17GPa
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 6.176e9
        shear_modulus = 2.17e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 1180 #kg/m^3
    []
    #SlipWeakeningMultifaults ONLY supports TRIA currently!
    [./czm_mat]
        type = SlipWeakeningMultifaultsGeneralizedPropCSV
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        nodal_area = nodal_area
        D_c = D_c
        mu_d = mu_d
        mu_s = mu_s
        ini_shear_sts = ini_shear_sts
        ini_normal_sts = ini_normal_sts
        boundary = 'Block100_Block200 Block200_Block300 Block100_Block300'
    [../]
[]
    
[Executioner]
    #dt = prefactor * dx / pressure wave speed = 0.0009 / 2985
    type = Transient
    dt = 3.015e-8 
    end_time = 10.0
    num_steps = 16000
    [TimeIntegrator]
        type = FarmsCentralDifference
        solve_type = lumped
    []
[]
    
[Outputs]
    [exodus]
        type = Exodus
        execute_on = 'timestep_end'
        time_step_interval = 40
        show = 'disp_slipweakening_x disp_slipweakening_y vel_slipweakening_x vel_slipweakening_y ini_shear_sts ini_normal_sts T1_aux T2_aux jump_x_aux jump_y_aux jump_x_rate_aux jump_y_rate_aux'
    []
    [csv]
        type = CSV
        execute_on = 'timestep_end'
        time_step_interval = 40
    []
[]

[MultiApps]
    [./sub_app]
        type = TransientMultiApp
        positions = '0 0 0'
        input_files = 'test_me_subtb.i'
        execute_on = 'TIMESTEP_BEGIN'
    [../]
[]

[Transfers]
    [pull_resid]
        type = MultiAppCopyTransfer
        from_multi_app = sub_app
        source_variable = 'resid_sub_x resid_sub_y'
        variable = 'resid_x resid_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
    [push_disp]
        type = MultiAppCopyTransfer
        to_multi_app = sub_app
        source_variable = 'disp_x disp_y'
        variable = 'disp_sub_x disp_sub_y'
        execute_on = 'TIMESTEP_BEGIN'
    []
[]

[Functions]
    [func_ini_shear_sts]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'properties'
        read_type = 'element'
        column_number = '0'
    []
    [func_ini_normal_sts]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'properties'
        read_type = 'element'
        column_number = '1'
    []
    [func_mus]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'properties'
        read_type = 'element'
        column_number = '2'
    []
    [func_mud]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'properties'
        read_type = 'element'
        column_number = '3'
    []
    [func_dc]
        type = PiecewiseConstantFromCSV
        read_prop_user_object = 'properties'
        read_type = 'element'
        column_number = '4'
    []
[]

[UserObjects]
    [./properties]
        type = PropertyReadFile
        prop_file_name = 'faults_elem_properties.csv'
        read_type = 'element'
        nprop = 5
    []
[]

[VectorPostprocessors]
    [main_fault]
        type = SideValueSampler
        variable = 'T1_aux T2_aux jump_x_aux jump_y_aux jump_x_rate_aux jump_y_rate_aux' 
        boundary = 'Block100_Block200 Block200_Block300'
        sort_by = x
    []
    # [main_fault_2]
    #     type = SideValueSampler
    #     variable = 'T1_aux T2_aux jump_x_aux jump_y_aux jump_x_rate_aux jump_y_rate_aux' 
    #     boundary = 'Block200_Block300'
    #     sort_by = x
    # []
    [branch_fault]
        type = SideValueSampler
        variable = 'T1_aux T2_aux jump_x_aux jump_y_aux jump_x_rate_aux jump_y_rate_aux' 
        boundary = 'Block100_Block300'
        sort_by = x
    []    
[]