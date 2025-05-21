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
        3 1535 2224 2421 2640 3059 3262 4308 4886 5013 5134 5842 6377 6542 6679 7247 7535 8835 9723 10248 10795 10971 11126 11288 11363 11464 12510 12902 13203 13549 15342 16536 17556 20032 20293 20543 20694 20782 21734 21840 22299 22528 22736 22987 23244 23304 23725 24217 24305 25436 26097 26309 26365 26474 26731 26864 27322 27407 28240 28835 29746 30644 30737 30957 31012 31047 31197 31221 31325 31399 31732 32584 33508 34580 35276 35965 36303 36453 36599 38444 39053 39245 39341 39459 39863 40982 41041 41218 41336 41375 41507 41636 41760 41776 41894 42629 42790 42928 43044 43134 43398 43951 44898 44905 45604 45626 45888 45982 47507 48374 48947 49291 49695 50123 50177 50479 51198 54113 54377 54405 54530 55873 56193 56578 56582 56945 57527 58791 59109 59164 59190 59325 59504 59535 60454 60546 60815 61536 61748 62235 62777 62980 63129 63188 63257 64040 64781 64875 65075 65102 67247 68148 68194 68637 68953 69349 69835 70677 71057 71920 72005 72162 72194 72673 72930 73133 73258 73461 73538 73604 73616 73703 74017 74665 74764 74819 74985 75127 75378 75797 76111 76985 77790 77882 78744 79091 79873 80073 80747 80889 80993 81336 81719 81743 82139 83258 83333 83531 84348 84468 84693 84750 86350 86693 86854 86870 86969 87570 89399 89514 89770 89916 90219 90579 90593 90772 90865 91929 92253 92531 92625 92699 93040 93088 93482 93653 93657 94819 95072 95325 95390 95473 95715 96408 96761 97935 98087 98455 98512 98625 98719 98890 99513 99538 100249 100347 100350 100511 101171 101675 101801 101979 102137 102384 102728 102736 102755 102799 102959 103156 103425 103971 104089 104671 104903 105375 105493 105604 106426 106546 106913 107157 107258 107521 108107 108438 108847 108910 109556 109700 109775 110051 110150 110492 110736 111004 111050 111463 111724 112693 113046 113130 113140 113246 113300 113341 113639 114415 114631 114735 114903 115227 116125 116628 116669 117425 117511 117995 118281 118804 119099 119420 119526 119935 120193 120541 120826 121155 122877 123584 124541 124865 125703 126173 126234 126312 126741 127533 127596 127926 128303 128635 128843 129630 129821 129857 130054 130076 130091 130277 130587 131338 132391 132685 133511 134061 134074 134268 134348 135150 135683 135888 135997 136061 136205 136390 136410 136427 136565 136713 136769 136802 136913 137951 138137 138775 138783 139539 141696 144684 145662 147332 149288 149953 150895 153655 155886 155983 157728 159830 162733 163984 164839 165625 165950 1030 1126 1143 2901 3035 3040 4295 4513 4788 5178 5298 6030 7306 7965 8393 8631 9339 9577 9625 10572 11465 11502 11571 11658 11949 12056 12478 13169 13850 15394 15432 15675 15907 16294 16351 16415 17214 18323 18594 19290 19313 20135 20445 20639 20869 21516 21723 22268 22355 22498 22804 24884 25644 25697 25735 25740 26218 26432 26763 27847 28152 28329 28950 29096 29245 29494 30522 30722 31716 31913 32589 32694 33467 33506 33925 34594 35586 36171 36615 36759 37132 38322 38577 38979 39237 39645 39915 40019 40343 42437 43422 43530 43597 43698 43945 44317 45598 45983 46541 47017 47045 47324 47746 48707 49188 49386 49532 50218 50428 50451 50776 51229 51284 51857 52432 53229 53333 54058 54324 54346 54349 55384 56584 56649 56882 56948 57319 57470 58020 58539 58666 58875 58927 59392 59741 59778 59868 60049 60167 60843 61431 62215 63975 63983 64788 65009 65147 65558 65677 65751 65785 66100 67617 68218 69736 70190 70568 70700 70925 71806 72804 73555 74028 74335 74462 74640 74736 75116 75860 76565 76857 76997 77021 77403 77677 78525 79498 79616 79970 80039 80399 80673 80809 80870 81941 83632 83787 83831 83848 84326 84403 84695 84899 84957 84985 85081 85130 85866 86041 86389 86562 86698 87210 87358 87470 87918 87937 88070 89223 89856 89881 90077 90263 90453 91199 91583 91692 91850 91994 92112 92162 92384 92848 92854 92935 93501 93805 93982 94521 94741 95167 95304 95963 96819 96981 97294 97755 97830 98040 99661 99940 100117 100123 100429 100455 101301 101598 101604 101798 102016 102460 102790 102886 102985 103015 103570 103585 103781 103902 104195 104264 104619 104690 105216 105304 105471 105838 106007 106027 106034 106750 107214 107590 107690 108482 108666 108678 108688 108746 109228 109478 109634 110079 110123 110175 110274 110530 111243 111376 111625 113033 113247 113896 113900 113999 114503 114755 115418 115726 115856 116294 116526 116621 117436 117641 117824 118065 118339 118375 119036 119497 119911 120717 121343 121726 121913 122066 122181 123144 123195 124080 124199 124858 125018 125576 125793 126068 126463 126532 126994 127056 128459 128481 129331 129813 130399 130583 130785 131042 131104 131248 131715 131795 131918 131923 132204 132218 133468 133901 134222 134483 134519 134525 134693 135079 135081 135437 135943 136058 136096 136285 136488 136692 136960 137293 137423 137795 138187 138702 138706 138773 139361 139688 139872 143547 145483 147639 149853 155334 157065 157538 157929 159678 160848 165142 1373 3363 4030 5226 6100 8741 9729 12751 12996 13026 13364 13503 13899 14877 16215 19625 22621 24547 26713 27563 27970 28504 32179 33632 35120 35485 36240 36943 38720 38852 38871 39045 39540 39824 41488 41841 43614 44904 45444 48241 49094 49296 49344 49582 50745 51645 52170 52717 53000 53588 53647 54193 56028 56702 56726 57603 58886 59380 61153 61252 61271 61364 62856 65452 67096 68695 69787 69875 70139 71794 73023 73124 74472 75323 79343 79423 79554 79909 80017 80294 81082 81772 81878 81896 82320 83755 84852 85525 86121 86512 87359 87630 88066 88211 90046 91321 91797 91819 92650 92995 93517 94633 95554 96160 96367 97139 97641 98270 98539 98676 99129 99872 101463 101466 102528 102612 102899 104783 104835 105926 105936 107855 108056 108110 111559 111572 112088 112573 113788 113970 114455 116822 116925 118015 118727 121153 121242 122659 122855 123393 124202 124966 126649 126717 127650 128642 129384 132983 132993 133081 133724 134907 135709 135850 135906 136086 136161 137323 137562 138629 138778 141116 141290 150383 157630 161577 164534 4681 5205 5219 5606 7824 7835 12926 13825 15751 17089 17185 17240 17877 18586 20036 20655 21531 21773 23453 23492 23608 24007 24215 26066 26620 27225 27818 28258 30457 31917 32247 32961 34050 34882 34910 35773 37535 37846 37877 38419 40801 41478 43257 44633 45988 48418 50675 52634 52704 53405 55094 57188 57605 59427 59824 59913 60229 60561 60998 61614 62707 65030 65093 66625 66952 67468 67932 67977 68671 68944 69653 70318 71118 71998 72980 73643 73796 73911 74115 76062 76250 76984 77723 78218 78805 78847 79154 80185 80338 82493 83824 83862 84083 87503 90018 92319 92760 93095 93261 93383 95210 96116 96310 97980 98937 99266 100519 100710 100719 101513 101813 101848 104013 104046 104053 104975 105701 107516 107694 107917 111048 111119 111734 113458 113526 115299 117001 117638 117696 117958 118147 118440 120137 120514 122030 122144 126172 127040 127203 127897 129759 130823 131241 132097 133440 134040 134603 135251 136667 137150 137435 137597 137713 138219 138304 138587 145544 146519 148595 151338 156303 158491 162968 164696 164761 4582 8085 10372 10799 12722 13586 13742 16808 17333 17496 18147 18364 21263 22157 22590 24180 31860 34600 39777 40416 42179 46842 47363 61519 64128 68201 69306 70095 74415 74631 80926 83622 84598 85567 89309 91719 91940 92945 93385 94747 95088 96470 96989 97458 98579 100118 103986 104509 107601 107635 108737 109018 110925 111995 113695 114462 115030 117155 118156 119527 120255 122701 124340 130986 132357 132494 133711 138439 138604 140594 145040 149616 153636 162328 20099 22380 25892 33180 35508 42047 42944 44709 50096 50562 52273 54365 62147 62201 63474 64241 64662 65065 65066 68865 70501 71110 74409 74464 75328 75553 76936 77088 77464 79859 82174 85340 88894 88954 91674 92312 92700 94264 94782 97581 98982 100520 100755 101101 103204 105376 106814 107003 108958 110002 111782 112164 113480 114359 117556 122154 123150 124028 125134 127413 129900 133389 135029 136071 136367 137407 137542 137611 138766 146299 151912 156830 157133 161843 166096 166424 
        '
        subdomain_ids = '
        100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 
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
        show = 'vel_slipweakening_x vel_slipweakening_y ini_shear_sts ini_normal_sts T1_aux T2_aux jump_x_aux jump_y_aux jump_x_rate_aux jump_y_rate_aux'
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