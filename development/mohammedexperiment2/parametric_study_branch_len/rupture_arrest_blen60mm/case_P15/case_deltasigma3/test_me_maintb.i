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
        618 923 1600 1626 1697 1959 3944 4091 4838 4961 5814 7143 7409 7638 8193 9219 9545 10081 11559 12601 13035 13133 13182 13248 13651 13802 14363 14430 14598 14625 14905 14992 15343 15438 17084 17729 18090 18352 18612 20258 20290 20308 20485 20852 20941 22440 23224 23269 25612 25688 28896 29056 29774 29885 30322 31791 31852 32368 33029 33331 33373 33533 33831 34153 34159 35512 36271 36298 36960 37322 37347 38327 38518 38575 38774 38963 39913 40005 40079 41332 41625 41865 42251 42796 42886 43216 43771 43934 44169 44208 45304 45412 46171 47730 48114 48263 49603 49783 49802 50187 51379 51630 51672 51978 52101 52159 52189 52346 53424 54212 54390 54799 54978 55189 55789 56033 56487 56712 56777 57473 58120 58809 60119 60294 61325 61982 62675 63245 64049 64398 64483 64538 64569 64577 64619 65391 66045 66744 67183 67957 68669 68868 69190 69235 69724 70805 71080 71255 71363 71414 71536 71551 71737 71741 72114 72650 72716 73093 73221 73431 74022 74095 74221 75948 76327 77214 77585 77702 78274 78364 78622 78727 79471 79675 79938 80368 80803 81023 81169 81801 82230 83060 84207 84664 84936 85741 86275 86665 87992 88460 88513 89030 90057 90472 90567 90663 90730 90865 91554 91789 92154 94047 94737 95830 96037 96865 97337 97768 97898 99091 99989 100222 100259 100516 100518 100603 101454 102504 102599 102834 102883 103008 103170 103910 103924 104157 104306 104600 104686 104721 105486 105830 106177 106617 106648 106717 107444 107502 107764 108138 108393 108860 109282 109538 109601 109805 110579 110880 111792 112604 113174 113358 114187 114746 114797 115190 115375 115406 115697 115851 115913 116132 116363 116747 116778 117214 117499 117795 118471 118666 118933 119363 119979 120696 120800 121123 121492 121689 121915 122009 122231 122665 123659 123961 124011 124044 124454 124797 124920 125321 126058 126171 126201 126359 126514 126725 126802 126950 127581 127791 127827 128401 128530 129240 129322 129806 130128 130331 130856 130910 131288 132223 132765 133690 133781 133796 133957 134488 134677 134850 134933 135029 135124 135501 135540 135835 136145 136152 136315 136699 137359 137403 137651 137767 138046 138158 138426 138527 141935 142442 143300 143911 144471 144892 145132 145197 146213 147218 148836 149140 149617 150147 151846 152709 153381 153952 154711 155530 155706 156010 156350 157378 158408 158657 161331 161550 163906 165093 165905 991 1093 1108 2842 2978 2984 4237 4451 4724 5108 5254 6026 7274 7938 8369 8604 9541 10512 11501 11529 11585 11661 11916 12026 12511 13174 13862 15403 15430 15666 15922 16307 16379 16437 17253 18329 18602 19310 20177 20454 20641 20864 21545 21774 22308 22398 22525 22829 24787 25587 25643 25692 25699 26123 26335 26678 27814 28147 28316 29025 29156 29321 30600 30775 31676 31885 32494 32604 33479 33547 33947 34630 36311 36745 36899 37285 38528 38807 39202 39447 39833 40099 40206 40538 42609 43625 43734 43816 43915 44114 44527 46180 46767 47244 47515 47898 48865 49315 49510 49689 50386 50615 50631 50901 51340 51393 51933 52479 53175 54029 54306 54330 54334 55335 56474 56534 56775 56849 57261 57403 57978 58493 58618 58847 58911 59724 60038 60146 60813 61380 62209 63985 63996 64807 65003 65117 65508 65628 65698 65731 66060 67562 68131 69668 70102 70444 70568 70761 71599 72581 73377 73854 74429 74543 74929 75632 76348 76620 76753 76770 77087 77357 78217 79164 79281 79674 79743 80051 80317 80465 80522 81659 83373 83537 83576 83588 84100 84160 84436 84687 84743 84769 84872 85760 86231 86559 87082 87227 87342 87780 87794 87921 89673 89711 89919 90134 90335 91084 91446 91546 91693 91832 91966 92023 92263 92682 92686 92754 93328 93653 93802 94317 94951 95094 95755 96618 96804 97103 97580 97647 97847 99508 99816 100002 100006 100356 100391 101227 101539 101543 101748 101975 102404 102723 102819 102912 102945 103502 103514 103719 103863 104178 104243 104597 104662 105225 105312 105488 105895 106075 106084 106094 106771 107229 107604 107708 108511 108689 108699 108711 108761 109229 109461 109630 110138 110193 110296 110547 111230 111376 111619 113015 113225 113827 113832 113928 114439 114701 115390 115707 115821 116291 116489 116595 117448 117663 117833 118084 118379 118407 119070 120007 120721 121343 121781 121962 122094 122205 123166 124128 124793 124972 125481 125681 125968 126356 126422 126874 126935 128207 128233 129090 129583 130112 130321 130513 130776 130832 130975 131461 131531 131640 131644 131922 133242 133638 133913 134180 134229 134417 134787 134789 135136 135766 135800 135985 136186 136398 136679 137023 137146 137515 137891 138376 138384 138428 139003 139343 139524 143011 144859 147057 149221 154715 156450 156888 157272 159043 160212 164500 138 159 1266 2799 2884 3006 3066 4068 4131 4161 5794 6086 6572 9777 9858 10755 12190 12734 14510 16964 20957 25763 26007 26121 27755 28836 28931 33405 33671 34577 35702 36148 37023 38033 38176 38248 38294 40118 41089 41218 41256 41549 41926 42209 42827 45152 46994 47062 47082 47422 47951 48278 48923 49741 51062 51221 51714 53198 53368 53793 54806 55434 55797 56069 56271 59098 60228 60967 61840 62590 62759 63269 63742 64435 64785 65313 66386 69084 69164 69208 70308 71787 71880 73491 73754 75144 75686 77590 77776 78221 78876 79752 83379 84040 85090 85992 86433 87260 87776 87983 88789 88846 88954 89301 90171 90713 92196 92672 94974 95150 95593 96161 97901 98089 98435 100108 100459 100953 101036 101209 102027 102360 102442 102946 102999 103261 104076 104189 104458 106461 107832 108808 111455 112010 112287 112314 112747 112873 113470 114080 114446 116668 118113 119017 120324 122727 122978 123812 124499 125607 126660 129259 129850 129956 130947 131570 133039 135711 135859 136340 137585 138238 138325 138451 139357 141384 142586 147185 160477 160735 163309 163964 4613 5146 5162 5580 7786 7798 12953 13831 15750 17139 17230 17270 17848 18593 20078 20654 21569 21813 23414 23445 23566 23923 24126 26006 26548 27151 27782 28244 30533 31888 32181 32926 34079 34938 34964 35891 37696 38003 38026 38639 40972 41627 43452 44863 46190 48581 50819 52647 52704 53347 55021 57128 57546 59439 59822 59916 60209 60548 60955 61567 62668 65033 65076 66583 66872 67424 67854 67897 68600 68853 69582 70208 70959 71778 72780 73455 73605 73733 73923 75829 76015 76735 77392 77892 78468 78519 78839 79877 79999 82198 83572 83608 83852 87370 89854 92215 92604 92913 93077 93206 94998 95920 96130 97793 98719 99050 100445 100631 100637 101437 101764 101810 103972 104009 104013 104941 105740 107537 107712 107948 111064 111123 111744 113419 113483 115265 117006 117657 117709 117975 118174 118463 120210 120552 122066 122166 126106 126922 127053 127701 129527 130548 130972 131813 133211 133754 134315 134958 136380 136863 137155 137309 137432 137915 137989 138258 144901 145917 148032 150704 155716 157902 162310 164046 164114 669 1669 2145 8422 11503 13227 13703 18940 19512 19940 20612 21177 27556 28581 30873 31583 32247 33026 33137 33735 34682 34685 36166 36678 37938 38558 40638 42135 42268 42854 42881 43936 45233 46044 46924 47345 47446 48575 49334 49997 51060 52022 52407 52705 53270 54517 55591 57879 58315 59229 62262 66997 68308 69149 69250 69289 69674 71387 72064 72492 73711 74919 75052 75657 76270 76378 76992 77583 78197 78696 80110 80603 81362 84216 84537 84955 84978 86225 86351 86464 87139 87343 87994 88116 88309 88440 89223 89723 92740 93290 93376 94981 95674 95869 97095 97674 97737 98024 98189 98407 99349 99417 100612 101220 101249 101351 102437 104079 104516 104638 104923 106217 107970 112548 117085 118314 119003 119216 119945 122612 122661 125004 125518 125904 127753 129871 131683 132969 133676 134304 134742 135359 135860 137882 137937 138402 138424 138439 139986 141562 141968 144612 153134 158681 810 1063 5102 9016 10775 12173 13830 14550 14841 17177 17414 17525 20994 23576 23645 24329 25487 27672 30282 30841 31302 36911 36918 37263 39928 40258 41993 42732 45585 46092 48538 48779 48828 48905 49021 49425 49515 50273 51975 52201 53349 53856 55501 56015 57171 59224 59444 59911 60126 60783 62813 64537 65409 68232 70057 70339 71117 71311 71381 71731 72304 74555 75043 75326 78137 78769 79127 79625 80154 80556 81363 84050 84617 84891 85091 85583 92908 94481 94507 96722 97089 97204 97380 97463 97562 101239 101242 101639 101874 102964 103166 104398 105251 105732 108292 109174 109391 109941 111514 111817 111981 113049 113055 113468 113699 114568 115429 116326 116508 116924 118560 118814 121183 122277 122716 122883 122945 123596 124874 125409 126493 128930 129700 130124 131858 131917 132531 132871 133608 133760 134479 135357 137022 137319 137503 137766 138233 138415 141997 147691 149364 150363 158495 160427 
        '
        subdomain_ids = '
        100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 
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
    [checkpoint]
        type = Checkpoint
        time_step_interval = 200
        num_files = 2
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