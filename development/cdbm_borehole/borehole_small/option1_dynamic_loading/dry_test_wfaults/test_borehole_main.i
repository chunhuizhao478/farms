##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.4 [Check!]
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677

# under the stress state, 
# -135MPa 70MPa; 70MPa -120e6
# the principal stress is 198e6 

# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.5
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.4) * 198e6) = 233.6m

# Use mesh size = 25m, resolved by more than 7 elements [Check!]

# Use 1% overstress

# Use dt = 2e-4
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../../../meshgenerator/cdbm/borehole_cracks/borehole_wofaults.msh'
  []
  [./subdomain_id] 
  input = msh
  type = ElementSubdomainIDGenerator 
  element_ids = '
  92 1468 1724 4484 6016 7104 8340 8472 10656 10832 11052 12644 12956 14104 21348 22760 24556 29000 30412 32564 32592 33076 33480 35652 36364 37544 39040 39212 40420 40660 41292 43532 43604 44316 47160 47208 47676 47912 49836 50124 50584 51524 51556 51768 53496 54572 55972 56292 57780 62900 63292 64512 1204 1516 1552 1568 2856 5496 5884 6708 8484 10900 10928 13840 14116 14992 17496 21056 21444 22120 22468 26216 28420 28660 29804 29832 30632 30884 31544 32436 33948 34084 35068 38416 39168 40648 41040 41788 41988 42380 45632 47052 47944 48612 48960 49808 50000 50816 52608 52704 53864 54080 55488 56544 64292 412 1584 2736 6640 7984 8476 9324 11008 12544 14108 18540 21324 22672 23856 24108 26836 28176 29244 32516 32544 33756 34968 35004 36060 36988 38020 39528 40496 40960 41112 41224 41744 45956 46460 46516 47132 48132 48232 49224 49744 50956 51220 51864 54356 55208 55676 56024 56332 57476 62492 63688 96 1472 1728 7108 8464 10652 11056 12640 12952 14096 16028 17236 19484 22756 24252 24552 26140 27040 27448 28996 29636 32596 33072 33484 35648 36360 37548 37736 37816 39044 39208 40308 41288 42440 42904 47164 47204 47672 47916 50120 50580 51528 51552 53380 54568 55468 56288 56352 57776 62904 64808 368 1208 1556 2572 2784 3152 5500 5888 6532 6912 8468 9024 10904 10924 11312 11872 14100 17492 21780 22464 23412 26172 26212 28416 28880 29828 30064 30280 30880 32432 33448 33944 34088 35056 35980 38264 38468 39576 40652 41036 43168 45636 47940 48616 48956 50004 52564 52612 53860 55824 56540 416 1588 2740 6636 7988 8480 9328 11004 12540 14112 18536 21320 22668 23860 24104 26832 28172 29332 32517 32545 33752 34964 35000 36056 36992 38016 39524 40492 40956 41108 41220 41820 45960 46456 46552 47128 48136 48228 49228 49740 50960 51216 51860 54360 55212 55672 56028 56296 57472 62488 64400 1205 1517 1553 1569 2857 5497 5885 6709 8485 10901 10929 13841 14117 14993 17497 21057 21445 22121 22469 26217 28421 28661 29805 29833 30633 30885 31545 32437 33949 34085 35069 38417 39169 40649 41041 41789 41989 42381 45633 47053 47945 48613 48961 49809 50001 50817 52609 52705 53865 54081 55489 56545 64293 93 1469 1725 4485 6017 7105 8341 8473 10657 10833 11053 12645 12957 14105 21349 22761 24557 29001 30413 32565 32593 33077 33481 35653 36365 37545 39041 39213 40421 40661 41293 43533 43605 44317 47161 47209 47677 47913 49837 50125 50585 51525 51557 51769 53497 54573 55973 56293 57781 62901 63293 64513 97 1473 1729 7109 8465 10653 11057 12641 12953 14097 16029 17237 19485 22757 24253 24553 26141 27041 27449 28997 29637 32597 33073 33485 35649 36361 37549 37737 37817 39045 39209 40309 41289 42441 42905 47165 47205 47673 47917 50121 50581 51529 51553 53381 54569 55469 56289 56353 57777 62905 64809 413 1585 2737 6641 7985 8477 9325 11009 12545 14109 18541 21325 22673 23857 24109 26837 28177 29245 32518 32546 33757 34969 35005 36061 36989 38021 39529 40497 40961 41113 41225 41745 45957 46461 46517 47133 48133 48233 49225 49745 50957 51221 51865 54357 55209 55677 56025 56333 57477 62493 63689 417 1589 2741 6637 7989 8481 9329 11005 12541 14113 18537 21321 22669 23861 24105 26833 28173 29333 32519 32547 33753 34965 35001 36057 36993 38017 39525 40493 40957 41109 41221 41821 45961 46457 46553 47129 48137 48229 49229 49741 50961 51217 51861 54361 55213 55673 56029 56297 57473 62489 64401 369 1209 1557 2573 2785 3153 5501 5889 6533 6913 8469 9025 10905 10925 11313 11873 14101 17493 21781 22465 23413 26173 26213 28417 28881 29829 30065 30281 30881 32433 33449 33945 34089 35057 35981 38265 38469 39577 40653 41037 43169 45637 47941 48617 48957 50005 52565 52613 53861 55825 56541 1206 1518 1554 1570 2858 5498 5886 6710 10902 10930 13842 14118 14994 17498 21058 21446 22122 22470 26218 28422 28662 29806 29834 30634 30886 31546 32438 33950 34086 35070 38418 39170 40650 41042 41790 41990 42382 45634 47054 47946 48614 48962 49810 50002 50818 52610 52706 53866 54082 55490 56546 62446 64294 94 1470 1726 4486 6018 7106 8342 10658 10834 11054 12646 12958 14106 21350 22762 24558 29002 30414 32566 32594 33078 33482 35654 36366 37546 39042 39214 40422 40662 41294 43534 43606 44318 47162 47210 47678 47914 49838 50126 50586 51526 51558 51770 53498 54574 55974 56294 57782 60482 62902 63294 64514 98 1474 1730 7110 8466 10654 11058 12642 12954 14098 16030 17238 19486 22758 24254 24554 26142 27042 27450 28998 29638 32598 33074 33486 35650 36362 37550 37738 37818 39046 39210 40310 41290 42442 42906 47166 47206 47674 47918 50122 50582 51530 51554 53382 54570 55470 56290 56354 57778 62906 64810 414 1586 2738 6642 7986 8478 9326 11010 12546 14110 18542 21326 22674 23858 24110 26838 28178 29246 32520 32548 33758 34970 35006 36062 36990 38022 39530 40498 40962 41114 41226 41746 45958 46462 46518 47134 48134 48234 49226 49746 50958 51222 51866 54358 55210 55678 56026 56334 57478 62494 63690 418 1590 2742 6638 7990 8482 9330 11006 12542 14114 18538 21322 22670 23862 24106 26834 28174 29334 32521 32549 33754 34966 35002 36058 36994 38018 39526 40494 40958 41110 41222 41822 45962 46458 46554 47130 48138 48230 49230 49742 50962 51218 51862 54362 55214 55674 56030 56298 57474 62490 64402 370 1210 1558 2574 2786 3154 5502 5890 6534 6914 8470 9026 10906 10926 11314 11874 14102 17494 21782 22466 23414 26174 26214 28418 28882 29830 30066 30282 30882 32434 33450 33946 34090 35058 35982 38266 38470 39578 40654 41038 43170 45638 47942 48618 48958 50006 52566 52614 53862 55826 56542 95 1471 1727 4487 6019 7107 8343 10659 10835 11055 12647 12959 14107 21351 22763 24559 29003 30415 32567 32595 33079 33483 35655 36367 37547 39043 39215 40423 40663 41295 43535 43607 44319 47163 47211 47679 47915 49839 50127 50587 51527 51559 51771 53499 54575 55975 56295 57783 60483 62903 63295 64515 1207 1519 1555 1571 2859 5499 5887 6711 10903 10931 13843 14119 14995 17499 21059 21447 22123 22471 26219 28423 28663 29807 29835 30635 30887 31547 32439 33951 34087 35071 38419 39171 40651 41043 41791 41991 42383 45635 47055 47947 48615 48963 49811 50003 50819 52611 52707 53867 54083 55491 56547 62447 64295 415 1587 2739 6643 7987 8479 9327 11011 12547 14111 18543 21327 22675 23859 24111 26839 28179 29247 32522 32550 33759 34971 35007 36063 36991 38023 39531 40499 40963 41115 41227 41747 45959 46463 46519 47135 48135 48235 49227 49747 50959 51223 51867 54359 55211 55679 56027 56335 57479 62495 63691 99 1475 1731 7111 8467 10655 11059 12643 12955 14099 16031 17239 19487 22759 24255 24555 26143 27043 27451 28999 29639 32599 33075 33487 35651 36363 37551 37739 37819 39047 39211 40311 41291 42443 42907 47167 47207 47675 47919 50123 50583 51531 51555 53383 54571 55471 56291 56355 57779 62907 64811 371 1211 1559 2575 2787 3155 5503 5891 6535 6915 8471 9027 10907 10927 11315 11875 14103 17495 21783 22467 23415 26175 26215 28419 28883 29831 30067 30283 30883 32435 33451 33947 34091 35059 35983 38267 38471 39579 40655 41039 43171 45639 47943 48619 48959 50007 52567 52615 53863 55827 56543 419 1591 2743 6639 7991 8483 9331 11007 12543 14115 18539 21323 22671 23863 24107 26835 28175 29335 32523 32551 33755 34967 35003 36059 36995 38019 39527 40495 40959 41111 41223 41823 45963 46459 46555 47131 48139 48231 49231 49743 50963 51219 51863 54363 55215 55675 56031 56299 57475 62491 64403 
  '
  subdomain_ids = '
  101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 101 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 201 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 102 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 202 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 103 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 203 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 104 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 204 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 105 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 205 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 106 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 206 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 107 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 207 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 108 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 208 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 109 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 209 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 110 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 210 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 111 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 211 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 112 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212 212
  '
[]
[./split]
  input = subdomain_id
  type = BreakMeshByBlockGenerator
  block_pairs = '101 201; 102 202; 103 203; 104 204; 105 205; 106 206; 107 207; 108 208; 109 209; 110 210; 111 211; 112 212;'
  split_interface = true
[]
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
    #damping ratio
    q = 0.2
  
    #characteristic length (m)
    Dc = 0.4
  
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 3.204e10
      
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 3.204e10
  
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.6
  
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.6
  
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
  
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8
  
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
  
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
  
    ##Compute gamma_damaged_r, xi_1
    #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
    #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
    #check struct_param.m 
  
    #coefficient of damage solid modulus
    gamma_damaged_r = 4.0464e+10
  
    #critical point of three phases (strain invariants ratio vs damage)
    xi_1 = 0.936229
  
    ##Compute parameters in granular states
    #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
    #check struct_param.m
  
    #--------------------------------------------------------------------------------#
    #Note: "computeAlphaCr" needs to change every time the related parameters changed
    #--------------------------------------------------------------------------------#
  
    # #coefficients
    # chi = 0.75
    a0 = 6.37464e+09
    a1 = -2.63965e+10
    a2 = 3.45711e+10
    a3 = -1.45789e+10
  
    #diffusion coefficient #for structural stress coupling
    D = 0
  
    #pressure analytical solution
    #Reference: Injection-induced seismicity: Poroelastic and earthquake nucleation effects (P. Segall1 and S. Lu2)
    # effec_sts_coeff = 0.31 #
    # flux_q = 5e2 #kg/s
    # density_rho_0 = 1e3 #kg/m^3
    # permeability_k = 3e-12 #m^2
    # viscosity_eta = 0.4e-3 #Pa s
    biotcoeff_alpha = 0.31 #-
    # undrained_nu_u = 0.3  #-
    # shear_modulus_mu = 32.04e9 #Pa
    # drained_nu = 0.25 #-
  
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
  [./mu_s]
      order = CONSTANT
      family = MONOMIAL
  []
  [./mu_d]
      order = CONSTANT
      family = MONOMIAL
  []
  [./ini_shear_stress]
      order = FIRST
      family = LAGRANGE
  []
  [./tangent_jump_rate]
      order = CONSTANT
      family = MONOMIAL
  []
  [./tria_area_aux]
      order = CONSTANT
      family = MONOMIAL
  []
  [./nodal_area]
      order = FIRST
      family = LAGRANGE
  [../]
  #obtain parameters from MaterialRealAux, pass parameters to subApp
  # [./alpha_old]
  #   order = FIRST
  #   family = LAGRANGE
  # []
  [./B_old]
    order = FIRST
    family = LAGRANGE
  []
  [./xi_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./I2_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./mu_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./lambda_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [./gamma_old]
      order = CONSTANT
      family = MONOMIAL
  []
  #updated alpha, B
  [./alpha_in]
    order = CONSTANT
    family = MONOMIAL
  []
  [./B_in]
    order = CONSTANT
    family = MONOMIAL
  []
  #high-order-dummy
  [./alpha_in_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [./B_in_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #grad_alpha
  [./alpha_grad_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [./alpha_grad_y]
    order = CONSTANT
    family = MONOMIAL
  []
  #mechanical strain rate
  [./mechanical_strain_rate]
    order = CONSTANT
    family = MONOMIAL
  []
  #track Cd
  [./track_Cd]
    order = CONSTANT
    family = MONOMIAL
  []
  #pressure
  [pressure]
    order = FIRST
    family = LAGRANGE
  []
  #
  [./initial_stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [./initial_stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [./initial_stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  #
  [arctanyx]
    order = FIRST
    family = LAGRANGE
  []
  [tangent_jump]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm_ik]
    boundary = 'Block101_Block201 Block102_Block202 Block103_Block203 Block104_Block204 Block105_Block205 Block106_Block206
                Block107_Block207 Block108_Block208 Block109_Block209 Block110_Block210 Block111_Block211 Block112_Block212'
    strain = SMALL
    generate_output='tangent_jump normal_jump jump_x jump_y traction_x traction_y'
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
      [../]
    [../]
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
  [restore_x]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_x'
    variable = 'resid_x'
    execute_on = 'TIMESTEP_END'
  []
  [restore_y]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_y'
    variable = 'resid_y'
    execute_on = 'TIMESTEP_END'
  []
  [StaticFricCoeff]
    type = FunctionAux
    variable = mu_s
    function = func_static_friction_coeff_mus
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [DynamicFricCoeff]
      type = FunctionAux
      variable = mu_d
      function = func_dynamic_friction_coeff_mud
      execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
  # [StrikeShearStress]
  #   type = FunctionAux
  #   variable = ini_shear_stress
  #   function = func_initial_strike_shear_stress
  #   execute_on = 'INITIAL TIMESTEP_BEGIN'
  # []
  [TJump_rate]
    type = FDCompVarRate
    variable = tangent_jump_rate
    coupled = tangent_jump
    execute_on = 'TIMESTEP_BEGIN'
  []
  [arctanyx]
      type = CompArctanyx
      variable = arctanyx
      execute_on = 'INITIAL'
  []
  [get_xi_old]
      type = MaterialRealAux
      property = xi
      variable = xi_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_I2_old]
      type = MaterialRealAux
      property = I2
      variable = I2_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_mu_old]
      type = MaterialRealAux
      property = shear_modulus
      variable = mu_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_lambda_old]
      type = MaterialRealAux
      property = lambda
      variable = lambda_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [get_gamma_old]
      type = MaterialRealAux
      property = gamma_damaged
      variable = gamma_old
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  #define shear strain material property (elastic) inside damage stress 
  #and compute its rate using "MaterialRateRealAux"
  [get_shear_strain_rate]
      type = MaterialRateRealAux
      property = principal_strain
      variable = mechanical_strain_rate
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  #fault length
  [fault_len]
      type = ConstantAux
      variable = nodal_area
      value = 25
      execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = true
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = true
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
  #damage breakage model
  [stress_medium]
    type = ComputeDamageBreakageStressBorehole
    alpha_in = alpha_in_dummy
    B_in = B_in_dummy
    alpha_grad_x = alpha_grad_x
    alpha_grad_y = alpha_grad_y
    pressure = pressure
    output_properties = 'eps_p eps_e eps_total I1 sts_total'
    outputs = exodus
  []
  [density]
      type = GenericConstantMaterial
      prop_names = density
      prop_values = 2670
  []
  #SlipWeakeningMultifaults ONLY supports TRIA currently!
  [./czm_mat2]
    type = SlipWeakeningMultifaults
    disp_slipweakening_x     = disp_slipweakening_x
    disp_slipweakening_y     = disp_slipweakening_y
    reaction_slipweakening_x = resid_slipweakening_x
    reaction_slipweakening_y = resid_slipweakening_y
    nodal_area = nodal_area
    mu_d = mu_d
    mu_s = mu_s
    tria_area = tria_area_aux
    boundary = 'Block101_Block201 Block102_Block202 Block103_Block203 Block104_Block204 Block105_Block205 Block106_Block206
                Block107_Block207 Block108_Block208 Block109_Block209 Block110_Block210 Block111_Block211 Block112_Block212'
  [../]
  [./static_initial_stress_tensor_slipweakening]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor_slipweakening
      tensor_functions = 'func_stress_xx              func_stress_xy               func_initial_stress_00 
                          func_stress_xy              func_stress_yy               func_initial_stress_00
                          func_initial_stress_00      func_initial_stress_00       func_initial_stress_00'
  [../]
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_stress_xx              func_stress_xy             func_initial_stress_00 
                          func_stress_xy              func_stress_yy             func_initial_stress_00
                          func_initial_stress_00      func_initial_stress_00     func_initial_stress_00'
  [../]
[]

[Functions]
  #-----------slipweakening---------------#
   #mus constant value: 0.7
   [func_static_friction_coeff_mus]
    type = ConstantFunction
    value = 0.677
  []
  #mud constant value: 0.4
  [func_dynamic_friction_coeff_mud]
    type = ConstantFunction
    value = 0.4
  []
  #----------------------------------------#
  #pressure
  [func_pressure_x]
    type = InitialStressXYPressureBorehole_x
  []
  [func_pressure_y]
    type = InitialStressXYPressureBorehole_y
  []
  [func_stress_xx]
    type = SolutionFunction
    solution = load_stress_xx
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_stress_xy]
    type = SolutionFunction
    solution = load_stress_xy
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_stress_yy]
    type = SolutionFunction
    solution = load_stress_yy
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  [../]
  [func_initial_stress_00]
    type = ConstantFunction
    value = 0
  []
  [func_frank_case]
    type = InitialStressXYfrankcase
  []
[]

[Executioner]
  type = Transient
  dt = 2e-4
  end_time = 100.0
  # num_steps = 10
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 200
  [sample_snapshots]
    type = Exodus
    interval = 2000
  []
  [snapshots]
    type = Exodus
    interval = 2000
    overwrite = true
  []
  [checkpoints]
    type = Checkpoint
    interval = 2000
    num_files = 2
  []
[]

[UserObjects]
  [recompute_residual_tag]
      type = ResidualEvaluationUserObject
      vector_tag = 'restore_tag'
      force_preaux = true
      execute_on = 'TIMESTEP_END'
  []
  [load_stress_xx]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xx_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_xy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_xy_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
  [load_stress_yy]
    type = SolutionUserObject
    mesh = ./static_solve/borehole_static_solve_out.e
    system_variables = stress_yy_saved
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    timestep = LATEST
    force_preaux = true
  []
[]

[BCs]
  # [Pressure_left_x]
  #   type = Pressure
  #   boundary = left
  #   function = func_pressure_x
  #   variable = disp_x
  # []
  # [Pressure_right_x]
  #   type = Pressure
  #   boundary = right
  #   function = func_pressure_x
  #   variable = disp_x
  # []
  [Pressure_top_y]
    type = Pressure
    boundary = top
    function = func_pressure_y
    variable = disp_y
  []
  [Pressure_bottom_y]
    type = Pressure
    boundary = bottom
    function = func_pressure_y
    variable = disp_y
  []
  [./dashpot_top_x]
    type = NonReflectDashpotBC
    component = 0
    variable = disp_x
    disp_x = disp_x
    disp_y = disp_y
    p_wave_speed = 6000
    shear_wave_speed = 3464
    boundary = top
  []
  [./dashpot_top_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = top
  []
  [./dashpot_bottom_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_bottom_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = bottom
  []
  [./dashpot_left_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_left_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = left
  []
  [./dashpot_right_x]
      type = NonReflectDashpotBC
      component = 0
      variable = disp_x
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
  [./dashpot_right_y]
      type = NonReflectDashpotBC
      component = 1
      variable = disp_y
      disp_x = disp_x
      disp_y = disp_y
      p_wave_speed = 6000
      shear_wave_speed = 3464
      boundary = right
  []
[]

[MultiApps]
  [./sub_app]
      type = TransientMultiApp
      positions = '0 0 0'
      input_files = 'test_borehole_sub.i'
      execute_on = 'TIMESTEP_BEGIN'
  [../]
[]

[Transfers]
  [pull_resid]
      type = MultiAppCopyTransfer
      from_multi_app = sub_app
      source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub track_Cd alpha_checked_dummy B_checked_dummy'
      variable = 'alpha_in B_in alpha_grad_x alpha_grad_y track_Cd alpha_in_dummy B_in_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
  #we actually don't need to pass alpha and B
  [push_disp]
      type = MultiAppCopyTransfer
      to_multi_app = sub_app
      source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate alpha_in_dummy B_in_dummy'
      variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate_sub alpha_old_dummy B_old_dummy'
      execute_on = 'TIMESTEP_BEGIN'
  []
[]