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
      file =  '../../../meshgenerator/cdbm/borehole_wofaults/borehole_wofaults.msh'
  []
  [./subdomain_id] 
    input = msh
    type = ElementSubdomainIDGenerator 
    element_ids = '
    241 717 1229 1239 1257 1415 1609 1893 2383 3539 4285 5153 6673 7863 7901 8329 9015 9235 10369 10955 10979 12583 13191 14303 14783 15455 15681 16909 17721 17741 18011 18147 18697 18873 19429 19789 20091 21107 21567 21589 22561 23441 24675 25403 25415 25633 25889 26351 27063 27541 27905 29003 29439 29535 29663 29859 30037 30197 31123 31883 31969 32049 32089 32439 32503 32547 33043 33299 33413 34487 34707 34955 35365 35483 35635 36823 37313 37409 37491 37605 37785 38183 38193 38397 38605 38775 38943 38999 39611 39857 40513 40615 40823 41149 41201 42299 42321 42493 42535 42585 43013 43033 43427 43441 43489 43619 43761 44323 44475 44827 45509 46043 46117 46703 47511 47785 48081 48571 48655 48967 49639 49811 50577 50667 50997 51385 51417 51875 52587 52771 52849 53237 53503 53639 53643 54211 54241 54295 54355 54455 54613 55003 55075 55195 55585 56025 56079 56231 56309 56575 56667 56691 56692 56765 56773 56785 56799 56815 56845 56853 240 716 1228 1238 1256 1414 1608 1892 2382 3538 4284 5152 6672 7862 7900 8328 9014 9234 10368 10954 10978 12582 13190 14302 14782 15454 15680 16908 17720 17740 18010 18146 18696 18872 19428 19788 20090 21106 21566 21588 22560 23440 24674 25402 25414 25632 25888 26350 27062 27540 27904 29002 29438 29534 29662 29858 30036 30196 31122 31882 31968 32048 32088 32438 32502 32546 33042 33298 33412 34486 34706 34954 35364 35482 35634 36822 37312 37408 37490 37604 37784 38182 38192 38396 38604 38774 38942 38998 39610 39856 40512 40614 40822 41148 41200 42298 42320 42492 42534 42584 43012 43032 43426 43440 43488 43618 43760 44322 44474 44826 45508 46042 46116 46702 47510 47784 48080 48570 48654 48966 49638 49810 50576 50666 50996 51384 51416 51874 52586 52770 52848 53236 53502 53638 53642 54210 54240 54294 54354 54454 54612 55002 55074 55194 55584 56024 56078 56230 56308 56574 56666 56690 56693 56764 56772 56784 56798 56814 56844 56852 243 719 1227 1241 1255 1413 1611 1895 2381 3541 4287 5151 6675 7865 7903 8331 9013 9237 10367 10953 10981 12581 13193 14305 14785 15457 15679 16911 17723 17739 18009 18149 18695 18875 19431 19791 20093 21105 21569 21587 22559 23439 24677 25401 25413 25631 25891 26349 27065 27539 27903 29001 29437 29537 29661 29861 30039 30199 31121 31881 31971 32047 32091 32437 32505 32549 33041 33297 33415 34485 34705 34953 35367 35485 35633 36825 37311 37411 37489 37607 37783 38181 38191 38399 38607 38773 38941 38997 39613 39859 40511 40613 40825 41151 41203 42301 42319 42495 42533 42583 43015 43035 43429 43443 43487 43621 43763 44325 44477 44825 45507 46045 46119 46701 47509 47787 48079 48573 48653 48965 49641 49809 50579 50669 50995 51383 51419 51877 52585 52769 52847 53235 53505 53637 53641 54209 54243 54293 54353 54457 54615 55001 55073 55193 55587 56023 56081 56229 56307 56577 56669 56689 56695 56763 56771 56783 56797 56813 56847 56855 242 718 1226 1240 1254 1412 1610 1894 2380 3540 4286 5150 6674 7864 7902 8330 9012 9236 10366 10952 10980 12580 13192 14304 14784 15456 15678 16910 17722 17738 18008 18148 18694 18874 19430 19790 20092 21104 21568 21586 22558 23438 24676 25400 25412 25630 25890 26348 27064 27538 27902 29000 29436 29536 29660 29860 30038 30198 31120 31880 31970 32046 32090 32436 32504 32548 33040 33296 33414 34484 34704 34952 35366 35484 35632 36824 37310 37410 37488 37606 37782 38180 38190 38398 38606 38772 38940 38996 39612 39858 40510 40612 40824 41150 41202 42300 42318 42494 42532 42582 43014 43034 43428 43442 43486 43620 43762 44324 44476 44824 45506 46044 46118 46700 47508 47786 48078 48572 48652 48964 49640 49808 50578 50668 50994 51382 51418 51876 52584 52768 52846 53234 53504 53636 53640 54208 54242 54292 54352 54456 54614 55000 55072 55192 55586 56022 56080 56228 56306 56576 56668 56688 56694 56762 56770 56782 56796 56812 56846 56854 3000 3001 4298 4299 4548 4549 7584 7585 8736 8737 11308 11309 12830 12831 13088 13089 13178 13179 13670 13671 14652 14653 14734 14735 14736 14737 14824 14825 16372 16373 16414 16415 17310 17311 18272 18273 18310 18311 19702 19703 20300 20301 20492 20493 21842 21843 21916 21917 22210 22211 22488 22489 22578 22579 23852 23853 24054 24055 24494 24495 25652 25653 26862 26863 27422 27423 28322 28323 28634 28635 29500 29501 31790 31791 31808 31809 31988 31989 33078 33079 34878 34879 35092 35093 35536 35537 36384 36385 36464 36465 36600 36601 37664 37665 37728 37729 40202 40203 42412 42413 43864 43865 43916 43917 46626 46627 46680 46681 46688 46689 46744 46745 49768 49769 49816 49817 50682 50683 50734 50735 50992 50993 52170 52171 54266 54267 54340 54341 55872 55873 55894 55895 56004 56005 56296 56297 56686 56687 56758 56759 56818 56819 
    '
    subdomain_ids = '
    1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 
    '
  []
  [./sidesets_leftbc]
    input = subdomain_id
    type = SideSetsAroundSubdomainGenerator
    block = 1000
    new_boundary = 'left'
    external_only = true
  []
  [./sidesets_rightbc]
    input = sidesets_leftbc
    type = SideSetsAroundSubdomainGenerator
    block = 1001
    new_boundary = 'right'
    external_only = true
  []
  [./sidesets_topbc]
    input = sidesets_rightbc
    type = SideSetsAroundSubdomainGenerator
    block = 1002
    new_boundary = 'top'
    external_only = true
  []
  [./sidesets_bottombc]
    input = sidesets_topbc
    type = SideSetsAroundSubdomainGenerator
    block = 1003
    new_boundary = 'bottom'
    external_only = true
  []
  [./sidesets_boreholebc]
    input = sidesets_bottombc
    type = SideSetsAroundSubdomainGenerator
    block = 1004
    new_boundary = 'borehole'
    external_only = true
  []
  [./extranodeset1]
    type = ExtraNodesetGenerator
    coord = '0  -2000  0'
    new_boundary = corner_ptr1
    input = sidesets_boreholebc
  []
  [./extranodeset2]
      type = ExtraNodesetGenerator
      coord = '2000  0  0'
      new_boundary = corner_ptr2
      input = extranodeset1
  []
  [./extranodeset3]
      type = ExtraNodesetGenerator
      coord = '0   2000  0'
      new_boundary = corner_ptr3
      input = extranodeset2
  []
  [./extranodeset4]
      type = ExtraNodesetGenerator
      coord = '-2000  0  0'
      new_boundary = corner_ptr4
      input = extranodeset3
  []
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  # #damping ratio
  # q = 0.2
  
  # #characteristic length (m)
  # # Dc = 0.4

  # ##----continuum damage breakage model----##
  # #initial lambda value (first lame constant) [Pa]
  # lambda_o = 3.204e10
    
  # #initial shear modulus value (second lame constant) [Pa]
  # shear_modulus_o = 3.204e10

  # #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  # xi_0 = -1.0

  # #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  # xi_d = -1.1

  # #<strain invariants ratio: maximum allowable value>: set boundary
  # #Xu_etal_P15-2D
  # #may need a bit space, use 1.5 as boundary
  # xi_max = 1.8

  # #<strain invariants ratio: minimum allowable value>: set boundary
  # #Xu_etal_P15-2D
  # xi_min = -1.8

  # #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  # C_g = 1e-10

  # #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  # m1 = 10

  # #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  # m2 = 1

  # ##Compute gamma_damaged_r, xi_1
  # #Determine two parameters using convexity of Hessian matrix, positivity of eigenvalues
  # #two equations [15a] = 0 [15b] = 0 solves gamma_damaged_r and xi_1 
  # #check struct_param.m 

  # #coefficient of damage solid modulus
  # gamma_damaged_r = 3.49061e10

  # #critical point of three phases (strain invariants ratio vs damage)
  # xi_1 = 0.73205

  # ##Compute parameters in granular states
  # #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  # #check struct_param.m

  # #--------------------------------------------------------------------------------#
  # #Note: "computeAlphaCr" needs to change every time the related parameters changed
  # #--------------------------------------------------------------------------------#

  # # #coefficients
  # # chi = 0.75
  # a0 = 3.92486e9
  # a1 = -1.51257e10
  # a2 = 1.93525e10
  # a3 = -8.21570e9

  # #diffusion coefficient #for structural stress coupling
  # D = 0

  # #pressure analytical solution
  # #Reference: Injection-induced seismicity: Poroelastic and earthquake nucleation effects (P. Segall1 and S. Lu2)
  # # effec_sts_coeff = 0.31 #
  # # flux_q = 5e2 #kg/s
  # # density_rho_0 = 1e3 #kg/m^3
  # # permeability_k = 3e-12 #m^2
  # # viscosity_eta = 0.4e-3 #Pa s
  # biotcoeff_alpha = 0.31 #-
  # # undrained_nu_u = 0.3  #-
  # # shear_modulus_mu = 32.04e9 #Pa
  # # drained_nu = 0.25 #-
  
[]

[AuxVariables]
  [./stress_xx_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy_saved]
    order = CONSTANT
    family = MONOMIAL
  []
  [./stress_yy_saved]
    order = CONSTANT
    family = MONOMIAL
  [../]
#   [./resid_x]
#     order = FIRST
#     family = LAGRANGE
#   [../]
#   [./resid_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./resid_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   [./resid_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   [./disp_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./disp_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./vel_slipweakening_x]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./vel_slipweakening_y]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./mu_s]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./mu_d]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./ini_shear_stress]
#       order = FIRST
#       family = LAGRANGE
#   []
#   [./tangent_jump_rate]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./tria_area_aux]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./nodal_area]
#       order = FIRST
#       family = LAGRANGE
#   [../]
#   #obtain parameters from MaterialRealAux, pass parameters to subApp
#   # [./alpha_old]
#   #   order = FIRST
#   #   family = LAGRANGE
#   # []
#   [./B_old]
#     order = FIRST
#     family = LAGRANGE
#   []
#   [./xi_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./I2_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./mu_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./lambda_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   [./gamma_old]
#       order = CONSTANT
#       family = MONOMIAL
#   []
#   #updated alpha, B
#   [./alpha_in]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./B_in]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #high-order-dummy
#   [./alpha_in_dummy]
#     order = FIRST
#     family = MONOMIAL
#   []
#   [./B_in_dummy]
#     order = FIRST
#     family = MONOMIAL
#   []
#   #grad_alpha
#   [./alpha_grad_x]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./alpha_grad_y]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #mechanical strain rate
#   [./mechanical_strain_rate]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #track Cd
#   [./track_Cd]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   #pressure
#   [pressure]
#     order = FIRST
#     family = LAGRANGE
#   []
# []

# [Modules/TensorMechanics/CohesiveZoneMaster]
#   [./czm_ik]
#     boundary = 'czm'
#     strain = SMALL
#     generate_output='tangent_jump normal_jump'
#   [../]
[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = SMALL
        add_variables = true
        planar_formulation = PLANE_STRAIN
        generate_output = 'strain_xx strain_xy strain_yy'
        # extra_vector_tags = 'restore_tag'
      [../]
    [../]
  [../]
[]

# [Problem]
#   extra_tag_vectors = 'restore_tag'
# []

[AuxKernels]
  [get_stress_xx]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 0
    variable = stress_xx_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_xy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 0
    j = 1
    variable = stress_xy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
  [get_stress_yy]
    type = MaterialRankTwoTensorAux
    property = stress
    i = 1
    j = 1
    variable = stress_yy_saved
    execute_on = 'LINEAR TIMESTEP_END'
  []
#   [Displacment_x]
#     type = ProjectionAux
#     variable = disp_slipweakening_x
#     v = disp_x
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Displacement_y]
#     type = ProjectionAux
#     variable = disp_slipweakening_y
#     v = disp_y
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Vel_x]
#       type = CompVarRate
#       variable = vel_slipweakening_x
#       coupled = disp_x
#       execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Vel_y]
#       type = CompVarRate
#       variable = vel_slipweakening_y
#       coupled = disp_y
#       execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Residual_x]
#     type = ProjectionAux
#     variable = resid_slipweakening_x
#     v = resid_x
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [Residual_y]
#     type = ProjectionAux
#     variable = resid_slipweakening_y
#     v = resid_y
#     execute_on = 'TIMESTEP_BEGIN'
#   []
#   [restore_x]
#     type = TagVectorAux
#     vector_tag = 'restore_tag'
#     v = 'disp_x'
#     variable = 'resid_x'
#     execute_on = 'TIMESTEP_END'
#   []
#   [restore_y]
#     type = TagVectorAux
#     vector_tag = 'restore_tag'
#     v = 'disp_y'
#     variable = 'resid_y'
#     execute_on = 'TIMESTEP_END'
#   []
#   [StaticFricCoeff]
#     type = FunctionAux
#     variable = mu_s
#     function = func_static_friction_coeff_mus
#     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   [DynamicFricCoeff]
#       type = FunctionAux
#       variable = mu_d
#       function = func_dynamic_friction_coeff_mud
#       execute_on = 'INITIAL TIMESTEP_BEGIN'
#     []
#   [StrikeShearStress]
#     type = FunctionAux
#     variable = ini_shear_stress
#     function = func_initial_strike_shear_stress
#     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   # [TJump_rate]
#   #   type = FDCompVarRate
#   #   variable = tangent_jump_rate
#   #   coupled = tangent_jump
#   #   execute_on = 'TIMESTEP_BEGIN'
#   # []
#   #obtain parameters from MaterialRealAux
#   # [get_alpha_old]
#   #     type = MaterialRealAux
#   #     property = alpha_damagedvar
#   #     variable = alpha_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_B_old]
#   #     type = MaterialRealAux
#   #     property = B
#   #     variable = B_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_xi_old]
#   #     type = MaterialRealAux
#   #     property = xi
#   #     variable = xi_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_I2_old]
#   #     type = MaterialRealAux
#   #     property = I2
#   #     variable = I2_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_mu_old]
#   #     type = MaterialRealAux
#   #     property = shear_modulus
#   #     variable = mu_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_lambda_old]
#   #     type = MaterialRealAux
#   #     property = lambda
#   #     variable = lambda_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # [get_gamma_old]
#   #     type = MaterialRealAux
#   #     property = gamma_damaged
#   #     variable = gamma_old
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   # #define shear strain material property (elastic) inside damage stress 
#   # #and compute its rate using "MaterialRateRealAux"
#   # [get_shear_strain_rate]
#   #     type = MaterialRateRealAux
#   #     property = principal_strain
#   #     variable = mechanical_strain_rate
#   #     execute_on = 'INITIAL TIMESTEP_BEGIN'
#   # []
#   #fault length
#   [fault_len]
#       type = ConstantAux
#       variable = nodal_area
#       value = 25
#       execute_on = 'INITIAL TIMESTEP_BEGIN'
#   []
#   #initial is zero, DO NOT add initial on execute_on
#   [./aux_func_pressure]
#     type = ConstantAux
#     variable = pressure
#     value = 0
#     execute_on = 'TIMESTEP_BEGIN'
#   [../]
[]

# [Kernels]
#   [./inertia_x]
#     type = InertialForce
#     use_displaced_mesh = false
#     variable = disp_x
#   []
#   [./inertia_y]
#     type = InertialForce
#     use_displaced_mesh = false
#     variable = disp_y
#   []
#   [./Reactionx]
#     type = StiffPropDamping
#     variable = 'disp_x'
#     component = '0'
#   []
#   [./Reactiony]
#     type = StiffPropDamping
#     variable = 'disp_y'
#     component = '1'
#   []
# []

[Materials]
  # #damage breakage model
  # [stress_medium]
  #   type = ComputeDamageBreakageStressv3pressurev2
  #   alpha_in = alpha_in_dummy
  #   B_in = B_in_dummy
  #   alpha_grad_x = alpha_grad_x
  #   alpha_grad_y = alpha_grad_y
  #   pressure = pressure
  #   output_properties = 'eps_p eps_e eps_total I1'
  #   outputs = exodus
  # []
  [stress_medium]
      type = ComputeLinearElasticStress
  []
  [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
  []
  # [density]
  #     type = GenericConstantMaterial
  #     prop_names = density
  #     prop_values = 2670
  # []
  #SlipWeakeningMultifaults ONLY supports TRIA currently!
  # [./czm_mat]
  #     type = SlipWeakeningMultifaultsPressure
  #     disp_slipweakening_x     = disp_slipweakening_x
  #     disp_slipweakening_y     = disp_slipweakening_y
  #     reaction_slipweakening_x = resid_slipweakening_x
  #     reaction_slipweakening_y = resid_slipweakening_y
  #     nodal_area = nodal_area
  #     mu_d = mu_d
  #     mu_s = mu_s
  #     tria_area = tria_area_aux
  #     pressure = pressure
  #     boundary = 'czm'
  # [../]
  # [./static_initial_stress_tensor_slipweakening]
  #     type = GenericFunctionRankTwoTensor
  #     tensor_name = static_initial_stress_tensor_slipweakening
  #     tensor_functions = 'func_initial_stress_xx         func_initial_strike_shear_stress           func_initial_stress_00 
  #     func_initial_strike_shear_stress         func_initial_stress_yy           func_initial_stress_00
  #                         func_initial_stress_00         func_initial_stress_00           func_initial_stress_00'
  # [../]
  # [./static_initial_stress_tensor]
  #     type = GenericFunctionRankTwoTensor
  #     tensor_name = static_initial_stress_tensor
  #     tensor_functions = 'func_initial_stress_xx             func_initial_stress_xy_const        func_initial_stress_00 
  #                         func_initial_stress_xy_const       func_initial_stress_yy              func_initial_stress_00
  #                         func_initial_stress_00             func_initial_stress_00              func_initial_stress_00'
  # [../]
[]

[Functions]
  # #mus constant value: 0.7
  # [func_static_friction_coeff_mus]
  #     type = ConstantFunction
  #     value = 0.677
  # []
  # #mud constant value: 0.4
  # [func_dynamic_friction_coeff_mud]
  #     type = ConstantFunction
  #     value = 0.4
  # []
  # #Note:restrict stress variation along the fault only
  # #this function is used in czm only
  # [func_initial_strike_shear_stress]
  #   # type = InitialStressXYnetwork_old
  #   type = ConstantFunction
  #   value = 0.0
  # []
  # #this function is used in medimum
  # [func_initial_stress_xy_const]
  #   type = ConstantFunction
  #   value = 0.0
  # []
  # [./func_initial_stress_00]
  #   type = ConstantFunction
  #   value = 0.0
  # []
  [./func_initial_stress_neg]
    type = ConstantFunction
    value = -100e6
  []
  #In problems with inelasticity, the sigma11 is important
  #This is different from pure tpv205 
  [./func_initial_stress_pos]
    type = ConstantFunction
    value = 100e6
  []
  # #pressure
  # [func_pressure_x]
  #   type = InitialStressXYPressureBorehole_x_fast
  # []
  # [func_pressure_y]
  #   type = InitialStressXYPressureBorehole_y_fast
  # []
[]

# [UserObjects]
#   [recompute_residual_tag]
#       type = ResidualEvaluationUserObject
#       vector_tag = 'restore_tag'
#       force_preaux = true
#       execute_on = 'TIMESTEP_END'
#   []
# []

[Preconditioning]
  [smp]
    full = true
    type = SMP
  []
[]

[Executioner]
  type = Steady
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       NONZERO               1e-15'
  solve_type = NEWTON
[]

#for cluster run
[Outputs]
  exodus = true
  interval = 1
  # [sample_snapshots]
  #   type = Exodus
  #   interval = 2000
  # []
  # [snapshots]
  #   type = Exodus
  #   interval = 2000
  #   overwrite = true
  # []
  # [checkpoints]
  #   type = Checkpoint
  #   interval = 2000
  #   num_files = 2
  # []
[]

[BCs]
  [Pressure_top]
    type = Pressure
    boundary = top
    variable = disp_y
    factor = 100e6
  []
  [Pressure_bottom]
    type = Pressure
    boundary = bottom
    variable = disp_y
    factor = 100e6
  []
  [Pressure_left]
    type = Pressure
    boundary = left
    variable = disp_x
    factor = 100e6
  []
  [Pressure_right]
    type = Pressure
    boundary = right
    variable = disp_x
    factor = 100e6
  []
  [./fix_cptr1_x]
    type = DirichletBC
    variable = disp_x
    boundary = corner_ptr1
    value = 0
  []
  # [./fix_cptr1_y]
  #     type = DirichletBC
  #     variable = disp_y
  #     boundary = corner_ptr1
  #     value = 0
  # []
  # [./fix_cptr2_x]
  #     type = DirichletBC
  #     variable = disp_x
  #     boundary = corner_ptr2
  #     value = 0
  # []
  [./fix_cptr2_y]
      type = DirichletBC
      variable = disp_y
      boundary = corner_ptr2
      value = 0
  []
  [./fix_cptr3_x]
      type = DirichletBC
      variable = disp_x
      boundary = corner_ptr3
      value = 0
  []
  # [./fix_cptr3_y]
  #     type = DirichletBC
  #     variable = disp_y
  #     boundary = corner_ptr3
  #     value = 0
  # []
  # [./fix_cptr4_x]
  #     type = DirichletBC
  #     variable = disp_x
  #     boundary = corner_ptr4
  #     value = 0
  # []
  [./fix_cptr4_y]
      type = DirichletBC
      variable = disp_y
      boundary = corner_ptr4
      value = 0
  []
  # [fix_top_x]
  #   type = DirichletBC
  #   boundary = top
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_left_y]
  #   type = DirichletBC
  #   boundary = left
  #   variable = disp_y
  #   value = 0
  # []
  # [fix_bottom_x]
  #   type = DirichletBC
  #   boundary = bottom
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_bottom_y]
  #   type = DirichletBC
  #   boundary = bottom
  #   variable = disp_y
  #   value = 0
  # []
  # [fix_right_x]
  #   type = DirichletBC
  #   boundary = right
  #   variable = disp_x
  #   value = 0
  # []
  # [fix_right_y]
  #   type = DirichletBC
  #   boundary = right
  #   variable = disp_y
  #   value = 0
  # []
[]

# [MultiApps]
#   [./sub_app]
#       type = TransientMultiApp
#       positions = '0 0 0'
#       input_files = 'test_borehole_sub.i'
#       execute_on = 'TIMESTEP_BEGIN'
#   [../]
# []

# [Transfers]
#   [pull_resid]
#       type = MultiAppCopyTransfer
#       from_multi_app = sub_app
#       source_variable = 'alpha_checked B_checked alpha_grad_x_sub alpha_grad_y_sub track_Cd alpha_checked_dummy B_checked_dummy'
#       variable = 'alpha_in B_in alpha_grad_x alpha_grad_y track_Cd alpha_in_dummy B_in_dummy'
#       execute_on = 'TIMESTEP_BEGIN'
#   []
#   #we actually don't need to pass alpha and B
#   [push_disp]
#       type = MultiAppCopyTransfer
#       to_multi_app = sub_app
#       source_variable = 'alpha_in B_in xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate alpha_in_dummy B_in_dummy'
#       variable = 'alpha_old B_old xi_old I2_old mu_old lambda_old gamma_old mechanical_strain_rate_sub alpha_old_dummy B_old_dummy'
#       execute_on = 'TIMESTEP_BEGIN'
#   []
# []