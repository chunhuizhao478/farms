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
      coord = '-2000   -2000  0'
      new_boundary = corner_ptr1
      input = sidesets_boreholebc
    []
    [./extranodeset2]
        type = ExtraNodesetGenerator
        coord = '2000  -2000  0'
        new_boundary = corner_ptr2
        input = extranodeset1
    []
    [./extranodeset3]
        type = ExtraNodesetGenerator
        coord = '2000   2000  0'
        new_boundary = corner_ptr3
        input = extranodeset2
    []
    [./extranodeset4]
        type = ExtraNodesetGenerator
        coord = '-2000  2000  0'
        new_boundary = corner_ptr4
        input = extranodeset3
    []
  []

  [GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y'
  []

  [AuxVariables]
    [./initial_stress_xx]
      order = CONSTANT
      family = MONOMIAL
    [../]
    [./initial_stress_xy]
      order = CONSTANT
      family = MONOMIAL
    []
    [./initial_stress_yy]
      order = CONSTANT
      family = MONOMIAL
    [../]
[]

[AuxKernels]
    [initial_stress_xx]
        type = SolutionAux
        variable = initial_stress_xx
        solution = load_stress_xx
        execute_on = 'INITIAL'
    []
    [initial_stress_xy]
        type = SolutionAux
        variable = initial_stress_xy
        solution = load_stress_xy
        execute_on = 'INITIAL'
    []
    [initial_stress_yy]
        type = SolutionAux
        variable = initial_stress_yy
        solution = load_stress_yy
        execute_on = 'INITIAL'
    []
[]

[Modules]
    [./TensorMechanics]
      [./Master]
        [./all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          # generate_output = 'stress_xx stress_yy stress_xy strain_xx strain_xy strain_yy'
          # extra_vector_tags = 'restore_tag'
        [../]
      [../]
    [../]
[]

[UserObjects]
    [load_stress_xx]
        type = SolutionUserObject
        mesh = test_borehole_main_out.e
        system_variables = stress_xx_saved
        execute_on = 'INITIAL'
        timestep = LATEST
        force_preaux = true
    []
    [load_stress_xy]
        type = SolutionUserObject
        mesh = test_borehole_main_out.e
        system_variables = stress_xy_saved
        execute_on = 'INITIAL'
        timestep = LATEST
        force_preaux = true
    []
    [load_stress_yy]
        type = SolutionUserObject
        mesh = test_borehole_main_out.e
        system_variables = stress_yy_saved
        execute_on = 'INITIAL'
        timestep = LATEST
        force_preaux = true
    []
[]

[Executioner]
    type = Transient
    dt = 2e-4
    end_time = 40.0
    num_steps = 1
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
[]

[Outputs]
    exodus = true
    interval = 1
[]

[Materials]
    [stress_medium]
      type = ComputeLinearElasticStress
    []
    [elasticity]
      type = ComputeIsotropicElasticityTensor
      lambda = 32.04e9
      shear_modulus = 32.04e9
    []
[]