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
      file =  '../../../../../meshgenerator/cdbm/borehole_wofaults_small_refined/borehole_wofaults.msh'
  []
  [./subdomain_id] 
    input = msh
    type = ElementSubdomainIDGenerator 
    element_ids = '
    443 941 2180 3243 3421 3697 4456 5317 5801 6365 6429 7413 7447 7717 8257 9459 9794 9875 9928 10173 10301 10442 11591 12521 13423 13725 14093 14987 16159 16313 16689 16719 17941 18034 19387 19515 20287 20831 21282 21641 23329 26383 26458 27277 27908 27921 28481 29738 29868 30936 30991 31550 31709 33055 33186 33711 35422 35429 35681 35732 36391 36992 37181 38457 38487 38542 38691 39801 40450 40991 42633 42651 43013 43817 44436 44565 44979 45412 47024 47125 47450 47870 49009 49407 49901 50533 51321 51679 51861 52117 52348 52586 52723 52821 52867 53267 53537 53673 53899 54164 54215 54222 55195 55600 56946 57122 57291 57475 57773 57817 58151 58335 58761 58957 59002 59070 59104 59193 59382 59431 59535 59659 60502 61377 61404 61567 61768 62353 62602 63528 63835 64293 64462 64793 65571 65587 66547 66901 67321 67465 67691 67955 68162 68615 69058 69995 70121 70615 71063 71387 71907 72873 73070 73141 73390 73943 74938 75117 75176 75275 75285 75289 75344 75414 75430 75514 75956 76472 76520 76561 77497 78564 78885 79025 79137 79633 80767 80863 80931 81503 81650 81669 82261 82913 83137 83219 83237 83273 83421 83476 83577 83809 84191 84328 84441 84447 84695 84762 84927 84943 442 940 2181 3242 3420 3696 4457 5316 5800 6364 6428 7412 7446 7716 8256 9458 9795 9874 9929 10172 10300 10443 11590 12520 13422 13724 14092 14986 16158 16312 16688 16718 17940 18035 19386 19514 20286 20830 21283 21640 23328 26382 26459 27276 27909 27920 28480 29739 29869 30937 30990 31551 31708 33054 33187 33710 35423 35428 35680 35733 36390 36993 37180 38456 38486 38543 38690 39800 40451 40990 42632 42650 43012 43816 44437 44564 44978 45413 47025 47124 47451 47871 49008 49406 49900 50532 51320 51678 51860 52116 52349 52587 52722 52820 52866 53266 53536 53672 53898 54165 54214 54223 55194 55601 56947 57123 57290 57474 57772 57816 58150 58334 58760 58956 59003 59071 59105 59192 59383 59430 59534 59658 60503 61376 61405 61566 61769 62352 62603 63529 63834 64292 64463 64792 65570 65586 66546 66900 67320 67464 67690 67954 68163 68614 69059 69994 70120 70614 71062 71386 71906 72872 73071 73140 73391 73942 74939 75116 75177 75274 75287 75288 75345 75415 75431 75515 75957 76473 76521 76560 77496 78565 78884 79024 79136 79632 80766 80862 80930 81502 81651 81668 82260 82912 83136 83218 83236 83272 83420 83477 83576 83808 84190 84329 84440 84446 84694 84763 84926 84942 533 939 3230 3315 3671 4425 4498 5209 6241 6302 7491 7596 7981 8237 9344 9763 9832 9948 10280 10305 10419 11670 12629 13394 13709 13943 15181 15946 16028 16560 16603 17779 18142 19364 19536 20265 20763 21571 21713 23523 26249 27028 27837 27934 28397 29834 29887 30996 31265 31417 31525 33028 33149 33592 35362 35504 35711 35750 36029 37019 37167 38385 38419 38613 38625 39578 40461 41204 42558 42843 42889 43658 44493 44773 44997 45529 47117 47221 47521 47810 49058 49233 49855 50541 51206 51359 51716 52141 52331 52632 52685 52809 53273 53520 53657 53799 54141 54296 54366 55015 55632 56969 57097 57341 57451 57622 58075 58148 58371 58693 59009 59019 59068 59085 59161 59169 59533 59691 59805 60478 61318 61333 61375 61626 62437 62537 63655 63787 64274 64339 64781 65581 65599 65841 66668 66733 67319 67324 67673 67695 67819 68152 68727 69092 70074 70088 70565 70995 71241 71841 72911 73049 73059 73398 74073 74930 75083 75207 75255 75286 75291 75382 75384 75462 75476 75960 76332 76500 76574 77539 78606 78862 79077 79154 79553 80769 80848 80915 80942 81500 81620 81673 82249 82832 82908 83081 83238 83296 83429 83458 83590 83777 84225 84320 84424 84431 84693 84771 84924 84950 532 938 3231 3314 3670 4424 4499 5208 6240 6303 7490 7597 7980 8236 9345 9762 9833 9949 10281 10304 10418 11671 12628 13395 13708 13942 15180 15947 16029 16561 16602 17778 18143 19365 19537 20264 20762 21570 21712 23522 26248 27029 27836 27935 28396 29835 29886 30997 31264 31416 31524 33029 33148 33593 35363 35505 35710 35751 36028 37018 37166 38384 38418 38612 38624 39579 40460 41205 42559 42842 42888 43659 44492 44772 44996 45528 47116 47220 47520 47811 49059 49232 49854 50540 51207 51358 51717 52140 52330 52633 52684 52808 53272 53521 53656 53798 54140 54297 54367 55014 55633 56968 57096 57340 57450 57623 58074 58149 58370 58692 59008 59018 59069 59084 59160 59168 59532 59690 59804 60479 61319 61332 61374 61627 62436 62536 63654 63786 64275 64338 64780 65580 65598 65840 66669 66732 67318 67325 67672 67694 67818 68153 68726 69093 70075 70089 70564 70994 71240 71840 72910 73048 73058 73399 74072 74931 75082 75206 75254 75284 75290 75383 75385 75463 75477 75961 76333 76501 76575 77538 78607 78863 79076 79155 79552 80768 80849 80914 80943 81501 81621 81672 82248 82833 82909 83080 83239 83297 83428 83459 83591 83776 84224 84321 84425 84430 84692 84770 84925 84951 780 781 2462 2463 3116 3117 4652 4653 8100 8101 8660 8661 9004 9005 10154 10155 10658 10659 11108 11109 11578 11579 12354 12355 13688 13689 14080 14081 14496 14497 14654 14655 14870 14871 14890 14891 14906 14907 14908 14909 15002 15003 15060 15061 15264 15265 15860 15861 16134 16135 17010 17011 18506 18507 19570 19571 20906 20907 21114 21115 21592 21593 22934 22935 23230 23231 23320 23321 23508 23509 23654 23655 23910 23911 24254 24255 24288 24289 26244 26245 27014 27015 28066 28067 28962 28963 30152 30153 30544 30545 31052 31053 31388 31389 32778 32779 33180 33181 34004 34005 34514 34515 34870 34871 35652 35653 36068 36069 36104 36105 36828 36829 37038 37039 37090 37091 37160 37161 37202 37203 37348 37349 37380 37381 39164 39165 39328 39329 39348 39349 39488 39489 39668 39669 40970 40971 41176 41177 41200 41201 41650 41651 41806 41807 42280 42281 42834 42835 43220 43221 47278 47279 47522 47523 47664 47665 48040 48041 48550 48551 48654 48655 48690 48691 48966 48967 49046 49047 49084 49085 50358 50359 50498 50499 50872 50873 51906 51907 52030 52031 53514 53515 54138 54139 54772 54773 55074 55075 56294 56295 56566 56567 56784 56785 57352 57353 57784 57785 58144 58145 58388 58389 58428 58429 58518 58519 58748 58749 59120 59121 59390 59391 59456 59457 59530 59531 59910 59911 60250 60251 61078 61079 62268 62269 62346 62347 63290 63291 64228 64229 64840 64841 65194 65195 66020 66021 68194 68195 69004 69005 69086 69087 69380 69381 70450 70451 72136 72137 73292 73293 74088 74089 74200 74201 74586 74587 74884 74885 77036 77037 77452 77453 77980 77981 78248 78249 78384 78385 78662 78663 78688 78689 78818 78819 78828 78829 78956 78957 78998 78999 82536 82537 83050 83051 83998 83999 84026 84027 84940 84941 85318 85319 86602 86603 88580 88581 90818 90819 92336 92337 92678 92679 93028 93029 93030 93031 93092 93093 93428 93429 94666 94667 
    '
    subdomain_ids = '
    1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1001 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1002 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1003 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 1004 
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
    coord = '0  -1000  0'
    new_boundary = corner_ptr1
    input = sidesets_boreholebc
  []
  [./extranodeset2]
      type = ExtraNodesetGenerator
      coord = '1000  0  0'
      new_boundary = corner_ptr2
      input = extranodeset1
  []
  [./extranodeset3]
      type = ExtraNodesetGenerator
      coord = '0   1000  0'
      new_boundary = corner_ptr3
      input = extranodeset2
  []
  [./extranodeset4]
      type = ExtraNodesetGenerator
      coord = '-1000  0  0'
      new_boundary = corner_ptr4
      input = extranodeset3
  []
  use_displaced_mesh = true
  displacements = 'disp_x disp_y'
[]

[GlobalParams]
  ##------------slip weakening------------##
  displacements = 'disp_x disp_y'
  
  #damping ratio
  q = 0.2
  
  #characteristic length (m)
  # Dc = 0.4

  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = 3.204e10
    
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = 3.204e10

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
  gamma_damaged_r = 3.7150e10

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.8248

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  # #coefficients
  # chi = 0.75
  a0 = 7.4289e9
  a1 = -2.214e10
  a2 = 2.0929e10
  a3 = -6.0672e9

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
  #arctanyx
  [./arctanyx]
    order = FIRST
    family = LAGRANGE    
  []
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

[AuxKernels]
  [Arctanyx]
    type = CompArctanyx
    variable = arctanyx
    execute_on = 'INITIAL'
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
  # #fault length
  # [fault_len]
  #     type = ConstantAux
  #     variable = nodal_area
  #     value = 10
  #     execute_on = 'INITIAL TIMESTEP_BEGIN'
  # []
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
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_stress_xx              func_stress_xy             func_initial_stress_00 
                          func_stress_xy              func_stress_yy             func_initial_stress_00
                          func_initial_stress_00      func_initial_stress_00     func_initial_stress_00'
  [../]
[]

[Functions]
  #pressure
  [func_pressure_x]
    type = InitialStressXYPressureBorehole_x_fast
  []
  [func_pressure_y]
    type = InitialStressXYPressureBorehole_y_fast
  []
  [func_stress_xx]
    type = ConstantFunction
    value = -200e6
  [../]
  [func_stress_xy]
    type = ConstantFunction
    value = 0
  [../]
  [func_stress_yy]
    type = ConstantFunction
    value = -40e6
  [../]
  [func_initial_stress_00]
    type = ConstantFunction
    value = 0
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

[BCs]
  [Pressure_borehole_x]
    type = Pressure
    boundary = borehole
    function = func_pressure_x
    variable = disp_x
  []
  [Pressure_borehole_y]
    type = Pressure
    boundary = borehole
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