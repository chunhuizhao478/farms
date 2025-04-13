
# Verification of Benchmark Problem TPV14-2D from the SCEC Dynamic Rupture Validation exercises #
# Reference: #
# Harris, R. M.-P.-A. (2009). The SCEC/USGS Dynamic Earthquake Rupture Code Verification Exercise. Seismological Research Letters, vol. 80, no. 1, pages 119-126. #
# [Note]: This serves as a test file, to run the full problem, please extend the domain size by modifying nx, ny, xmin, xmax, ymin, ymax

[Mesh]
    [./msh]
      type = FileMeshGenerator
      file = '../../../meshgenerator/tpv142d/adaptive_mesh_v2/tpv142d_100m_adaptive.msh'
    []
    [./subdomain_id] # copy from the outputs of meshio_run_tria.py
        input = msh
        type = SubdomainPerElementGenerator 
        element_ids = '
        40 49 93 128 143 148 206 215 245 321 327 400 407 408 457 507 516 588 594 625 690 716 740 762 843 896 910 963 985 997 1091 1094 1145 1175 1217 1289 1295 1468 1493 1613 1634 1701 1706 1735 1750 1924 1974 2019 2053 2076 2129 2147 2168 2261 2295 2306 2405 2408 2424 2460 2549 2586 2590 2644 2649 2696 2945 2979 3006 3051 3075 3110 3175 3278 3285 3311 3319 3421 3484 3491 3500 3821 3823 3826 3910 3914 3949 3991 4082 4136 4143 4239 4277 4308 4351 4384 4421 4501 4663 4688 4799 4831 4860 5054 5065 5093 5100 5213 5241 5262 5274 5341 5352 5453 5481 5630 5687 5757 5767 5785 5879 5930 5934 5936 5943 5951 6013 6158 6185 6218 6252 6292 6314 6325 6371 6511 6568 6585 6737 6745 6870 6874 7213 7470 7539 7577 7586 7653 7992 8070 8161 8403 8521 8688 8718 8799 8864 8972 9042 9047 9384 9503 9530 9587 9678 9817 9855 9930 10533 10924 11129 11230 11320 11387 11529 11558 11629 11855 11976 12032 12066 12104 12213 13174 13291 13451 13534 13667 13819 13878 13895 14124 14190 14272 14661 14834 14905 15169 15196 15328 15502 15618 15674 15801 15841 16043 16048 16464 16502 16708 17034 17645 17661 17715 18003 18421 18782 18831 19122 19228 19935 20309 20343 22024 22406 22719 23336 24007 25121 26157 27988 28338 28704 28817 29882 30810 31026 31654 31859 32612 33003 33965 34403 34564 36395 37492 37579 39141 41400 41415 41466 41476 41487 41531 41609 41748 41825 41826 41828 41894 41922 41928 41933 41954 41957 41970 41997 43663 45484 45763 46581 46985 47073 47378 49504 51260 51663 51667 51950 51951 52041 52344 52361 52379 52514 52783 52842 53265 53954 54014 54653 55045 55190 55393 55735 55824 57307 57742 58033 58245 58567 59901 59908 60430 60984 61024 61198 61420 61624 62399 62497 62664 63060 63116 63602 63737 63756 63884 64141 64251 64353 64641 64754 64927 64964 65025 66184 66651 66719 67442 67465 67673 67678 67817 67894 67912 68401 68649 68890 69206 69371 69766 69802 69819 70080 70109 70424 70444 70572 70712 70860 71027 71412 71464 71586 71707 71917 71971 72121 72244 72299 72305 72346 72473 72559 72835 72877 73054 73061 73121 73173 73265 73369 73409 73518 73715 73751 73796 73940 74001 74044 74082 74224 74310 74323 74348 74426 74447 74456 74489 74522 74549 74576 74633 74705 74879 74880 75032 75039 75072 75086 75240 75307 75409 75430 75664 75683 75803 76022 76146 76236 76263 76493 76554 76707 76778 76956 77063 77133 77203 77245 77347 77507 77757 77770 77773 77901 78019 78081 78219 78228 78308 78363 78399 78421 78519 78537 78643 78694 78762 78765 78810 78874 78895 78958 78959 78987 79245 79255 79273 79588 79657 79680 79707 79713 79854 79989 80042 80059 80177 80199 80283 80337 80364 80445 80453 80513 80537 80552 80639 80659 80660 80807 80845 80861 80894 80952 80970 80980 80988 81115 81212 81234 81368 81426 81525 81540 81550 81582 81791 81883 81965 81983 81995 82014 82107 82123 82175 82252 82271 82284 82352 82405 82408 82448 82523 82599 82647 82697 82768 82886 82895 82918 82932 83128 83181 83197 83286 83424 83456 83472 83483 83492 83544 83562 83593 83608 83616 83677 83721 83748 83819 83828 83869 83884 83893 83898 83919 84029 84030 84033 84068 84119 84132 84144 84174 84186 84201 84209 84220 84249 84254 84283 84355 84377 84390 84393 84413 84426 84445 84454 84613 84769 84850 85098 85304 85346 85356 85398 85423 85477 85536 85642 85712 85730 86022 86341 7 39 74 152 160 186 243 254 326 418 492 569 613 664 674 733 745 766 817 862 893 923 1019 1037 1039 1052 1058 1067 1076 1096 1101 1119 1139 1151 1168 1184 1317 1342 1352 1366 1380 1408 1475 1506 1596 1621 1788 1791 1822 1839 1887 1939 1954 2012 2014 2027 2038 2065 2120 2122 2177 2221 2232 2302 2357 2401 2520 2545 2686 2707 2735 2774 2781 2782 2810 2830 2857 2888 2900 2978 2989 2996 3045 3046 3233 3262 3283 3318 3352 3406 3470 3499 3502 3539 3617 3620 3645 3693 3712 3725 3733 3789 3798 3804 3860 3877 3880 3887 3917 3927 3935 4008 4098 4144 4148 4165 4172 4212 4220 4258 4315 4362 4429 4474 4491 4505 4540 4557 4583 4600 4615 4635 4680 4697 4764 4792 4807 4829 4893 4903 4908 5002 5006 5026 5110 5125 5373 5435 5495 5516 5556 5586 5734 5745 5774 5842 5898 5918 5935 5959 5978 5986 5992 6095 6113 6215 6294 6311 6316 6340 6398 6549 6592 6644 6684 6730 7051 7118 7121 7156 7286 7308 7359 7376 7443 7672 7728 7754 7872 7926 8024 8064 8447 8471 8535 8589 8605 8714 8803 8837 8888 8955 8986 9088 9111 9127 9129 9209 9231 9253 9324 9474 9613 9665 9679 9910 9911 9966 9974 10001 10321 10483 10575 10590 10939 10959 11096 11132 11274 11419 11420 11433 11618 11687 11758 11807 12290 12305 12415 12426 12462 12646 12688 12738 12741 12912 13172 13248 13289 13364 13375 13629 13715 13752 13884 14018 14089 14141 14365 14418 15553 16050 16165 16298 16462 16795 17072 17122 17218 17243 17244 17960 18145 20172 21059 21295 23674 24851 25258 25306 25585 26065 26668 26863 28637 29537 30315 30709 31733 31800 32665 32666 33241 33952 34805 35234 35823 36842 38399 41210 41409 41797 41802 41824 42229 42252 42286 42418 43252 43436 43455 44149 44329 45927 47677 47689 48421 48430 48665 50027 51647 51959 52263 52291 52299 52345 52712 52794 53061 53194 53480 53489 53526 53923 53995 54022 54083 54172 54256 55230 55318 55424 55625 58535 58679 58721 59161 59990 60047 60646 60924 61690 61713 62184 62502 62537 62616 62951 63021 63079 63088 63138 63151 63233 63307 63512 63615 65191 65262 66050 66052 66208 66287 66476 66486 66586 66653 66918 66946 67542 67568 67621 67854 68302 68354 68510 68929 69006 69038 69051 69245 69383 69649 69824 69908 70025 70578 70631 70710 70717 70955 71036 71608 71644 71740 71748 71809 71861 72122 72184 72201 72269 72838 72870 72878 72915 72921 73230 73328 73509 73527 73607 74109 74120 74217 74322 74457 74507 74598 74634 74673 74738 74749 75130 75170 75370 75720 75852 76056 76089 76358 76464 76578 76627 76799 76920 76931 77001 77192 77227 77280 77364 77501 77583 77814 77839 77863 77883 78141 78301 78307 78639 78644 78665 78700 78793 78823 78997 79118 79181 79209 79223 79246 79399 79484 79520 79541 79768 79810 79849 79992 80294 80362 80428 80520 80539 80560 80567 80650 80761 80796 80806 80848 81057 81077 81319 81603 81606 81638 81677 81696 81702 81719 81882 82072 82210 82358 82369 82391 82432 82546 82753 82776 82836 82846 82881 82975 83044 83152 83198 83258 83287 83307 83412 83413 83504 83565 83581 83595 83599 83626 83795 83820 83905 83913 83922 83958 83975 83980 84002 84005 84079 84135 84139 84146 84238 84266 84299 84340 84345 84358 84385 84391 84406 84421 84477 84622 85154 85200 85259 85460 85654 85885 86002 86041 86121 86167 86197 86254 86340
        12 172 396 431 624 629 739 780 1087 1197 1488 1650 1689 1755 1967 2264 2581 2603 2747 2842 2889 2911 2983 3354 3819 4032 4432 4562 4696 4837 5144 5251 5270 5799 5938 5984 6268 6285 6355 6813 6994 7327 7403 7680 8195 8343 8556 8602 8618 8783 9686 10197 10241 10330 10926 11160 11385 11583 11728 11967 12006 12301 12419 12522 13113 14195 14313 14896 15428 15773 15990 16246 17028 17447 18063 18109 18375 22846 25416 26714 27225 27243 28814 30185 30909 35060 36544 37719 37771 37783 38159 42166 42167 42173 42183 42184 42258 42285 42289 42321 42380 42474 43085 43151 43658 44212 44435 45044 45211 45573 46653 47701 48717 49245 49604 50609 52942 54391 54708 56135 56336 57169 59186 59730 60809 61175 61337 62219 62308 63185 63380 63439 63550 63779 63835 64542 64572 65446 65678 65879 66998 67455 67810 68220 68296 68513 68645 69050 69169 69202 69361 69625 70364 70395 70408 70715 70770 71049 71262 71311 71414 71424 71457 71571 71765 71852 71853 72104 72183 72451 72478 72574 72823 73024 73066 73090 73119 73462 73858 73931 74011 74406 74438 74462 74923 75063 75160 75196 75219 75333 75376 75645 75998 76180 76400 76499 76562 76664 76708 76977 77273 77285 77456 77509 77608 78711 78768 78894 79034 79219 79260 79673 80069 80184 80361 80573 80683 80764 80869 80955 81199 81269 81273 81721 81775 82148 82185 82194 82394 82547 82558 82740 82781 83054 83324 83530 83532 83773 83863 83909 83982 84121 84143 85468 85726 85732 85840 85883 86204 89 122 373 392 542 558 754 1038 1075 1384 1409 1578 1585 2074 2078 2126 2185 2423 2501 2571 2572 2776 2923 3148 3517 3523 3562 3637 3638 3727 3889 3942 4017 4115 4169 4355 4370 4460 4493 4580 4734 4898 5038 5201 5231 5264 5404 5455 5497 5557 5724 5829 5894 6224 6308 6425 6495 6535 6656 7009 7160 8084 8148 8299 8302 9350 9486 9636 9731 10250 10658 10838 10841 11344 11431 11501 11677 12094 12634 12967 13751 13837 13871 14931 15493 17077 17377 17514 19225 20266 20328 20735 21323 22653 22715 22883 23578 24571 26794 28122 33342 34313 34363 36864 39319 46442 50784 50910 51384 51430 51493 51895 54751 56278 58149 59497 59585 61477 61896 62596 63074 63522 63643 63700 64296 64950 65304 65596 66365 67859 67922 68350 69009 69018 69039 69636 69769 69804 70354 70454 70580 70796 71008 71150 71265 71332 71640 73397 73525 73601 73606 73852 74084 74237 74340 74745 75277 75368 75930 76064 76157 76164 76609 76879 77123 77274 77596 77725 77853 77888 77904 78015 78045 78116 78173 78194 78466 78754 78860 78869 79037 79281 79417 79518 79545 79562 80147 80322 80446 80476 80607 80827 80916 81092 81271 81286 81290 81599 81658 82028 82121 82671 82682 82744 82795 82885 82913 83063 83143 83200 83215 83338 83373 83402 83433 83529 83645 83702 83754 83803 83904 84009 84041 84099 84125 84141 84297 84322 84370 84373 84400 84422 84594 84640 84701 84732 84932 84999 85031 85224 85498 85540 85671 85835 85857 86126 86241 86342 
        '
        subdomain_ids = '
        100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 
        300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 300 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400 400
        '
    []
    [./split]
      type = BreakMeshByBlockGenerator
      input = subdomain_id
      split_interface = true
      block_pairs = '100 200; 300 400'
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
    #primary variables
    displacements = 'disp_x disp_y'
    #damping ratio
    q = 0.1
    #characteristic length (m)
    Dc = 0.4
    #initial normal stress (Pa)
    T2_o = 120e6
    #dynamic friction coefficient
    mu_d = 0.525
    #element edge length (m)
    len = 100
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
    [./ini_shear_stress]
      order = FIRST
      family = MONOMIAL
    []
    [./ini_normal_stress]
      order = FIRST
      family = MONOMIAL
    []
    #global quantities
    [./global_jump_x]
      order = FIRST
      family = MONOMIAL
    []
    [./global_jump_y]
      order = FIRST
      family = MONOMIAL
    []
    [./global_traction_x]
      order = FIRST
      family = MONOMIAL
    []
    [./global_traction_y]
      order = FIRST
      family = MONOMIAL
    []
    #local quantities
    [./local_shear_jump]
        order = FIRST
        family = MONOMIAL
    []
    [./local_shear_jump_rate]
        order = FIRST
        family = MONOMIAL
    []
    [./local_normal_jump]
      order = FIRST
      family = MONOMIAL
    []
    [./local_normal_jump_rate]
      order = FIRST
      family = MONOMIAL
    []
    [./local_shear_traction]
        order = FIRST
        family = MONOMIAL
    []
    [./local_normal_traction]
      order = FIRST
      family = MONOMIAL
    []
    [./normal_x]
      order = FIRST
      family = MONOMIAL
    []
    [./normal_y]
      order = FIRST
      family = MONOMIAL
    []
    [./tangent_x]
      order = FIRST
      family = MONOMIAL
    []
    [./tangent_y]
      order = FIRST
      family = MONOMIAL
    []
    [./fractal_shear_stress_aux]
      order = FIRST
      family = MONOMIAL      
    []
  []

  [Physics/SolidMechanics/CohesiveZone]
    [./czm_ik]
      boundary = 'Block100_Block200 Block300_Block400'
      strain = SMALL
      generate_output='traction_x traction_y jump_x jump_y'
    [../]
  []

  [Physics]
    [SolidMechanics]
      [QuasiStatic]
        [all]
          strain = SMALL
          add_variables = true
          planar_formulation = PLANE_STRAIN
          generate_output = 'stress_xx stress_yy stress_xy'
          extra_vector_tags = 'restore_tag'
        []
      []
    []
  []

  [Problem]
    extra_tag_vectors = 'restore_tag'
  []

  [Functions]
    [func_static_friction_coeff_mus]
      type = ConstantFunction
      value = 0.677
    []
    #the center is (-8000,0)
    #expand 1.5e3 on both side (-9500,0) and (-6500,0)
    [func_initial_strike_shear_stress]
      type = PiecewiseConstant
      axis=x
      x = '-1000e3 -9.5e3  -6.5e3'
      y = ' 70.0e6 81.6e6 70.0e6'
    []
    [func_initial_normal_stress]
      type = ConstantFunction
      value = -120e6
    []
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
      execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
      type = CompVarRate
      variable = vel_slipweakening_y
      coupled = disp_y
      execute_on = 'TIMESTEP_END'
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
      execute_on = 'LINEAR TIMESTEP_BEGIN'
    []
    ##
    [StrikeShearStress]
      type = FunctionAux
      variable = ini_shear_stress
      function = func_initial_strike_shear_stress
      execute_on = 'LINEAR TIMESTEP_BEGIN'
    []
    [NormalStress]
      type = FunctionAux
      variable = ini_normal_stress
      function = func_initial_normal_stress
      execute_on = 'LINEAR TIMESTEP_BEGIN'
    []
    ##
    [GetShearJump]
      type = FarmsMaterialRealAux
      material_property_name = 'local_shear_jump'
      variable = local_shear_jump
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetNormalJump]
      type = FarmsMaterialRealAux
      material_property_name = 'local_normal_jump'
      variable = local_normal_jump
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetShearJumpRate]
      type = FarmsMaterialRealAux
      material_property_name = 'local_shear_jump_rate'
      variable = local_shear_jump_rate
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetNormalJumpRate]
      type = FarmsMaterialRealAux
      material_property_name = 'local_normal_jump_rate'
      variable = local_normal_jump_rate
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetShearTraction]
      type = FarmsMaterialRealAux
      material_property_name = 'local_shear_traction'
      use_fractal_shear_stress = true
      variable = local_shear_traction
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetNormalTraction]
      type = FarmsMaterialRealAux
      material_property_name = 'local_normal_traction'
      variable = local_normal_traction
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetNormalX]
      type = FarmsMaterialRealAux
      material_property_name = 'normal_x'
      variable = normal_x
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetNormalY]
      type = FarmsMaterialRealAux
      material_property_name = 'normal_y'
      variable = normal_y
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetTangentX]
      type = FarmsMaterialRealAux
      material_property_name = 'tangent_x'
      variable = tangent_x
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetTangentY]
      type = FarmsMaterialRealAux
      material_property_name = 'tangent_y'
      variable = tangent_y
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
    []
    [GetFractalShearStress]
      type = FarmsMaterialRealAux
      material_property_name = 'fractal_shear_stress'
      use_fractal_shear_stress = true
      variable = fractal_shear_stress_aux
      boundary = 'Block100_Block200 Block300_Block400'
      ini_normal_sts = ini_normal_stress
      ini_shear_sts = ini_shear_stress
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
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
        use_displaced_mesh = false
    []
    [stress]
        type = ComputeLinearElasticStress
    []
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = 2670
    []
    [./czm_mat]
        type = SlipWeakeningFrictionczm2dParametricStudy
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        reaction_slipweakening_x = resid_slipweakening_x
        reaction_slipweakening_y = resid_slipweakening_y
        mu_s = mu_s
        ini_shear_sts = ini_shear_stress
        boundary = 'Block100_Block200 Block300_Block400'
        use_fractal_shear_stress = true
        peak_shear_stress = 81.6e6
        nucl_center = '-8e3 0 0'
        nucl_radius = 1.5e3
    [../]
    #fractal shear stress
    [./fractal_shear_stress_mainfault]
        type = FractalShearStress
        boundary = 'Block100_Block200'
        csv_file = 'examples/benchmark_tpv142D/tpv142D_parametric_study/fractal_stress_generation/main_fault_fractal_stress_85.csv'
    []
    [./fractal_shear_stress_branchfault]
        type = FractalShearStress
        boundary = 'Block300_Block400'
        csv_file = 'examples/benchmark_tpv142D/tpv142D_parametric_study/fractal_stress_generation/branch_fault_fractal_stress_85.csv'
    []
  []

  [UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
  []

  [Executioner]
    type = Transient
    dt = 0.005
    end_time = 12
    # num_steps = 10
    [TimeIntegrator]
      type = CentralDifference
      solve_type = lumped
    []
  []

  [Outputs]
    exodus = false
    [csv]
      type = CSV
      execute_on = 'timestep_end'
      time_step_interval = 20
    []
  []

  [BCs]
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

  [VectorPostprocessors]
    [main_fault]
      type = SideValueSampler
      variable = 'local_shear_jump local_normal_jump local_shear_jump_rate local_normal_jump_rate local_shear_traction local_normal_traction normal_x normal_y tangent_x tangent_y fractal_shear_stress_aux' 
      boundary = 'Block100_Block200'
      sort_by = x
    []
    [branch_fault]
      type = SideValueSampler
      variable = 'local_shear_jump local_normal_jump local_shear_jump_rate local_normal_jump_rate local_shear_traction local_normal_traction normal_x normal_y tangent_x tangent_y fractal_shear_stress_aux'
      boundary = 'Block300_Block400'
      sort_by = x
    []
  []