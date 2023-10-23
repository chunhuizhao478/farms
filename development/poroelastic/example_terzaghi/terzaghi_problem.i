# Terzaghi problem 

[Mesh]
    [./msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 80
    xmin = 0
    xmax = 0.1
    ymin = 0
    ymax = 1
    []
[]
[GlobalParams]
    displacements = 'disp_x disp_y'
    pressure = 'p'
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [p]
  []
[]


[Executioner]
   type = Transient
   solve_type = Newton
   start_time = 0
   dt = 0.1
   end_time = 1
[]
  
[Outputs]
    exodus = true
    interval = 10
[]