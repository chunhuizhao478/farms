##########################################################
# Unified Parameter Choice For CBDM Complex Network Problem
# mu_d = 0.1
# For Main Fault, 
# mu = shear stress / normal stress = 70e6 / 120e6 = 0.583
# mu_s = 0.677
# S = ( mu_s - mu ) / ( mu - mu_d ) = ( 0.677 - 0.583 ) / ( 0.583 - 0.4 ) = 0.514
# Frictional Length Scale L = G Dc / ( ( mu_s - mu_d ) sigma_yy ) = 32.04e9 * 0.4 / (( 0.677 - 0.1) * 120e6) = 185m
# Use mesh size = 25m

# Diffusion Length Scale D = 5e5
# sqrt(5e5*185/3464) = 163. using 6~7, 25m mesh to resolve it

# Check CFL condition 0.1 * 25 / 6000 ~ 0.0001s
##########################################################

[Mesh]
  [./msh]
      type = FileMeshGenerator
      file =  '../../../meshgenerator/cdbm/network_oldgeo_30m/network_fine.msh'
  []
  [./subdomain_id] 
  input = msh
  type = ElementSubdomainIDGenerator 
  element_ids = '
  '
  subdomain_ids = '
  '
[]
[./split]
  input = subdomain_id
  type = BreakMeshByBlockGenerator
  surrounding_blocks = 
  '
  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 
  '
  interface_name = czm
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

  #<coefficient gives positive damage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #under slow strain rate < low strain rate threshold
  # C_d_min = 10

  #if option 2, use Cd_constant
  Cd_constant = 1e6

  #power-law correction
  #index
  m = 0.9

  #low strain rate threshold
  # mechanical_strain_rate_threshold = -1e-4

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = 10 

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  CBCBH_multiplier = 10

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = 300

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = 0.05

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = 0.03 #1e-3

  #critical point of three phases (strain invariants ratio vs damage)
  xi_1 = 0.8248

  ##Compute parameters in granular states
  #see note_mar25 for detailed setup for solving coefficients a0 a1 a2 a3
  #check struct_param.m

  #--------------------------------------------------------------------------------#
  #Note: "computeAlphaCr" needs to change every time the related parameters changed
  #--------------------------------------------------------------------------------#

  #coefficients
  # chi = 0.75
  a0 = 7.4289e9
  a1 = -2.214e10
  a2 = 2.0929e10
  a3 = -6.0672e9

  #diffusion coefficient #for structural stress coupling
  D = 0
  
[]

[Variables]
  [alpha_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_sub_dummy]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxVariables]
  [alpha_old]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_old]
    order = CONSTANT
    family = MONOMIAL
  []
  #-high-order-dummy-#
  [alpha_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_old_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #
  [xi_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [I2_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [mu_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [lambda_old]
      order = CONSTANT
      family = MONOMIAL
  []
  [gamma_old]
      order = CONSTANT
      family = MONOMIAL
  []
  #checked
  [alpha_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  [B_checked]
    order = CONSTANT
    family = MONOMIAL
  []
  #high-order-dummy
  [alpha_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  [B_checked_dummy]
    order = FIRST
    family = MONOMIAL
  []
  #grad_alpha
  [alpha_grad_x_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  [alpha_grad_y_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #mechanical strain rate sub
  [mechanical_strain_rate_sub]
    order = CONSTANT
    family = MONOMIAL
  []
  #Track Cd
  [track_Cd]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [./timederivative_alpha]
      type = TimeDerivative
      variable = alpha_sub
  []
  [./alpha_forcing_func]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub
      healing = true
  []
  [./timederivative_B]
    type = TimeDerivative
    variable = B_sub
  []
  [./B_forcing_func]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
  #-high-order-dummy-#
  [./timederivative_alpha_dummy]
    type = TimeDerivative
    variable = alpha_sub_dummy
  []
  [./alpha_forcing_func_dummy]
      type = DamageVarForcingFuncDev
      option = 2
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      variable = alpha_sub_dummy
      healing = true
  []
  [./timederivative_B_dummy]
    type = TimeDerivative
    variable = B_sub_dummy
  []
  [./B_forcing_func_dummy]
      type = BreakageVarForcingFuncDevOld
      option = 2
      variable = B_sub_dummy
      alpha_old = alpha_old
      B_old = B_old
      xi_old = xi_old
      I2_old = I2_old
      mu_old = mu_old
      gamma_old = gamma_old
      lambda_old = lambda_old
      healing = true
  []
[]

[AuxKernels]
  [check_alpha]
    type = CheckAlphaB
    coupled = alpha_sub
    variable = alpha_checked
    execute_on = 'TIMESTEP_END'
  []
  [check_B]
    type = CheckAlphaB
    coupled = B_sub
    variable = B_checked
    execute_on = 'TIMESTEP_END'
  []
  #-high-order-dummy-#
  [check_alpha_dummy]
    type = CheckAlphaB
    coupled = alpha_sub_dummy
    variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  [check_B_dummy]
    type = CheckAlphaB
    coupled = B_sub_dummy
    variable = B_checked_dummy
    execute_on = 'TIMESTEP_END'
  []
  #
  [alpha_grad_x_sub]
    type = VariableGradientComponent
    variable = alpha_grad_x_sub
    component = x
    gradient_variable = alpha_checked_dummy
    execute_on = 'TIMESTEP_END'
 []
 [alpha_grad_y_sub]
   type = VariableGradientComponent
   variable = alpha_grad_y_sub
   component = y
   gradient_variable = alpha_checked_dummy
   execute_on = 'TIMESTEP_END'
  []
[]

#by default, subApp is using the same time step as mainApp
[Executioner]
  type = Transient
  [TimeIntegrator]
    type = ExplicitSSPRungeKutta
    order = 3
  []
[]