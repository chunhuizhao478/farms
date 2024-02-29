/*
Update Breakage Variable Using Runge-Kutta own solver

- 10/5/2023 Chunhui Zhao

Include power-law correction on Cd (function of strain rate e)
if (e < 1e-4){ Cd = 10 };
else{ Cd = pow(10, log10(1+m*log10(e/1e-4)) ) }

Cb = CdCb_multiplier * Cd

*/

#include "BreakageVarForcingFuncThermalEffects.h"

registerMooseObject("farmsApp", BreakageVarForcingFuncThermalEffects);

InputParameters
BreakageVarForcingFuncThermalEffects::validParams()
{
  InputParameters params = Kernel::validParams();

  //constant parameters
  params.addParam<Real>(   "C_d_min", 10.0, "coefficient gives positive damage evolution (small strain e < 1e-4 threshold value)");
  params.addRequiredParam<Real>(  "CBCBH_multiplier", "coefficient of healing for breakage evolution");
  params.addRequiredParam<Real>(         "a0", "parameters in granular states");
  params.addRequiredParam<Real>(         "a1", "parameters in granular states");
  params.addRequiredParam<Real>(         "a2", "parameters in granular states");
  params.addRequiredParam<Real>(         "a3", "parameters in granular states");
  params.addRequiredParam<Real>(     "xi_min", "strain invariant ratio at minimum value");
  params.addRequiredParam<Real>(     "xi_max", "strain invariant ratio at maximum value");
  params.addRequiredParam<Real>(       "xi_d", "strain invariant ratio at maximum value");
  params.addRequiredParam<Real>(       "xi_0", "strain invariants ratio: onset of damage evolution");
  params.addRequiredParam<Real>(       "xi_1", "strain invariants ratio: onset of breakage healing");
  params.addRequiredParam<Real>( "beta_width", "coefficient gives width of transitional region");
  params.addRequiredParam<Real>(     "m", "Cd power-law correction index");
  params.addParam<Real>("mechanical_strain_rate_threshold", 0, "threshold value for strain rate such that Cd takes constant value Cd_min if strain rate below this value.");
  params.addRequiredParam<Real>("CdCb_multiplier", "multiplier between Cd and Cb");
  params.addParam<Real>( "scale", 1.0, "scale the Cd power-law");
  params.addRequiredParam<Real>( "temperature", "temperature value");

  //variable parameters
  params.addRequiredCoupledVar(  "alpha_old", "damage variable at previous time step");
  params.addRequiredCoupledVar(      "B_old", "breakage variable at previous time step");
  params.addRequiredCoupledVar(     "xi_old", "strain invariant ratio at previous time step");
  params.addRequiredCoupledVar(     "I2_old", "second strain invariant at previous time step");
  params.addRequiredCoupledVar(     "mu_old", "shear modulus at previous time step");
  params.addRequiredCoupledVar( "lambda_old", "first lame constant at previous time step");
  params.addRequiredCoupledVar(  "gamma_old", "damage modulus at previous time step");
  params.addCoupledVar("mechanical_strain_rate", 0.0, "strain rate");

  //add options
  params.addRequiredParam<int>( "option", "option 1 : Cd power-law; option 2 : use constant Cd");
  params.addParam<Real>( "Cd_constant", 0.0, "constant Cd value for option 2 only");

  //add healing
  params.addParam<bool>("healing", false, "if turn on healing, true = on, false = off, default is false = off");

  return params;
}

BreakageVarForcingFuncThermalEffects::BreakageVarForcingFuncThermalEffects(const InputParameters & parameters)
 : Kernel(parameters),
  _Cd_min(getParam<Real>("C_d_min")),
  _CBCBH_multiplier(getParam<Real>("CBCBH_multiplier")),
  _a0(getParam<Real>("a0")),
  _a1(getParam<Real>("a1")),
  _a2(getParam<Real>("a2")),
  _a3(getParam<Real>("a3")),
  _xi_min(getParam<Real>("xi_min")),
  _xi_max(getParam<Real>("xi_max")),
  _xi_d(getParam<Real>("xi_d")),
  _xi_0(getParam<Real>("xi_0")),
  _xi_1(getParam<Real>("xi_1")),
  _beta_width(getParam<Real>("beta_width")),
  _scale(getParam<Real>("scale")),
  _m(getParam<Real>("m")),
  _mechanical_strain_rate_threshold(getParam<Real>("mechanical_strain_rate_threshold")),
  _CdCb_multiplier(getParam<Real>("CdCb_multiplier")),
  _alpha_old(coupledValue("alpha_old")),
  _B_old(coupledValue("B_old")),
  _xi_old(coupledValue("xi_old")),
  _I2_old(coupledValue("I2_old")),
  _mu_old(coupledValue("mu_old")),
  _lambda_old(coupledValue("lambda_old")),
  _gamma_old(coupledValue("gamma_old")),
  _mechanical_strain_rate(coupledValue("mechanical_strain_rate")),
  _option(getParam<int>("option")),
  _Cd_constant(getParam<Real>("Cd_constant")),
  _healing(getParam<bool>("healing")),
  _temp(getParam<Real>("temperature"))
{
}

Real
BreakageVarForcingFuncThermalEffects::computeQpResidual()
{ 
  
 //get parameters
    Real alpha = _alpha_old[_qp];
    Real B = _B_old[_qp];
    Real I2 = _I2_old[_qp];
    Real xi = _xi_old[_qp];
    // Real mu = _mu_old[_qp];
    // Real gamma_damaged = _gamma_old[_qp];
    // Real lambda = _lambda_old[_qp];

    //Power-law correction
    //Initialize Cd
    Real Cd = 0;
    //Check options
    if ( _option == 1 ){

      //close _option == 1
      //mooseError("Option 1 is NOT available!");

      //power-law correction on coefficient Cd(function of strain rate)
      if ( _mechanical_strain_rate[_qp] < _mechanical_strain_rate_threshold ) //Cd follows power-law
      {
        Cd = _scale * pow(10, 1 + _m * log10( _mechanical_strain_rate[_qp] / _mechanical_strain_rate_threshold ) ) * _Cd_min;
      }
      else if ( _mechanical_strain_rate[_qp] < 0 ){ //Cd remains constant
        Cd = _Cd_min;
      }

    }
    else if ( _option == 2 ){

      if ( _Cd_constant == 0.0 ){
        mooseError("For option 2, need to provide nonzero Cd_constant value !");
      }
      else{
        Cd = _Cd_constant;
      }

    }
    else{
      mooseError("Please provide valid option number!");
    }

    //Compute C_B
    Real C_B = _CdCb_multiplier * Cd;

    //Compute C_BH //keep constant 
    Real C_BH = 1e4;

    //
    Real Prob = 0;
    Real alphacr = computeAlphaCr(xi);
    Prob = 1.0 / ( exp( (alphacr - alpha) / _beta_width ) + 1.0 );

    //no healing //this formulation is used in the splitstrain article
    if ( xi >= _xi_0 && xi <= _xi_max ){
        return -1.0 * C_B * Prob * (1-B) * I2 * (xi - _xi_0) * _test[_i][_qp]; //could heal if xi < xi_0
    }
    else if ( xi < _xi_0 && xi >= _xi_min ){
      
      if ( _healing == true ){
        return -1.0 * C_BH * I2 * ( xi - _xi_0 ) * _test[_i][_qp];
      }
      else{
        return 0.0;
      }
    }
    else{
      //mooseError("xi_old is OUT-OF-RANGE!.");
      return 0;
    }

    // else{ //with healing
    //   //add pre-check to the value (cause problem if encountered large values)
    //   if ( xi >= _xi_d && xi <= _xi_max ){
    //     return -1.0 * C_B * Prob * (1-B) * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) ) * _test[_i][_qp];
    //   }
    //   else if ( xi < _xi_d && xi >= _xi_min ){
    //     if ( B != 1.0 ){ //no healing until B reaches its maximum
    //       return 0.0 * _test[_i][_qp];
    //     }
    //     else{
    //       return -1.0 * C_BH * I2 * ( ( mu - gamma_damaged * xi + 0.5 * lambda * pow(xi,2) ) -  ( _a0 + _a1 * xi + _a2 * pow(xi,2) + _a3 * pow(xi,3) ) ) * _test[_i][_qp];
    //     }
    //   }
    //   else{
    //     std::cout<<"xi: "<<xi<<std::endl;
    //     mooseError("xi_old is OUT-OF-RANGE!.");
    //     return 0;
    //   }
    // }
}

Real
BreakageVarForcingFuncThermalEffects::computeQpJacobian()
{
  return 0.0;
}

/// Function: Compute alpha_cr based on the current xi
Real 
BreakageVarForcingFuncThermalEffects::computeAlphaCr(Real xi)
{
  Real alphacr;
  if ( xi < _xi_0 )
  {
    alphacr = 1.0;
  } 
  else if ( xi > _xi_0 && xi <= _xi_1 )
  { 

    if ( _temp == 20.0 ){
      // alphacr = ((xi*2.76e5-7.100521107637101e2*xi*7.5e2-7.100521107637101e2*1.4e3-7.100521107637101e+2*pow(xi,3)*1.25e2+pow(xi,3)*4.6e4+sqrt((7.100521107637101e2*3.68e2-3.19799e5)*(xi*(-1.44e3)-pow(xi,2)*2.1e3+pow(xi,3)*5.6e2+pow(xi,4)*3.0e2+pow(xi,6)*2.5e1+3.576e3)*(-3.590922148807814e-1))*5.9e1+5.152e5)*(-5.9e1/4.0))/(xi*3.837588e7-7.100521107637101e2*xi*4.416e4+7.100521107637101e2*4.048e3-7.100521107637101e2*pow(xi,2)*2.76e4+pow(xi,2)*2.3984925e7-3.517789e6);
      alphacr = ((xi*9.9453523665e18-1.137709870639449e10*xi*1.5e9-1.137709870639449e10*3.2515e9+sqrt((1.137709870639449e10*6.630234911e9-8.669919497504318e19)*(xi*-1.6722e7-pow(xi,2)*2.1e7+pow(xi,3)*6.503e6+pow(xi,4)*3.0e6+pow(xi,6)*2.5e5+3.7767369e7)*(-1.09490798927339e-1))*2.136959e6-1.137709870639449e10*pow(xi,3)*2.5e8+pow(xi,3)*1.65755872775e18+2.15582088131165e19)*(-1.0684795e6))/(xi*1.208153281977227e26-1.137709870639449e10*xi*9.2392323484785e15-1.137709870639449e10*7.49488384574351e14-1.137709870639449e10*pow(xi,2)*4.97267618325e15+pow(xi,2)*6.502439623128238e25+9.800563699173856e24);
    }
    else if ( _temp == 600.0 ){
      //alphacr = ((xi*1.137621028125e14+sqrt((1.749510772765205e8*1.21346243e8-2.266639506521505e16)*(xi*-1.3185e6-pow(xi,2)*1.3125e6+pow(xi,3)*5.1275e5+pow(xi,4)*1.875e5+pow(xi,6)*1.5625e4+2.647641e6)*(-3.024312823089825e-4))*1.01651e5-1.749510772765205e8*xi*9.375e5-1.749510772765205e8*2.56375e6+pow(xi,3)*1.896035046875e13-1.749510772765205e8*pow(xi,3)*1.5625e5+3.1110143049125e14)*(-1.2706375e+6))/(xi*2.490470157790504e21-1.749510772765205e8*xi*1.3332918449625e13-1.749510772765205e8*4.729348474682e12+pow(xi,2)*1.062487268681955e21-1.749510772765205e+8*pow(xi,2)*5.688105140625e12+8.834000812716913e20);
      alphacr = ((xi*8.779020701874656e94-9.242581079866561e63*xi*1.50667102635956e31-9.242581079866561e63*3.623202406887102e31+sqrt((9.242581079866561e63*5.826766791345719e63-5.96882581293183e127)*(xi*-1.871648259117607e63-pow(xi,2)*2.118720409559946e63+pow(xi,3)*7.278632118790692e62+pow(xi,4)*3.026743442228494e62+pow(xi,6)*2.522286201857078e61+3.991219859628105e63)*(-4.165254719705582e-149))*1.095630767443456e74+pow(xi,3)*1.463170116979109e94-9.242581079866561e63*pow(xi,3)*2.511118377265934e30+2.111155546277365e95)*-1.594352604942768e63)/(xi*1.871531268118702e159-9.242581079866561e63*xi*1.826988520658936e95-9.242581079866561e63*3.689272971941764e94+pow(xi,2)*9.079662749304941e158-9.242581079866561e63*pow(xi,2)*8.863565304527178e94+3.779218996474023e158);
    }
    else if ( _temp == 200.0 ){
      //alphacr = ((xi*2.363936726348856e76-8.676965362967033e59*pow(xi,3)*6.755399441055744e15+sqrt((8.676965362967033e59*2.07421602174062e48-1.943695954383818e108)*(xi*-4.101191936737136e32-pow(xi,2)*4.259306016766851e32+pow(xi,3)*1.594907975397775e32+pow(xi,4)*6.084722881095501e31+pow(xi,6)*5.070602400912918e30+8.388275514830981e32)*(-2.096432380634983e-135))*4.913439498690475e73-8.676965362967033e59*xi*4.053239664633446e16+pow(xi,3)*3.939894543914759e75-8.676965362967033e59*1.062422133866941e17+6.196274854033607e76)*-1.787498876616328e+62)/(xi*6.642891812989588e139-8.676965362967033e59*pow(xi,2)*3.155257421955008e79-8.676965362967033e59*xi*7.088965019511341e79+pow(xi,2)*2.956712811882994e139-8.676965362967033e59*2.153704555619364e79+2.018182671341275e139);
      alphacr = ((xi*8.6769917534725e63-1.037208690502572e48*pow(xi,3)*2.251799813685248e15+sqrt((1.037208690502572e48*1.284452803891704e47-1.488256619010376e95)*(xi*-3.673455462073407e32-pow(xi,2)*4.259306016766851e32+pow(xi,3)*1.428566013028547e32+pow(xi,4)*6.084722881095501e31+pow(xi,6)*5.070602400912918e30+7.932831564587895e32)*(-2.509558278984746e-118))*1.996187616876173e59-1.037208690502572e48*xi*1.351079888211149e16+pow(xi,3)*1.446165292245417e63-1.037208690502572e48*3.172053759722508e16+2.037176672885399e64)*(-4.088192239362401e61))/(xi*4.555870338392138e126-1.037208690502572e48*pow(xi,2)*1.953884841381781e78-1.037208690502572e48*xi*3.931983473526235e78+pow(xi,2)*2.263907275658566e126-1.037208690502572e48*6.836800257386762e77+7.921593697393311e125);
    }
    else if ( _temp == 400.0 ){
      //alphacr = ((xi*5.14920446970841e93+sqrt((4.277630321419082e63*2.924312638044248e63-1.342486278587462e127)*(xi*-4.561293264881058e62-pow(xi,2)*4.630089990897793e62+pow(xi,3)*1.773836269675967e62+pow(xi,4)*6.614414272711133e61+pow(xi,6)*5.512011893925945e60+9.23563355109902e62)*(-1.420337799758779e-146))*1.483300654147999e72-4.277630321419082e63*pow(xi,3)*2.934709284453792e29-4.277630321419082e63*xi*1.760825570672275e30+pow(xi,3)*8.582007449514017e92-4.277630321419082e63*4.722135827986369e30+1.380900148034208e94)*-5.525725536599753e+63)/(xi*4.694247859743426e158-4.277630321419082e63*pow(xi,2)*4.448408005046146e94-4.277630321419082e63*xi*1.022539191745306e95+pow(xi,2)*2.042164244219469e158-4.277630321419082e63*3.386502802513623e94+1.554667406496272e158);
      alphacr = ((xi*8.700514760536611e63-1.037039359043885e48*xi*1.351079888211149e16-1.037039359043885e48*3.184190194465384e16-1.037039359043885e48*pow(xi,3)*2.251799813685248e15+sqrt((1.037039359043885e48*6.439674542159184e47-7.450723571477409e95)*(xi*-3.687510284555424e32-pow(xi,2)*4.259306016766851e32+pow(xi,3)*1.434031777327109e32+pow(xi,4)*6.084722881095501e31+pow(xi,6)*5.070602400912918e30+7.947000529070308e32)*(-5.058549665865397e-119))*1.988391915343732e59+pow(xi,3)*1.450085793422769e63+2.050514853269164e64)*-2.036113321311982e62)/(xi*2.289551649766872e127-1.037039359043885e48*xi*1.978863841950144e79-1.037039359043885e48*3.529048975068337e78-1.037039359043885e48*pow(xi,2)*9.795908758371046e78+pow(xi,2)*1.133389704902155e127+4.083120693646676e126);
    }
    else if ( _temp == 800.0 ){
      alphacr = ((xi*3.501862902491418e61-1.616488623180232e46*xi*3.377699720527872e15-1.616488623180232e46*8.311880304874433e15-1.616488623180232e46*pow(xi,3)*5.62949953421312e14+sqrt((1.616488623180232e46*1.036759686247107e46-1.84395305794916e92)*(xi*-2.406431638528757e31-pow(xi,2)*2.662066260479282e31+pow(xi,3)*9.358345260945165e+30+pow(xi,4)*3.802951800684688e30+pow(xi,6)*3.169126500570574e29+5.071903204405169e31)*(-5.834225891325813e-113))*1.851498556120166e56+pow(xi,3)*5.836438170819029e60+8.617422417005125e61)*-7.583738085868201e59)/(xi*3.697789148842258e122-1.616488623180232e46*xi*2.079076092113488e76-1.616488623180232e46*4.760883979485085e75-1.616488623180232e46*pow(xi,2)*9.85686778922682e75+pow(xi,2)*1.753116150526449e122+8.46758095339395e121);
    }
    else{
      mooseError("please provide a valid temperature value!");
    }

  }
  else if ( xi > _xi_1 && xi <= _xi_max )
  {

    if (_temp == 20.0 ){
      // alphacr = 6.408e10/(7.100521107637101e2*1.737762711864407e8+xi*(7.100521107637101e2*1.086101694915254e8-3.996854237288136e10)-6.394966779661017e10);
      alphacr = 2.8e10/(1.137709870639449e10*6.086218780987375+xi*(1.137709870639449e10*3.275682874589545-2.171854695246844e10)-4.035306023768636e10);
    }
    else if ( _temp == 600.0 ){
      //alphacr = 1.86046e10/(xi*(1.749510772765205e8*1.830242693136319e2-2.220930745902942e10)+1.749510772765205e8*4.290088872711533e2-5.205861668396497e10);
      alphacr = 1.7016e10/(xi*(9.242581079866561e63/3.712141432200485e53-1.569651075468783e10)+9.242581079866561e63*5.552682678145896e-54-3.235418703190112e10);
    }
    else if ( _temp == 200.0 ){
      //alphacr = 2.222222222222222e10/(xi*(8.676965362967033e59*4.148400389812175e-50-2.419436512142137e10)+8.676965362967033e59*9.32027448717148e-50-5.435784948055721e10);
      alphacr = 3.0666e10/(xi*(1.037208690502572e48*4.222752322016632e-38-2.711963030077233e10)+1.037208690502572e48*8.497835691903819e-38-5.457534440738429e10);
    }
    else if ( _temp == 400.0 ){
      //alphacr = 3.181818181818182e10/(xi*(4.277630321419082e63/8.040988958382154e52-3.636757435160836e10)+4.277630321419082e63*2.858682237775646e-53-8.359680596079935e10);
      alphacr = 2.7638e10/(xi*(1.037039359043885e48*3.82071337823005e-38-2.460415067467507e10)+1.037039359043885e48*7.718193116257537e-38-4.970265172223192e10);
    }
    else if ( _temp == 800.0 ){
      alphacr = 1.094e10/(1.616488623180232e46*2.141140911968788e-36+xi*(1.616488623180232e46*1.01511161459832e-36-1.052426799056748e10)-2.219848580103605e10);
    }
  }
  else
  {
    std::cout<<"xi: "<<xi<<std::endl;
    mooseError("xi exceeds the maximum allowable range!");
  }
  return alphacr;

}