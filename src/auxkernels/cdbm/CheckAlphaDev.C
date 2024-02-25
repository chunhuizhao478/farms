#include "CheckAlphaDev.h"

//Here we set threshold on maximum damage variable which governed by the shear modulus equation
//factor * mu_init = mu_init + alpha_max * xi_o * gamma
//alpha_max = (1-factor) * mu_init / ( xi_o * gamma )

registerMooseObject("farmsApp", CheckAlphaDev);

InputParameters
CheckAlphaDev::validParams()
{
  InputParameters params = AuxKernel::validParams();

  //constant parameters
  params.addRequiredParam<Real>(   "shear_modulus_o", "initial shear modulus value");
  params.addRequiredParam<Real>(   "xi_0", "strain invariant ratio at onset of damage");
  params.addRequiredParam<Real>(   "gamma_damaged_r", "damage modulus at maximum alpha value");   
  params.addRequiredParam<Real>(   "factor", "the ratio of residual shear modulus");
  params.addRequiredCoupledVar("coupled","Either Alpha or B that needs to be checked");

  return params;
}

CheckAlphaDev::CheckAlphaDev(const InputParameters & parameters)
  : AuxKernel(parameters),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _xi_0(getParam<Real>("xi_0")),
  _gamma_damaged_r(getParam<Real>("gamma_damaged_r")),
  _factor(getParam<Real>("factor")),
  _coupled_val(coupledValue("coupled"))

{
}

Real
CheckAlphaDev::computeValue()
{

  Real alpha_maxallowed = ( _factor - 1 ) * _shear_modulus_o / ( _xi_0 * _gamma_damaged_r );

  // std::cout<<alpha_maxallowed<<std::endl;
  
  Real var = _coupled_val[_qp];
  Real var_out;
  if ( var < 0 )
  {
    var_out = 0;
  }
  else if ( var > alpha_maxallowed )
  {
    var_out = alpha_maxallowed;
  } 
  else{
    var_out = var;
  }
  return var_out;
}