#include "InitialStressTPV243D.h"

registerMooseObject("farmsApp", InitialStressTPV243D);

InputParameters
InitialStressTPV243D::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("i", "index");
  params.addRequiredParam<Real>("j", "index");
  return params;
}

InitialStressTPV243D::InitialStressTPV243D(const InputParameters & parameters)
  : Function(parameters),
  _i(getParam<Real>("i")),
  _j(getParam<Real>("j"))
{
}

Real
InitialStressTPV243D::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real z_coord = p(2); //along the dip direction
  Real To = 0;
  Real fluid_density = 1000; //kg/m^3 fluid density
  Real rock_density = 2670; //kg/m^3 rock density
  Real gravity = 9.8; //m/s^2
  Real b22 = 0.926793; //0.926793; modify the coefficient to make xi larger #MODIFIED TPV24
  Real b33 = 1.073206;
  Real b23 = -0.169029; //-0.169029; modify May20
  Real sigma11 = 0;
  Real sigma22 = 0;
  Real sigma33 = 0;
  Real sigma23 = 0;

  //Pf
  Real Pf = fluid_density * gravity * abs(z_coord);

  //sigma11
  sigma11 = -1 * rock_density * gravity * abs(z_coord);

  //sigma22
  if ( abs(z_coord) <= 15600 ) {
    sigma22 = b22 * ( sigma11 + Pf ) - Pf;
  }
  else{
    sigma22 = sigma11;
  }

  //sigma33
  if ( abs(z_coord) <= 15600 ) {
    sigma33 = b33 * ( sigma11 + Pf ) - Pf;
  }
  else{
    sigma33 = sigma11;
  } 

  //sigma23
  if ( abs(z_coord) <= 15600 ) {
    sigma23 = b23 * ( sigma11 + Pf );
  }
  else{
    sigma23 = 0;
  } 

  //convert benchmark coordinate to problem definition coordinate
  //benchmark description: 1: dip, 2: strike, 3: normal
  //problem setup: 1: strike, 2: normal, 3: dip
  //so just replace 1 -> 2 and 2 -> 3 and 3 -> 1 

  if ( _i == 1 && _j == 1 ){ To = -sigma22; } //we use pressure BC, so the stress is positive
  else if ( _i == 2 && _j == 2 ){ To = -sigma33; } //we use pressure BC, so the stress is positive
  else if ( _i == 3 && _j == 3 ){ To = -sigma11; } //we use pressure BC, so the stress is positive
  else if ( ( _i == 1 && _j == 2 ) || ( _i == 2 && _j == 1 ) ){ To = sigma23; }
  else{ To = 0.0; }

  return To;

}