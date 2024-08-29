#include "InitialStressTPV243Dcdbm.h"

registerMooseObject("farmsApp", InitialStressTPV243Dcdbm);

InputParameters
InitialStressTPV243Dcdbm::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<Real>("i", "index");
  params.addRequiredParam<Real>("j", "index");
  return params;
}

InitialStressTPV243Dcdbm::InitialStressTPV243Dcdbm(const InputParameters & parameters)
  : Function(parameters),
  _i(getParam<Real>("i")),
  _j(getParam<Real>("j"))
{
}

Real
InitialStressTPV243Dcdbm::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real y_coord = p(1); //along the dip direction
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
  Real Pf = fluid_density * gravity * abs(y_coord);

  //sigma11
  sigma11 = -1 * rock_density * gravity * abs(y_coord);

  //sigma22
  if ( abs(y_coord) <= 15600 ) {
    sigma22 = b22 * ( sigma11 + Pf ) - Pf;
  }
  else{
    sigma22 = sigma11;
  }

  //sigma33
  if ( abs(y_coord) <= 15600 ) {
    sigma33 = b33 * ( sigma11 + Pf ) - Pf;
  }
  else{
    sigma33 = sigma11;
  } 

  //sigma23
  if ( abs(y_coord) <= 15600 ) {
    sigma23 = b23 * ( sigma11 + Pf );
  }
  else{
    sigma23 = 0;
  } 

  //overstress within the region
  Real x_coord = p(0);
  if (x_coord >= -9500 and x_coord <= -6500 and y_coord >= -11500 and y_coord <= -8500){
    sigma23 = 5.7210e7;
  }

  //convert benchmark coordinate to problem definition coordinate
  //benchmark description: 1: dip, 2: strike, 3: normal
  //problem setup: 1: strike, 2: dip, 3: normal
  //so just replace 1 -> 2 and 2 -> 1 

  if ( _i == 1 && _j == 1 ){ To = sigma22; }
  else if ( _i == 2 && _j == 2 ){ To = sigma11; }
  else if ( _i == 3 && _j == 3 ){ To = sigma33; }
  else if ( ( _i == 1 && _j == 3 ) || ( _i == 3 && _j == 1 ) ){ To = sigma23; }
  else{ To = 0.0; }

  return To;

}