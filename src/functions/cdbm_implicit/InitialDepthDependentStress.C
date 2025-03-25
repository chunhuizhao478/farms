#include "InitialDepthDependentStress.h"

registerMooseObject("farmsApp", InitialDepthDependentStress);

InputParameters
InitialDepthDependentStress::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Initial Depth Dependent Stress: xx - 11, yy - 22, zz - 33, xy - 12, pf - 00");
  params.addRequiredParam<Real>("i", "index");
  params.addRequiredParam<Real>("j", "index");
  params.addRequiredParam<bool>("pos_sign", "sign of the stress"); //0 for positive, 1 for negative
  params.addRequiredParam<Real>("fluid_density", "kg/m^3 fluid density");
  params.addRequiredParam<Real>("rock_density", "kg/m^3 rock density");
  params.addRequiredParam<Real>("gravity", "m/s^2");
  params.addRequiredParam<Real>("bxx", "coefficient for sigma_xx");
  params.addRequiredParam<Real>("byy", "coefficient for sigma_yy");
  params.addRequiredParam<Real>("bxy", "coefficient for sigma_xy");
  params.addRequiredParam<Real>("linear_variation_cutoff_distance", "linear variation cutoff distance");
  params.addParam<Real>("constant_cohesion", 0, "constant cohesion");
  params.addParam<Real>("constant_cohesion_cutoff_distance", 0, "cutoff distance for constant cohesion");
  return params;
}

InitialDepthDependentStress::InitialDepthDependentStress(const InputParameters & parameters)
  : Function(parameters),
  _i(getParam<Real>("i")),
  _j(getParam<Real>("j")),
  _sign(getParam<bool>("pos_sign")),
  _fluid_density(getParam<Real>("fluid_density")),
  _rock_density(getParam<Real>("rock_density")),
  _gravity(getParam<Real>("gravity")),
  _bxx(getParam<Real>("bxx")),
  _byy(getParam<Real>("byy")),
  _bxy(getParam<Real>("bxy")),
  _linear_variation_cutoff_distance(getParam<Real>("linear_variation_cutoff_distance")),
  _constant_cohesion(getParam<Real>("constant_cohesion")),
  _constant_cohesion_cutoff_distance(getParam<Real>("constant_cohesion_cutoff_distance"))
{
}

Real
InitialDepthDependentStress::value(Real /*t*/, const Point & p) const
{
  
  //the coordinate follows benchmark
  Real z_coord = p(2); //along the dip direction
  Real To = 0;
  Real fluid_density = _fluid_density; //kg/m^3 fluid density
  Real rock_density = _rock_density; //kg/m^3 rock density
  Real gravity = _gravity; //m/s^2
  Real bxx = _bxx; //coefficient for sigma_xx
  Real byy = _byy; //coefficient for sigma_yy
  Real bxy = _bxy; //coefficient for sigma_xy
  Real sigma_zz = 0;
  Real sigma_xx = 0;
  Real sigma_yy = 0;
  Real sigma_xy = 0;

  //Pf
  Real Pf = fluid_density * gravity * abs(z_coord);

  //sigma_zz
  sigma_zz = -1 * rock_density * gravity * abs(z_coord);

  //sigma_xx
  if ( abs(z_coord) <= _linear_variation_cutoff_distance ) {
    sigma_xx = bxx * ( sigma_zz + Pf ) - Pf;
  }
  else{
    sigma_xx = sigma_zz;
  }

  //sigma_yy
  if ( abs(z_coord) <= _linear_variation_cutoff_distance ) {
    sigma_yy = byy * ( sigma_zz + Pf ) - Pf;
  }
  else{
    sigma_yy = sigma_zz;
  } 

  //sigma_xy
  if ( abs(z_coord) <= _linear_variation_cutoff_distance ) {
    sigma_xy = bxy * ( sigma_zz + Pf );
  }
  else{
    sigma_xy = 0;
  } 

  if (_sign){
  
    if ( _i == 1 && _j == 1 ){ To = -sigma_xx; } 
    else if ( _i == 2 && _j == 2 ){ To = -sigma_yy; } 
    else if ( _i == 3 && _j == 3 ){ To = -sigma_zz; } 
    else if ( ( _i == 1 && _j == 2 ) || ( _i == 2 && _j == 1 ) ){ To = -sigma_xy; }
    else{ To = 0.0; }
  
  }
  else{

    if ( _i == 1 && _j == 1 ){ To = sigma_xx; } 
    else if ( _i == 2 && _j == 2 ){ To = sigma_yy; } 
    else if ( _i == 3 && _j == 3 ){ To = sigma_zz; } 
    else if ( ( _i == 1 && _j == 2 ) || ( _i == 2 && _j == 1 ) ){ To = sigma_xy; }
    else{ To = 0.0; }

  }

  //return fluid_density * gravity * abs(z_coord) if _i = 0, _j = 0
  if ( _i == 0 && _j == 0 ){ To = Pf; }

  return To;

}