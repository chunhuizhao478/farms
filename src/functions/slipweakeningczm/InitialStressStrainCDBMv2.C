#include "InitialStressStrainCDBMv2.h"

registerMooseObject("farmsApp", InitialStressStrainCDBMv2);

InputParameters
InitialStressStrainCDBMv2::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Function to compute initial stress and strain based on the CDBM v2 model.");
  params.addRequiredParam<Real>("i", "index");
  params.addRequiredParam<Real>("j", "index");
  params.addRequiredParam<Real>("lambda_o", "initial lambda parameter for the CDBM model");
  params.addRequiredParam<Real>("shear_modulus_o", "initial shear modulus parameter for the CDBM model");
  params.addRequiredParam<Real>("fluid_density", "fluid density in kg/m^3");
  params.addRequiredParam<Real>("rock_density", "rock density in kg/m^3");
  params.addRequiredParam<Real>("gravity", "gravity in m/s^2");
  params.addRequiredParam<Real>("bxx", "coefficient for sigmaxx");
  params.addRequiredParam<Real>("byy", "coefficient for sigmayy");
  params.addRequiredParam<Real>("bxy", "coefficient for sigmaxy");
  params.addRequiredParam<Real>("cutoff_distance", "cutoff distance for the depth varying stress");
  params.addParam<Real>("peak_shear_value", -1 ,"initial shear stress perturbation peak value");
  params.addParam<Real>("nucl_center_x", -1, "nucleation center x coordinate");
  params.addParam<Real>("nucl_center_z", 1, "nucleation center z coordinate, the maximum value is 0");
  params.addParam<Real>("nucl_size", -1, "nucleation size");
  params.addParam<Real>("elem_size", -1, "element size for the simulation, used for determining the nucleation zone");
  params.addParam<bool>("get_initial_stress", false, "flag to get initial stress");
  params.addParam<bool>("get_initial_strain", false, "flag to get initial strain");
  params.addParam<bool>("get_shear_overstress", false, "flag to get initial overstress");
  params.addParam<bool>("get_fluid_pressure", false, "flag to get fluid pressure");
  return params;
}

InitialStressStrainCDBMv2::InitialStressStrainCDBMv2(const InputParameters & parameters)
  : Function(parameters),
  _i(getParam<Real>("i")),
  _j(getParam<Real>("j")),
  _lambda_o(getParam<Real>("lambda_o")),
  _shear_modulus_o(getParam<Real>("shear_modulus_o")),
  _fluid_density(getParam<Real>("fluid_density")),
  _rock_density(getParam<Real>("rock_density")),
  _gravity(getParam<Real>("gravity")),
  _bxx(getParam<Real>("bxx")),
  _byy(getParam<Real>("byy")),
  _bxy(getParam<Real>("bxy")),
  _peak_shear_value(getParam<Real>("peak_shear_value")),
  _nucl_center_x(getParam<Real>("nucl_center_x")),
  _nucl_center_z(getParam<Real>("nucl_center_z")),
  _nucl_size(getParam<Real>("nucl_size")),
  _elem_size(getParam<Real>("elem_size")),
  _cutoff_distance(getParam<Real>("cutoff_distance")),
  _get_initial_stress(getParam<bool>("get_initial_stress")),
  _get_initial_strain(getParam<bool>("get_initial_strain")),
  _get_shear_overstress(getParam<bool>("get_shear_overstress")),
  _get_fluid_pressure(getParam<bool>("get_fluid_pressure"))
{
  //some checks for parameters
  if (_get_shear_overstress && (_peak_shear_value < 0 || _nucl_center_x < 0 || _nucl_center_z > 0 || _nucl_size < 0 || _elem_size < 0)) {
    std::cout<<"peak_shear_value: "<<_peak_shear_value<<", nucl_center_x: "<<_nucl_center_x<<", nucl_center_z: "<<_nucl_center_z<<", nucl_size: "<<_nucl_size<<", elem_size: "<<_elem_size<<std::endl;
    mooseError("When get_shear_overstress is true, peak_shear_value, nucl_center_x, nucl_center_z, nucl_size, elem_size must be provided.");
  }
  if (_get_initial_stress && _get_initial_strain) {
    mooseError("Cannot get both initial stress and strain at the same time. Please choose one.");
  }
  if (_get_shear_overstress && !_get_initial_stress) {
    mooseError("When get_shear_overstress is true, get_initial_stress must also be true.");
  }
  if (_get_fluid_pressure && (_get_initial_stress || _get_initial_strain)) {
    mooseError("When get_fluid_pressure is true, get_initial_stress and get_initial_strain must be false.");
  }
}

Real
InitialStressStrainCDBMv2::value(Real /*t*/, const Point & p) const
{
  
  //Define variable takes the value
  Real var = 0.0; 

  //Compute the initial stress
  //the coordinate follows benchmark
  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction
  Real z_coord = p(2); //along the dip direction
  
  //define the parameters
  Real lambda_o = _lambda_o; //Pa
  Real shear_modulus_o = _shear_modulus_o; //Pa
  Real fluid_density = _fluid_density; //kg/m^3 fluid density
  Real rock_density = _rock_density; //kg/m^3 rock density
  Real gravity = _gravity; //m/s^2
  Real bxx = _bxx; 
  Real byy = _byy;
  Real bxy = _bxy; 
  Real peak_shear_value = _peak_shear_value; //Pa
  Real nucl_center_x = _nucl_center_x; //nucleation center x coordinate
  Real nucl_center_z = _nucl_center_z; //nucleation center z coordinate
  Real nucl_size = _nucl_size; //nucleation size
  Real elem_size = _elem_size; //element size for the simulation
  Real cutoff_distance = _cutoff_distance; //cutoff distance for the depth varying stress

  //define stress components
  Real sigmazz = 0;
  Real sigmaxx = 0;
  Real sigmayy = 0;
  Real sigmaxy = 0;
  Real sigmaxz = 0;
  Real sigmayz = 0; 

  //Pf
  Real Pf = fluid_density * gravity * abs(z_coord);

  //sigmazz
  sigmazz = -1 * rock_density * gravity * abs(z_coord);

  //sigmaxx
  if ( abs(z_coord) <= cutoff_distance ) {
    sigmaxx = bxx * ( sigmazz + Pf ) - Pf;
  }
  else{
    sigmaxx = sigmazz;
  }

  //sigmayy
  if ( abs(z_coord) <= cutoff_distance ) {
    sigmayy = byy * ( sigmazz + Pf ) - Pf;
  }
  else{
    sigmayy = sigmazz;
  } 

  //sigmaxy
  if ( abs(z_coord) <= cutoff_distance ) {
    sigmaxy = bxy * ( sigmazz + Pf );
  }
  else{
    sigmaxy = 0;
  } 

  //overstress within the region
  if (_get_shear_overstress)
  {
    //Check if the point is within the nucleation zone
    if ((x_coord <= (nucl_center_x + nucl_size * 0.5)) && (x_coord >= (nucl_center_x - nucl_size * 0.5)) &&
        (z_coord <= (nucl_center_z + nucl_size * 0.5)) && (z_coord >= (nucl_center_z - nucl_size * 0.5)) &&
        (y_coord >= -1 * elem_size) && (y_coord <= elem_size))
    {
      sigmaxy = peak_shear_value;
    }
    else
    {
      sigmaxy = sigmaxy; //use the background value if not in nucleation zone
    }
  }

  //Compute the initial strain components
  Real sigma_mean = (sigmaxx + sigmayy + sigmazz);
  Real first_hooke_law_factor = (1.0 / (2.0 * shear_modulus_o));
  Real second_hooke_law_factor = (lambda_o / (2 * shear_modulus_o * (3 * lambda_o + 2 * shear_modulus_o)));
  Real epsxx = first_hooke_law_factor * sigmaxx - second_hooke_law_factor * sigma_mean;
  Real epsyy = first_hooke_law_factor * sigmayy - second_hooke_law_factor * sigma_mean;
  Real epszz = first_hooke_law_factor * sigmazz - second_hooke_law_factor * sigma_mean;
  Real epsxy = first_hooke_law_factor * sigmaxy;
  Real epsxz = first_hooke_law_factor * sigmaxz;
  Real epsyz = first_hooke_law_factor * sigmayz;

  //output the properties
  if (_get_initial_stress)
  {
    //output the initial stress tensor
    if ( _i == 1 && _j == 1 ){ var = sigmaxx; }
    else if ( _i == 2 && _j == 2 ){ var = sigmayy; }
    else if ( _i == 3 && _j == 3 ){ var = sigmazz; }
    else if ( ( _i == 1 && _j == 2 ) || ( _i == 2 && _j == 1 ) ){ var = sigmaxy; }
    else if ( ( _i == 1 && _j == 3 ) || ( _i == 3 && _j == 1 ) ){ var = sigmaxz; }
    else if ( ( _i == 2 && _j == 3 ) || ( _i == 3 && _j == 2 ) ){ var = sigmayz; }
    else{ var = 0.0; }
  }
  else if (_get_initial_strain)
  {
    //output the initial strain tensor
    if ( _i == 1 && _j == 1 ){ var = epsxx; }
    else if ( _i == 2 && _j == 2 ){ var = epsyy; }
    else if ( _i == 3 && _j == 3 ){ var = epszz; }
    else if ( ( _i == 1 && _j == 2 ) || ( _i == 2 && _j == 1 ) ){ var = epsxy; }
    else if ( ( _i == 1 && _j == 3 ) || ( _i == 3 && _j == 1 ) ){ var = epsxz; }
    else if ( ( _i == 2 && _j == 3 ) || ( _i == 3 && _j == 2 ) ){ var = epsyz; }
    else{ var = 0.0; }
  }
  else if (_get_fluid_pressure)
  {
    var = Pf; //output the fluid pressure
  }

  return var;

}