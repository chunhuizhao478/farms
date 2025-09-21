#include "PerturbationRadialAdvanced.h"
#include <cmath>

registerMooseObject("farmsApp", PerturbationRadialAdvanced);

InputParameters PerturbationRadialAdvanced::validParams(){
  InputParameters params = Material::validParams();
  params.addClassDescription("Advanced radial perturbation with smooth ramps and optional mean stress");
  params.addRequiredParam<std::vector<Real>>("nucl_center","nucleation center (x,y,z)");
  params.addRequiredParam<Real>("peak_value","Peak amplitude");
  params.addRequiredParam<Real>("thickness","Thickness in y");
  params.addRequiredParam<Real>("length","Characteristic length in x,z");
  params.addRequiredParam<Real>("duration","Ramp duration");
  params.addRequiredParam<std::string>("perturbation_type","damage | shear_stress | mean_stress");
  params.addParam<Real>("sigma_divisor",2.0,"Sigma = length / sigma_divisor");
  params.addParam<std::string>("temporal_ramp","cosine","linear | cosine | smoothstep | quintic");
  params.addParam<bool>("hold_after_duration",true,"Hold peak after duration");
  params.addParam<std::string>("thickness_taper","cosine","none | cosine");
  params.addParam<bool>("zero_outside_band",true,"Zero outside thickness each step instead of carry old value");
  return params;
}

PerturbationRadialAdvanced::PerturbationRadialAdvanced(const InputParameters & p)
 : Material(p),
   _damage_perturbation(declareProperty<Real>("damage_perturbation")),
   _shear_stress_perturbation(declareProperty<Real>("shear_stress_perturbation")),
   _mean_stress_perturbation(declareProperty<Real>("mean_stress_perturbation")),
   _damage_perturbation_old(getMaterialPropertyOldByName<Real>("damage_perturbation")),
   _shear_stress_perturbation_old(getMaterialPropertyOldByName<Real>("shear_stress_perturbation")),
   _mean_stress_perturbation_old(getMaterialPropertyOldByName<Real>("mean_stress_perturbation")),
   _nucl_center_mat(declareProperty<std::vector<Real>>("nucl_center_mat")),
   _thickness_mat(declareProperty<Real>("thickness_mat")),
   _length_mat(declareProperty<Real>("length_mat")),
   _nucl_center(getParam<std::vector<Real>>("nucl_center")),
   _peak_value(getParam<Real>("peak_value")),
   _thickness(getParam<Real>("thickness")),
   _length(getParam<Real>("length")),
   _duration(getParam<Real>("duration")),
   _perturbation_type(getParam<std::string>("perturbation_type")),
   _sigma_divisor(getParam<Real>("sigma_divisor")),
   _temporal_ramp(getParam<std::string>("temporal_ramp")),
   _hold_after_duration(getParam<bool>("hold_after_duration")),
   _thickness_taper(getParam<std::string>("thickness_taper")),
   _zero_outside_band(getParam<bool>("zero_outside_band"))
{
  if (_nucl_center.size() != 3)
    mooseError("nucl_center must have size 3");
}

void PerturbationRadialAdvanced::initQpStatefulProperties(){
  _damage_perturbation[_qp] = 0.0;
  _shear_stress_perturbation[_qp] = 0.0;
  _mean_stress_perturbation[_qp] = 0.0;
  _nucl_center_mat[_qp] = _nucl_center;
  _thickness_mat[_qp] = _thickness;
  _length_mat[_qp] = _length;
}

Real PerturbationRadialAdvanced::rampFactor(Real t) const {
  if (t <= 0.0) return 0.0;
  if (t >= _duration) return 1.0;
  Real x = t / _duration;
  if (_temporal_ramp == "linear")
    return x;
  else if (_temporal_ramp == "cosine")
    return 0.5 * (1.0 - std::cos(M_PI * x));
  else if (_temporal_ramp == "smoothstep")
    return x * x * (3.0 - 2.0 * x);
  else if (_temporal_ramp == "quintic")
    return x*x*x*(10 + x*(-15 + 6*x));
  else
    return 0.5 * (1.0 - std::cos(M_PI * x));
}

Real PerturbationRadialAdvanced::thicknessFactor(Real y) const {
  Real ymin = _nucl_center[1] - 0.5 * _thickness;
  Real ymax = _nucl_center[1] + 0.5 * _thickness;
  if (y < ymin || y > ymax)
    return 0.0;
  if (_thickness_taper == "cosine") {
    Real s = (y - ymin) / _thickness; // 0..1
    return 0.5 * (1.0 - std::cos(M_PI * s));
  }
  return 1.0;
}

void PerturbationRadialAdvanced::computeQpProperties(){
  // reset each step
  _damage_perturbation[_qp] = 0.0;
  _shear_stress_perturbation[_qp] = 0.0;
  _mean_stress_perturbation[_qp] = 0.0;

  const Real x = _q_point[_qp](0);
  const Real y = _q_point[_qp](1);
  const Real z = _q_point[_qp](2);

  const Real sigma_x = _length / _sigma_divisor;
  const Real sigma_z = _length / _sigma_divisor;

  const Real dx = x - _nucl_center[0];
  const Real dz = z - _nucl_center[2];

  const Real gaussian = _peak_value * std::exp(-(dx*dx)/(2.0*sigma_x*sigma_x) - (dz*dz)/(2.0*sigma_z*sigma_z));

  Real amp = 0.0;
  if (_t <= _duration)
    amp = gaussian * rampFactor(_t);
  else if (_hold_after_duration)
    amp = gaussian;

  // spatial thickness taper
  amp *= thicknessFactor(y);

  if (_perturbation_type == "damage")
    _damage_perturbation[_qp] = amp;
  else if (_perturbation_type == "shear_stress")
    _shear_stress_perturbation[_qp] = amp;
  else if (_perturbation_type == "mean_stress")
    _mean_stress_perturbation[_qp] = amp;
  else
    mooseError("Invalid perturbation_type: " + _perturbation_type);

  _nucl_center_mat[_qp] = _nucl_center;
  _thickness_mat[_qp] = _thickness;
  _length_mat[_qp] = _length;
}
