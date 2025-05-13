//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NDPhaseFieldPsiDerivatives.h"

/**
 * Class to evaluate the phase field psi function and other related properties.
 * This class is used for non-AD implementation of the phase field model.
 * Created by Chunhui Zhao on 2025-5-11.
 */
registerMooseObject("farmsApp", NDPhaseFieldPsiDerivatives);

InputParameters
NDPhaseFieldPsiDerivatives::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "This is a non-AD implementation of the phase field psi function and other related properties.");

  params.addRequiredCoupledVar("phase_field", "The phase-field variable");
  params.addRequiredCoupledVar("psie_active", "The active strain energy density");

  params.addRequiredParam<std::string>(
      "model_type", "The type of the model: AT1, AT2, PF_CZM");

  params.addParam<MaterialPropertyName>("normalization_constant",
                                        "c0",
                                        "Name of the material to store the normalization constant, "
                                        "$4\\int_0^1 \\sqrt{\\alpha(s)} \\diff{s}$");
  params.addParam<MaterialPropertyName>("alpha", "alpha", "Name of the material to store alpha(d)");
  params.addParam<MaterialPropertyName>("dalpha_dd", "dalpha_dd", "Name of the material to store "
                                        "dalpha/dd");
  params.addParam<MaterialPropertyName>("d2alpha_dd2", "d2alpha_dd2", "Name of the material to store "
                                        "d2alpha/dd2");
  params.addParam<MaterialPropertyName>("g", "g", "Name of the material to store g(d)");
  params.addParam<MaterialPropertyName>("dg_dd", "dg_dd", "Name of the material to store dg/dd");
  params.addParam<MaterialPropertyName>(
      "d2g_dd2", "d2g_dd2", "Name of the material to store d2g/dd2");
  params.addParam<MaterialPropertyName>("psi", "psi", "Name of the material to store psi(d)");
  params.addParam<MaterialPropertyName>(
      "dpsi_dd", "dpsi_dd", "Name of the material to store dpsi/dd");
  params.addParam<MaterialPropertyName>(
      "d2psi_dd2", "d2psi_dd2", "Name of the material to store d2psi/dd2");
  params.addParam<MaterialPropertyName>(
      "Gc", "Gc", "Name of the material to store the energy release rate");
  params.addParam<MaterialPropertyName>(
      "l", "l", "Name of the material to store the regularization length");

  params.addRequiredParam<Real>("eta", "Parameter in the degradation function");
  params.addParam<Real>(
      "tolerance", 1e-8, "Tolerance of the numerically computed normalization constant");
  params.addParam<unsigned int>("maximum_iterations",
                                1e9,
                                "Maximum number of iterations allowed for the numerical "
                                "computation of the normalization constant");
  return params;
}

NDPhaseFieldPsiDerivatives::NDPhaseFieldPsiDerivatives(const InputParameters & parameters)
  : Material(parameters),
    _d(coupledValue("phase_field")),
    _psie_active(coupledValue("psie_active")),
    _model_type(getParam<std::string>("model_type")),
    _c0(declareProperty<Real>(getParam<MaterialPropertyName>("normalization_constant"))),
    _alpha(declareProperty<Real>(getParam<MaterialPropertyName>("alpha"))),
    _dalpha_dd(declareProperty<Real>(getParam<MaterialPropertyName>("dalpha_dd"))),
    _d2alpha_dd2(declareProperty<Real>(getParam<MaterialPropertyName>("d2alpha_dd2"))),
    _g(declareProperty<Real>(getParam<MaterialPropertyName>("g"))),
    _dg_dd(declareProperty<Real>(getParam<MaterialPropertyName>("dg_dd"))),
    _d2g_dd2(declareProperty<Real>(getParam<MaterialPropertyName>("d2g_dd2"))),
    _psi(declareProperty<Real>(getParam<MaterialPropertyName>("psi"))),
    _dpsi_dd(declareProperty<Real>(getParam<MaterialPropertyName>("dpsi_dd"))),
    _d2psi_dd2(declareProperty<Real>(getParam<MaterialPropertyName>("d2psi_dd2"))),
    _eta(getParam<Real>("eta")),
    _Gc(getMaterialProperty<Real>(getParam<MaterialPropertyName>("Gc"))),
    _l(getMaterialProperty<Real>(getParam<MaterialPropertyName>("l"))),
    _tolerance(getParam<Real>("tolerance")),
    _max_its(getParam<unsigned int>("maximum_iterations"))
{
  _c0_0 = computeNormalizationConstant();
}

void
NDPhaseFieldPsiDerivatives::computeQpProperties()
{
  // Compute the normalization constant
  _c0[_qp] = _c0_0;

  // Compute the derivatives of alpha, g, and psi
  computeAlphaDerivatives();
  computeGDerivatives();
  computePsiDerivatives();
}

void
NDPhaseFieldPsiDerivatives::computeAlphaDerivatives()
{ 
  if (_model_type == "AT2"){
    _alpha[_qp] = _d[_qp] * _d[_qp];
    _dalpha_dd[_qp] = 2 * _d[_qp];
    _d2alpha_dd2[_qp] = 2;
  }
  else if (_model_type == "AT1"){
    _alpha[_qp] = _d[_qp];
    _dalpha_dd[_qp] = 1;
    _d2alpha_dd2[_qp] = 0;
  }
  else if (_model_type == "PF_CZM"){
    _alpha[_qp] = 2 * _d[_qp] - _d[_qp] * _d[_qp];
    _dalpha_dd[_qp] = 2 - 2 * _d[_qp];
    _d2alpha_dd2[_qp] = -2;
  }
  else
    mooseError("Unknown model type: " + _model_type);
}

void
NDPhaseFieldPsiDerivatives::computeGDerivatives()
{
  if (_model_type == "AT2"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "AT1"){
    _g[_qp] = std::pow( (1 - _d[_qp]) , 2 ) * (1 - _eta) + _eta;
    _dg_dd[_qp] = -2 * (1 - _eta) * (1 - _d[_qp]);
    _d2g_dd2[_qp] = 2 * (1 - _eta);
  }
  else if (_model_type == "PF_CZM"){
    mooseError("PF_CZM model is not implemented yet.");
  }
  else
    mooseError("Unknown model type: " + _model_type);
}

void
NDPhaseFieldPsiDerivatives::computePsiDerivatives()
{
  /*
    * The free energy density psi(d) is computed as:
    * psi(d) = alpha(d) * Gc / ( c0 * l ) + g(d) * psie_active
    * We evaluate the first derivative of psi(d) with respect to d:
    * dpsi/dd = dalpha/dd * Gc / ( c0 * l ) + dg/dd * psie_active
    * This term is used in the PFFSource kernel to compute the source term
  */
  // Compute the psi and its derivatives
  _psi[_qp] = _alpha[_qp] * _Gc[_qp] / (_c0[_qp] * _l[_qp]) + _g[_qp] * _psie_active[_qp];
  // dpsi/dd
  _dpsi_dd[_qp] = _dalpha_dd[_qp] * _Gc[_qp] / (_c0[_qp] * _l[_qp]) + _dg_dd[_qp] * _psie_active[_qp];
  // d2psi/dd2
  _d2psi_dd2[_qp] = _d2alpha_dd2[_qp] * _Gc[_qp] / (_c0[_qp] * _l[_qp]) + _d2g_dd2[_qp] * _psie_active[_qp];
}

Real
NDPhaseFieldPsiDerivatives::computeNormalizationConstant()
{
  // We use an adaptive trapezoidal rule to integrate the normalization constant to a target
  // precision

  /*
    * It basically computes c0 = 4 * \int_0^1 \sqrt{alpha(d)} dd
  */
  std::stack<std::pair<Real, Real>> S;
  S.emplace(0, 1);
  Real I = 0;
  unsigned int its = 0;
  while (!S.empty())
  {
    auto interval = S.top();
    S.pop();

    Real a = interval.first;
    Real b = interval.second;
    Real m = (a + b) / 2;

    Real I1 = (normalizationIntegrand(a) + normalizationIntegrand(b)) * (b - a) / 2;
    Real I2 =
        (normalizationIntegrand(a) + 2 * normalizationIntegrand(m) + normalizationIntegrand(b)) *
        (b - a) / 4;

    if (std::abs(I1 - I2) < 3 * (b - a) * _tolerance)
      I += I2;
    else
    {
      S.emplace(a, m);
      S.emplace(m, b);
      its++;
    }

    if (its >= _max_its)
      mooseError("Maximum number of iterations reached, but the crack geometric function still "
                 "hasn't converge.");
  }
  return I;
}

Real
NDPhaseFieldPsiDerivatives::normalizationIntegrand(const Real & d)
{
  Real alpha_d;
  if (_model_type == "AT2")
    alpha_d = d * d;
  else if (_model_type == "AT1") 
    alpha_d = d;
  else if (_model_type == "PF_CZM")
    alpha_d = 2 * d - d * d;
  else
    mooseError("Unknown model type: " + _model_type);
    
  return 4.0 * std::sqrt(alpha_d);
}