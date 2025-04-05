#include "Assembly.h"
#include "CZMComputeLocalTractionBaseQDRSF2D.h"

InputParameters
CZMComputeLocalTractionBaseQDRSF2D::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constitutive material "
                             "models that can be formulated using the total displacement jump");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addRequiredParam<Real>("Vini","initial velocity");
  params.addRequiredParam<Real>("statevarini","initial state variable");
  params.addRequiredParam<Real>("V_o", "Reference slip rate");
  params.addRequiredParam<Real>("f_o", "Initial friction coefficient");
  params.addRequiredParam<Real>("a", "Direct effect parameter");
  params.addRequiredParam<Real>("b", "State variable evolution parameter");
  params.addRequiredParam<Real>("L", "State variable characteristic distance");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

CZMComputeLocalTractionBaseQDRSF2D::CZMComputeLocalTractionBaseQDRSF2D(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _V_o(getParam<Real>("V_o")),
    _f_o(getParam<Real>("f_o")),
    _a(getParam<Real>("a")),
    _b(getParam<Real>("b")),
    _L(getParam<Real>("L")),
    _statevarini(getParam<Real>("statevarini")),
    _Vini(getParam<Real>("Vini")),
    _sliprate(declarePropertyByName<Real>("sliprate")),
    _statevar(declarePropertyByName<Real>("statevar")),
    _statevar_dot(declarePropertyByName<Real>("statevar_dot")),
    _slip(declarePropertyByName<Real>("slip"))

{
}

void
CZMComputeLocalTractionBaseQDRSF2D::initQpStatefulProperties()
{
  
  _sliprate[_qp] = _Vini;
  _slip[_qp] = 0;

  _statevar[_qp] = _statevarini;
  // _statevar_dot[_qp] = 1 - (_statevarini * _Vini / _L);
   Real f_new = _a * std::asinh(_Vini/ (2.0 * _V_o) * std::exp(_statevarini/ _a));
  _statevar_dot[_qp] = - _Vini/ _L * (f_new - _f_o + (_b - _a) * std::log(_Vini / _V_o));;

  

}

void
CZMComputeLocalTractionBaseQDRSF2D::computeQpProperties()
{
  computeInterfaceTractionAndDerivatives();
}