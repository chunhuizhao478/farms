#include "Assembly.h"
#include "CZMComputeLocalTractionBaseRSF2D.h"

InputParameters
CZMComputeLocalTractionBaseRSF2D::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constitutive material "
                             "models that can be formulated using the total displacement jump");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addRequiredParam<Real>("Tn_o","initial normal traction");
  params.addRequiredParam<Real>("Ts_o","initial shear traction");
  params.addRequiredParam<Real>("Vini","initial velocity");
  params.addRequiredParam<Real>("statevarini","initial state variable");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

CZMComputeLocalTractionBaseRSF2D::CZMComputeLocalTractionBaseRSF2D(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _interface_traction(declarePropertyByName<RealVectorValue>(_base_name + "interface_traction")),
    _dinterface_traction_djump(
        declarePropertyByName<RankTwoTensor>(_base_name + "dinterface_traction_djump")),
    _interface_displacement_jump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _alongfaultvel_strike_plus(declarePropertyByName<Real>("alongfaultvel_strike_plus")),
    _alongfaultvel_strike_minus(declarePropertyByName<Real>("alongfaultvel_strike_minus")),
    _alongfaultvel_normal_plus(declarePropertyByName<Real>("alongfaultvel_normal_plus")),
    _alongfaultvel_normal_minus(declarePropertyByName<Real>("alongfaultvel_normal_minus")),
    _alongfaultdisp_strike_plus(declarePropertyByName<Real>("alongfaultdisp_strike_plus")),
    _alongfaultdisp_strike_minus(declarePropertyByName<Real>("alongfaultdisp_strike_minus")),
    _alongfaultdisp_normal_plus(declarePropertyByName<Real>("alongfaultdisp_normal_plus")),
    _alongfaultdisp_normal_minus(declarePropertyByName<Real>("alongfaultdisp_normal_minus")),
    _sliprate_strike(declarePropertyByName<Real>("sliprate_strike")),
    _slip_strike(declarePropertyByName<Real>("slip_strike")),
    _sliprate_normal(declarePropertyByName<Real>("sliprate_normal")),
    _slip_normal(declarePropertyByName<Real>("slip_normal")),
    _sliprate_mag(declarePropertyByName<Real>("sliprate_mag")),
    _slip_mag(declarePropertyByName<Real>("slip_mag")),
    _statevar(declarePropertyByName<Real>("statevar")),
    _traction_strike(declarePropertyByName<Real>("traction_strike")),
    _traction_normal(declarePropertyByName<Real>("traction_normal")),
    _sliprate_predict(declarePropertyByName<Real>("sliprate_predict")),
    _slip_predict(declarePropertyByName<Real>("slip_predict")),
    _Ts_o(getParam<Real>("Ts_o")),
    _Tn_o(getParam<Real>("Tn_o")),
    _Vini(getParam<Real>("Vini")),
    _statevarini(getParam<Real>("statevarini")),
    _alongfaultvel_x_plus(declarePropertyByName<Real>("alongfaultvel_x_plus")),
    _alongfaultvel_x_minus(declarePropertyByName<Real>("alongfaultvel_x_minus")),
    _alongfaultvel_y_plus(declarePropertyByName<Real>("alongfaultvel_y_plus")),
    _alongfaultvel_y_minus(declarePropertyByName<Real>("alongfaultvel_y_minus")),
    _alongfaultdisp_x_plus(declarePropertyByName<Real>("alongfaultdisp_x_plus")),
    _alongfaultdisp_x_minus(declarePropertyByName<Real>("alongfaultdisp_x_minus")),
    _alongfaultdisp_y_plus(declarePropertyByName<Real>("alongfaultdisp_y_plus")),
    _alongfaultdisp_y_minus(declarePropertyByName<Real>("alongfaultdisp_y_minus"))
{
}

void
CZMComputeLocalTractionBaseRSF2D::initQpStatefulProperties()
{
  _interface_traction[_qp] = 0;

  _alongfaultvel_strike_plus[_qp] = 0;
  _alongfaultvel_strike_minus[_qp] = 0;
  _alongfaultvel_normal_plus[_qp] = 0;
  _alongfaultvel_normal_minus[_qp] = 0;
  _alongfaultdisp_strike_plus[_qp] = 0;
  _alongfaultdisp_strike_minus[_qp] = 0;
  _alongfaultdisp_normal_plus[_qp] = 0;
  _alongfaultdisp_normal_minus[_qp] = 0;
  
  _sliprate_strike[_qp] = _Vini;
  _slip_strike[_qp] = 0;
  
  _sliprate_normal[_qp] = 0;
  _slip_normal[_qp] = 0;
  
  _sliprate_mag[_qp] = _Vini;
  _slip_mag[_qp] = 0;

  _statevar[_qp] = _statevarini;

  _traction_strike[_qp] = _Ts_o;
  _traction_normal[_qp] = _Tn_o; //positive number

  _sliprate_predict[_qp] = _Vini;

  _alongfaultvel_x_plus[_qp] = _Vini*0.5;
  _alongfaultvel_x_minus[_qp] = -_Vini*0.5;
  _alongfaultvel_y_plus[_qp] = 0;
  _alongfaultvel_y_minus[_qp] = 0;
  _alongfaultdisp_x_plus[_qp] = 0;
  _alongfaultdisp_x_minus[_qp] = 0;
  _alongfaultdisp_y_plus[_qp] = 0;
  _alongfaultdisp_y_minus[_qp] = 0;
}

void
CZMComputeLocalTractionBaseRSF2D::computeQpProperties()
{
  computeInterfaceTractionAndDerivatives();
}