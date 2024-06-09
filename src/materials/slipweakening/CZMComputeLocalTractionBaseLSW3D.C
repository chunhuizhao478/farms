#include "Assembly.h"
#include "CZMComputeLocalTractionBaseLSW3D.h"

InputParameters
CZMComputeLocalTractionBaseLSW3D::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constitutive material "
                             "models that can be formulated using the total displacement jump");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

CZMComputeLocalTractionBaseLSW3D::CZMComputeLocalTractionBaseLSW3D(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _interface_traction(declarePropertyByName<RealVectorValue>(_base_name + "interface_traction")),
    _dinterface_traction_djump(
        declarePropertyByName<RankTwoTensor>(_base_name + "dinterface_traction_djump")),
    _interface_displacement_jump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _accumulated_slip_along_normal(declarePropertyByName<Real>("accumulated_slip_along_normal")),
    _accumulated_slip_along_strike(declarePropertyByName<Real>("accumulated_slip_along_strike")),
    _accumulated_slip_along_dip(declarePropertyByName<Real>("accumulated_slip_along_dip")),
    _slip_along_normal(declarePropertyByName<Real>("slip_along_normal")),
    _slip_along_strike(declarePropertyByName<Real>("slip_along_strike")),
    _slip_along_dip(declarePropertyByName<Real>("slip_along_dip"))  
{
}

void
CZMComputeLocalTractionBaseLSW3D::initQpStatefulProperties()
{
  _interface_traction[_qp] = 0;

  _accumulated_slip_along_normal[_qp] = 0.0;
  _accumulated_slip_along_strike[_qp] = 0.0;
  _accumulated_slip_along_dip[_qp] = 0.0;

  _slip_along_normal[_qp] = 0.0;
  _slip_along_strike[_qp] = 0.0;
  _slip_along_dip[_qp] = 0.0;

}

void
CZMComputeLocalTractionBaseLSW3D::computeQpProperties()
{
  computeInterfaceTractionAndDerivatives();
}