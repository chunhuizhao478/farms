#include "Assembly.h"
#include "CZMComputeLocalTractionBaseSWF2D.h"

InputParameters
CZMComputeLocalTractionBaseSWF2D::validParams()
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

CZMComputeLocalTractionBaseSWF2D::CZMComputeLocalTractionBaseSWF2D(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _interface_traction(declarePropertyByName<RealVectorValue>(_base_name + "interface_traction")),
    _dinterface_traction_djump(
        declarePropertyByName<RankTwoTensor>(_base_name + "dinterface_traction_djump")),
    _interface_displacement_jump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _flag_track_opening(declarePropertyByName<Real>("flag_track_opening")),
    _flag_track_activecase(declarePropertyByName<Real>("flag_track_activecase")),
    _jump_track_opening(declarePropertyByName<Real>("jump_track_opening")),
    _jump_track_reversal(declarePropertyByName<Real>("jump_track_reversal")),
    _T1(declarePropertyByName<Real>("T1")),
    _T2(declarePropertyByName<Real>("T2")),
    _jump_effective(declarePropertyByName<Real>("jump_effective"))
{
}

void
CZMComputeLocalTractionBaseSWF2D::initQpStatefulProperties()
{
  _interface_traction[_qp] = 0;

  _flag_track_opening[_qp] = 0;
  _flag_track_activecase[_qp] = 0;
  _jump_track_opening[_qp] = 0;
  _jump_track_reversal[_qp] = 0;
  _T1[_qp] = 0;
  _T2[_qp] = 0;
  _jump_effective[_qp] = 0;
}

void
CZMComputeLocalTractionBaseSWF2D::computeQpProperties()
{
  computeInterfaceTractionAndDerivatives();
}