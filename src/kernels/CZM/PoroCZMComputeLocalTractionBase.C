#include "Assembly.h"
#include "PoroCZMComputeLocalTractionBase.h"

InputParameters
PoroCZMComputeLocalTractionBase::validParams()
{
  InputParameters params = InterfaceMaterial::validParams();

  params.addClassDescription("Base class for implementing cohesive zone constitutive material "
                             "models that can be formulated using the total fluid velocity jump");
  params.addRequiredCoupledVar("fluid_vel",
                               "The string of fluid velocity suitable for the problem statement");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addRequiredCoupledVar("porepressure",
                               "The string of porepressure suitable for the problem statement");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  return params;
}

PoroCZMComputeLocalTractionBase::PoroCZMComputeLocalTractionBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _interface_traction(declarePropertyByName<RealVectorValue>(_base_name + "interface_traction")),
    _dinterface_traction_djump(
        declarePropertyByName<RankTwoTensor>(_base_name + "dinterface_traction_djump")),
    _dinterface_traction_djump_vf(
        declarePropertyByName<RankTwoTensor>(_base_name + "dinterface_traction_djump_vf")),
    _dinterface_traction_dpressure(
        declarePropertyByName<RealVectorValue>(_base_name + "dinterface_traction_dpressure")),
    _interface_displacement_jump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_displacement_jump")),
    _interface_pressure(declarePropertyByName<Real>(_base_name + "interface_pressure")),
    _dinterface_pressure_djump(declarePropertyByName<RealVectorValue>(_base_name + "dinterface_pressure_djump")),
    _dinterface_pressure_djump_vf(
        declarePropertyByName<RealVectorValue>(_base_name + "dinterface_pressure_djump_vf")),
    // _interface_type(getParam<MooseEnum>("interfacetype").getEnum<InterfaceType>()),
    _interface_fluid_vel_jump(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "interface_fluid_vel_jump"))
{
}

void
PoroCZMComputeLocalTractionBase::initQpStatefulProperties()
{
    _interface_traction[_qp] = 0;
    _interface_pressure[_qp] = 0;
}

void
PoroCZMComputeLocalTractionBase::computeQpProperties()
{
     computeInterfaceTractionAndDerivatives(); 
}