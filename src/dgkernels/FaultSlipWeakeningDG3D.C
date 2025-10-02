#include "FaultSlipWeakeningDG3D.h"
using namespace ElasticDG3D;
registerMooseObject("farmsApp", FaultSlipWeakeningDG3D);

InputParameters FaultSlipWeakeningDG3D::validParams()
{
  return DGKernel::validParams();
}

FaultSlipWeakeningDG3D::FaultSlipWeakeningDG3D(const InputParameters & params)
  : DGKernel(params),
    _fx_sxx(getMaterialProperty<Real>("fault_flux_sxx")),
    _fx_sxy(getMaterialProperty<Real>("fault_flux_sxy")),
    _fx_sxz(getMaterialProperty<Real>("fault_flux_sxz")),
    _fx_syy(getMaterialProperty<Real>("fault_flux_syy")),
    _fx_syz(getMaterialProperty<Real>("fault_flux_syz")),
    _fx_szz(getMaterialProperty<Real>("fault_flux_szz")),
    _fx_ux(getMaterialProperty<Real>("fault_flux_ux")),
    _fx_uy(getMaterialProperty<Real>("fault_flux_uy")),
    _fx_uz(getMaterialProperty<Real>("fault_flux_uz")) {}

// Helper to fetch optional aux variable pointer
static const VariableValue * optionalCoupled(const InputParameters & p, const std::string & name)
{ return p.isParamValid(name) ? &p.get<VariableValue>(name) : nullptr; }

Real FaultSlipWeakeningDG3D::computeQpResidual(Moose::DGResidualType type)
{
  const std::string & var = _var.name();
  Real F = 0.0;
  if (var=="sxx") F=_fx_sxx[_qp]; else if (var=="sxy") F=_fx_sxy[_qp]; else if (var=="sxz") F=_fx_sxz[_qp];
  else if (var=="syy") F=_fx_syy[_qp]; else if (var=="syz") F=_fx_syz[_qp]; else if (var=="szz") F=_fx_szz[_qp];
  else if (var=="ux") F=_fx_ux[_qp]; else if (var=="uy") F=_fx_uy[_qp]; else if (var=="uz") F=_fx_uz[_qp];

  return (type==Moose::Element ? -_test[_i][_qp]*F : _test_neighbor[_i][_qp]*F);
}

Real FaultSlipWeakeningDG3D::computeQpJacobian(Moose::DGJacobianType /*type*/)
{
  // Placeholder: analytical Jacobian not yet implemented
  return 0.0;
}

Real FaultSlipWeakeningDG3D::computeQpOffDiagJacobian(Moose::DGJacobianType /*type*/, unsigned int /*jvar*/)
{
  // Placeholder: off-diagonal couplings not yet implemented
  return 0.0;
}
