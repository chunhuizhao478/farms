#include "PoroCZMInterfaceKernelFluxbase.h"

InputParameters
PoroCZMInterfaceKernelFluxbase::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<unsigned int>("component",
                                        "the component of the "
                                        "displacement or fluid velocity vector this kernel is working on:"
                                        " component == 0, ==> X"
                                        " component == 1, ==> Y"
                                        " component == 2, ==> Z");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("fluid_vel", "the string containing fluid velocity variables");
  return params;
}

PoroCZMInterfaceKernelFluxbase::PoroCZMInterfaceKernelFluxbase(const InputParameters & parameters)
  : JvarMapKernelInterface<InterfaceKernel>(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _component(getParam<unsigned int>("component")),
    _nvel(coupledComponents("fluid_vel")),
    _fluid_vel_var(_nvel),
    _fluid_vel_neighbor_var(_nvel),
    _fluid_vars(_nvel),
    _dtraction_djump_global_vf(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "dtraction_djump_global_vf")),
    _dpressure_djump_global_vf(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "dpressure_djump_global_vf"))
{


  // Enforce consistency
  if (_nvel != _mesh.dimension())
    paramError("fluid_vel", "Number of fluid velocity must match problem dimension.");

  if (_nvel > 3 || _nvel < 1)
    mooseError("the CZM material requires 1, 2 or 3 fluid velocity variables");

  for (unsigned int i = 0; i < _nvel ; ++i)
  {
    _fluid_vel_var[i] = coupled("fluid_vel", i);
    _fluid_vel_neighbor_var[i] = coupled("fluid_vel", i);
    _fluid_vars[i] = getVar("fluid_vel", i);
  }
}

Real
PoroCZMInterfaceKernelFluxbase::computeQpResidual(Moose::DGResidualType type)
{
  return 0.0;
}

Real
PoroCZMInterfaceKernelFluxbase::computeQpJacobian(Moose::DGJacobianType type)
{
  // retrieve the diagonal Jacobian coefficient dependning on the displacement
  // component (_component) this kernel is working on
  // return computeDResidualDFlux(_component, type);
  return 0.0;
}

Real
PoroCZMInterfaceKernelFluxbase::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  // bail out if jvar is not coupled
  if (getJvarMap()[jvar] < 0)
    return 0.0;
  // Jacobian of the residul[_component] w.r.t to the coupled fluid velocity
  // component[off_diag_component]
  for (unsigned int off_diag_component = 0; off_diag_component < _nvel; ++off_diag_component)
  {
    if (jvar == _fluid_vel_var[off_diag_component])
      return computeDResidualDFlux(off_diag_component, type);
  }
  return 0.0;
}
