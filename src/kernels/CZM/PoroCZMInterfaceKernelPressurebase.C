#include "PoroCZMInterfaceKernelPressurebase.h"

InputParameters
PoroCZMInterfaceKernelPressurebase::validParams()
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
  params.addRequiredCoupledVar("p", "the string containing pressure");
  return params;
}

PoroCZMInterfaceKernelPressurebase::PoroCZMInterfaceKernelPressurebase(const InputParameters & parameters)
  : JvarMapKernelInterface<InterfaceKernel>(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _component(getParam<unsigned int>("component")),
    _p_var_num(coupled("p")),
    _p_var_num_porepressure_neighbor_var(coupled("p")),
    _p_var(*getVar("p", 0)),
    _phi(_assembly.phiFace(_p_var)),
    _phi_neighbor(_assembly.phiFaceNeighbor(_p_var)),
    _dtraction_dpressure_global(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "dtraction_dpressure_global"))
{

}

Real
PoroCZMInterfaceKernelPressurebase::computeQpResidual(Moose::DGResidualType type)
{
  return 0.0;
}

Real
PoroCZMInterfaceKernelPressurebase::computeQpJacobian(Moose::DGJacobianType type)
{
  // retrieve the diagonal Jacobian coefficient dependning on the displacement
  // component (_component) this kernel is working on
  return 0.0;
  
}

Real
PoroCZMInterfaceKernelPressurebase::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  if (_p_var_num == jvar)
  {
     return computeDResidualDPressure(type);
  }
   
  return 0.0;
}
