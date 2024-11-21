#include "PoroCZMInterfaceKernelbase.h"

InputParameters
PoroCZMInterfaceKernelbase::validParams()
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
  params.addRequiredCoupledVar("displacements", "the string containing displacement variables");
  // params.addRequiredCoupledVar("fluid_vel", "the string containing fluid velocity variables");
  params.set<std::string>("traction_global_name") = "traction_global";
  params.set<std::string>("pressure_global_name") = "pressure_global";
  return params;
}

PoroCZMInterfaceKernelbase::PoroCZMInterfaceKernelbase(const InputParameters & parameters)
  : JvarMapKernelInterface<InterfaceKernel>(parameters),
    _base_name(isParamValid("base_name") && !getParam<std::string>("base_name").empty()
                   ? getParam<std::string>("base_name") + "_"
                   : ""),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _disp_neighbor_var(_ndisp),
    _vars(_ndisp),
    // _nvel(coupledComponents("fluid_vel")),
    // _fluid_vel_var(_nvel),
    // _fluid_vel_neighbor_var(_nvel),
    // _fluid_vars(_nvel),
    // _pp_var(coupled("porepressure")),        
    // _pp_neighbor_var(coupled("porepressure")), 
    // _pp_vars(nullptr), 
    _traction_global(getMaterialPropertyByName<RealVectorValue>(
        _base_name + getParam<std::string>("traction_global_name"))),
    _dtraction_djump_global(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "dtraction_djump_global")),
    // _dtraction_djump_global_vf(
    //     getMaterialPropertyByName<RankTwoTensor>(_base_name + "dtraction_djump_global_vf")),
    // _dtraction_dpressure_global(
    //     getMaterialPropertyByName<RealVectorValue>(_base_name + "dtraction_dpressure_global")),
    _pressure_global(getMaterialPropertyByName<Real>(
        _base_name + getParam<std::string>("pressure_global_name"))),
    _dpressure_djump_global(
        getMaterialPropertyByName<RealVectorValue>(_base_name + "dpressure_djump_global"))
    // _dpressure_djump_global_vf(
    //     getMaterialPropertyByName<RealVectorValue>(_base_name + "dpressure_djump_global_vf"))
{

    // Enforce consistency
  if (_ndisp != _mesh.dimension())
    paramError("displacements", "Number of displacements must match problem dimension.");

  if (_ndisp > 3 || _ndisp < 1)
    mooseError("the CZM material requires 1, 2 or 3 displacement variables");

  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp_var[i] = coupled("displacements", i);
    _disp_neighbor_var[i] = coupled("displacements", i);
    _vars[i] = getVar("displacements", i);
  }


  // // Enforce consistency
  // if (_nvel != _mesh.dimension())
  //   paramError("fluid_vel", "Number of fluid velocity must match problem dimension.");

  // if (_nvel > 3 || _nvel < 1)
  //   mooseError("the CZM material requires 1, 2 or 3 fluid velocity variables");

  // for (unsigned int i = 0; i < _nvel ; ++i)
  // {
  //   _fluid_vel_var[i] = coupled("fluid_vel", i);
  //   _fluid_vel_neighbor_var[i] = coupled("fluid_vel", i);
  //   _fluid_vars[i] = getVar("fluid_vel", i);
  // }

  //   _pp_var = coupled("porepressure");
  //   _pp_neighbor_var = coupled("porepressure");
  //   _pp_vars = getVar("porepressure",0);
}

Real
PoroCZMInterfaceKernelbase::computeQpResidual(Moose::DGResidualType type)
{
  
  Real r = 0.0;
  
  if (_component == 0)
     r = _traction_global[_qp](_component) - _pressure_global[_qp];
  else
     r = _traction_global[_qp](_component);

  switch (type)
  {
  // [test_secondary-test_primary]*(T-p) where (T-p) represents the effective traction.
    case Moose::Element:
      r *= -_test[_i][_qp];
      break;
    case Moose::Neighbor:
      r *= _test_neighbor[_i][_qp];
      break;
  }
  return r;
  }

Real
PoroCZMInterfaceKernelbase::computeQpJacobian(Moose::DGJacobianType type)
{
  // retrieve the diagonal Jacobian coefficient dependning on the displacement
  // component (_component) this kernel is working on
  return computeDResidualDDisplacement(_component, type);
}

Real
PoroCZMInterfaceKernelbase::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
  // bail out if jvar is not coupled
  if (getJvarMap()[jvar] < 0)
    return 0.0;

  // Jacobian of the residul[_component] w.r.t to the coupled displacement
  // component[off_diag_component]
  for (unsigned int off_diag_component = 0; off_diag_component < _ndisp; ++off_diag_component)
  {
    if (jvar == _disp_var[off_diag_component])
      return computeDResidualDDisplacement(off_diag_component, type);
  }
  // // Jacobian of the residul[_component] w.r.t to the coupled fluid velocity
  // // component[off_diag_component]
  // for (unsigned int off_diag_component = 0; off_diag_component < _nvel; ++off_diag_component)
  // {
  //   if (jvar == _fluid_vel_var[off_diag_component])
  //     return computeDResidualDFlux(off_diag_component, type);
  // }

  // for (unsigned int off_diag_component = 0; off_diag_component < _nvel; ++off_diag_component)
  // {
  //   if (jvar == _fluid_vel_var[off_diag_component])
  //     return computeDResidualDFlux(off_diag_component, type);
  // }

  // if (jvar == _pp_var)
  //     return computeDResidualDPressure(type);
  // // this is the place where one should implement derivatives of the residual w.r.t. other variables
  return 0.0;
}
