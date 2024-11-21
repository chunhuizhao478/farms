
#include "ComputeFluidVelocityDivergence.h"
#include "MooseMesh.h"

#include "libmesh/quadrature.h"

registerMooseObject("farmsApp", ComputeFluidVelocityDivergence);

InputParameters
ComputeFluidVelocityDivergence::validParams()
{
  InputParameters params = PorousFlowMaterialVectorBase::validParams();
  params.addParam<std::string>("base_name",
                               "This should be the same base_name as given to the TensorMechanics "
                               "object that computes strain");
  params.addRequiredCoupledVar(
      "darcy_velocities",
      "The darcy_velocities appropriate for the simulation geometry and coordinate system");
  params.addClassDescription(
      "Compute divergence of darcy velocities, for use in PorousFlow.");
  params.set<std::string>("pf_material_type") = "divergence_darcy_vel";
  params.set<bool>("stateful_darcy_velocities") = true;
  params.set<bool>("at_nodes") = false;
  return params;
}

ComputeFluidVelocityDivergence::ComputeFluidVelocityDivergence(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _darcy_vel_grad(getMaterialProperty<RankTwoTensor>(_base_name + "darcy_vel_grad")),
    _nvel(coupledComponents("darcy_velocities")),
    _vel_var_num(coupledIndices("darcy_velocities")),
    _darcy_vel_div_qp(declareProperty<Real>("PorousFlow_divergence_darcy_velocity_qp")),
    _ddarcy_vel_div_qp_dvar(
        declareProperty<std::vector<RealGradient>>("dPorousFlow_divergence_darcy_velocity_qp"))
{
  if (_nvel != _mesh.dimension())
    paramError("darcy_velocities", "The number of variables supplied must match the mesh dimension.");

  if (_nodal_material)
    mooseError("DarcyFluidVelocity classes are only defined for at_nodes = false");
}

void
ComputeFluidVelocityDivergence::initQpStatefulProperties()
{
  _darcy_vel_div_qp[_qp] = 0.0;
}

void
ComputeFluidVelocityDivergence::computeQpProperties()
{
  _darcy_vel_div_qp[_qp] = _darcy_vel_grad[_qp].trace();

  // prepare the derivatives with zeroes
  _ddarcy_vel_div_qp_dvar[_qp].resize(_num_var, RealGradient());
  for (unsigned i = 0; i < _nvel; ++i)
    if (_dictator.isPorousFlowVariable(_vel_var_num[i]))
    {
      // the i_th displacement is a PorousFlow variable
      const unsigned int pvar = _dictator.porousFlowVariableNum(_vel_var_num[i]);
      _ddarcy_vel_div_qp_dvar[_qp][pvar](i) = 1.0;
    }
}