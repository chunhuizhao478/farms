
#include "DarcyFluidVelocity.h"
#include "MooseMesh.h"
#include "Assembly.h"

registerMooseObject("farmsApp", DarcyFluidVelocity);

InputParameters
DarcyFluidVelocity::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar(
      "darcy_velocities",
      "The darcy fluid velocity appropriate for the simulation geometry and coordinate system");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

DarcyFluidVelocity::DarcyFluidVelocity(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _nvel(coupledComponents("darcy_velocities")),
    _vel(coupledValues("darcy_velocities")),
    _grad_vel(coupledGradients("darcy_velocities")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _darcy_vel_grad(declareProperty<RankTwoTensor>(_base_name + "darcy_vel_grad")),
    _current_elem_volume(_assembly.elemVolume())
{

  // set unused dimensions to zero
  _vel.resize(3, &_zero);
  _grad_vel.resize(3, &_grad_zero);

  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "The strain calculator needs to run on the undisplaced mesh.");
}

void
DarcyFluidVelocity::initialSetup()
{
  VelocityIntegrityCheck();
}

void
DarcyFluidVelocity::VelocityIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_nvel != _mesh.dimension())
    paramError(
        "darcy_velocities",
        "The number of variables supplied in 'darcy_velocities' must match the mesh dimension.");
}

void
DarcyFluidVelocity::initQpStatefulProperties()
{
  _darcy_vel_grad[_qp].zero();
}