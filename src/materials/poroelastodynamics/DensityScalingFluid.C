

#include "DensityScalingFluid.h"
#include "libmesh/utility.h"

registerMooseObject("farmsApp", DensityScalingFluid);

InputParameters
DensityScalingFluid::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>(
      "density_fluid",
      "Name of Material Property or a constant real number defining the fluid density of the material.");
  params.addRequiredParam<Real>("desired_time_step", "Time step to achieve.");
  params.addParam<Real>(
      "factor",
      1.0,
      "Factor to multiply to the critical time step. This factor is typically less than one to be "
      "on the conservative side due to two types of approximation: Time step lagging and the "
      "approximation of the critical time step formula.");
  return params;
}

DensityScalingFluid::DensityScalingFluid(const InputParameters & parameters)
  : Material(parameters),
    _desired_time_step(getParam<Real>("desired_time_step")),
    _density_scaling_fluid(declareProperty<Real>("density_scaling_fluid")),
    _material_density_fluid(getMaterialPropertyByName<Real>("density_fluid")),
    _effective_stiffness(getMaterialPropertyByName<Real>("effective_stiffness")),
    _factor(getParam<Real>("factor"))
{
  mooseInfo("Since it can change key simulation results, usage of selective density (mass) scaling "
            "is only recommended for advanced users.");
}

void
DensityScalingFluid::computeQpProperties()
{
  const Real critical = _factor * _current_elem->hmin() * std::sqrt(_material_density_fluid[_qp]) /
                        (_effective_stiffness[0]);

  if (critical < _desired_time_step)
  {
    const Real desired_density = std::pow(_effective_stiffness[_qp] * _desired_time_step, 2) /
                                 std::pow(_factor * _current_elem->hmin(), 2);

    const Real density_to_add =
        desired_density > _material_density_fluid[_qp] ? desired_density - _material_density_fluid[_qp] : 0.0;

    _density_scaling_fluid[_qp] = density_to_add;
  }
  else
    _density_scaling_fluid[_qp] = 0.0;
}