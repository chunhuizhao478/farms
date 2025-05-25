#include "InitialDamageCycleSim3DPlaneFunction.h"

registerMooseObject("farmsApp", InitialDamageCycleSim3DPlaneFunction);

InputParameters
InitialDamageCycleSim3DPlaneFunction::validParams()
{
  InputParameters params = Function::validParams();
  params.addClassDescription("Function for an initial damage profile in 3D with an "
                             "exponential decay away from a rectangular fault plane in x-z.");
  params.addRequiredParam<Real>("len_of_fault_strike",
                                "Length of the fault in the x-direction");
  params.addRequiredParam<Real>("len_of_fault_dip",
                                "Length of the fault in the z-direction (the 'dip' extent)");
  params.addRequiredParam<Real>("sigma",
                                "Decay rate (controls radial spread of damage)");
  params.addRequiredParam<Real>("peak_val",
                                "Peak value of the initial damage at the plane");
  params.addRequiredParam<std::vector<Real>>("nucl_center", "nucleation center (x,y,z)");
  return params;
}

InitialDamageCycleSim3DPlaneFunction::InitialDamageCycleSim3DPlaneFunction(const InputParameters & parameters)
  : Function(parameters),
    _len_of_fault_strike(getParam<Real>("len_of_fault_strike")),
    _len_of_fault_dip(getParam<Real>("len_of_fault_dip")),
    _sigma(getParam<Real>("sigma")),
    _peak_val(getParam<Real>("peak_val")),
    _nucl_center(getParam<std::vector<Real>>("nucl_center"))
{
}

Real
InitialDamageCycleSim3DPlaneFunction::value(Real /*t*/, const Point & p) const
{
    // Coordinates of the current quadrature point
    const Real x_coord = p(0); //along the dip direction
    const Real y_coord = p(1); //along the dip direction
    const Real z_coord = p(2); //along the dip direction

    // 1) Clamp x to the fault rectangle in the x-direction
    const Real x_min = _nucl_center[0] - 0.5 * _len_of_fault_strike;
    const Real x_max = _nucl_center[0] + 0.5 * _len_of_fault_strike;
    const Real x_clamped = std::max(x_min, std::min(x_coord, x_max));;

    // 2) Clamp z to the fault rectangle in the z-direction
    const Real z_min = _nucl_center[2] - 0.5 * _len_of_fault_dip;
    const Real z_max = _nucl_center[2] + 0.5 * _len_of_fault_dip;
    const Real z_clamped = std::max(z_min, std::min(z_coord, z_max));

    // The plane is at y = 0, so the closest point on the plane to (x, y, z)
    // is (x_clamped, 0, z_clamped). Compute the distance from that point.
    const Real dx = x_coord - x_clamped;
    const Real dy = y_coord - _nucl_center[1];           // because plane is at y=0
    const Real dz = z_coord - z_clamped;
    const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);

    // Exponential decay
    Real alpha_o = _peak_val * std::exp(-1.0 * (r * r) / (_sigma * _sigma));

    // Final initial-damage value
    return alpha_o;

}