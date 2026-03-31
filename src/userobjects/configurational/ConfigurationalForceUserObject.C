//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConfigurationalForceUserObject.h"

registerMooseObject("farmsApp", ConfigurationalForceUserObject);

InputParameters
ConfigurationalForceUserObject::validParams()
{
  InputParameters params = ElementUserObject::validParams();

  params.addClassDescription("Computes configurational forces at nodes using the Eshelby stress "
                             "tensor for small deformation fracture mechanics. "
                             "Formula: Sigma = Psi*I - sigma, F_CNF = integral(Sigma dot grad_N)");

  // Required material properties
  params.addRequiredParam<MaterialPropertyName>(
      "energy_density", "Total energy density material property (Psi)");
  params.addRequiredParam<MaterialPropertyName>("stress",
                                                 "Total stress tensor material property (sigma)");

  // Optional displacement variables for more accurate formulation
  params.addParam<bool>("use_displacement_gradient",
                        false,
                        "Use displacement gradient in Eshelby stress: Sigma = Psi*I - (grad_u)^T "
                        "* sigma (more accurate than simple Sigma = Psi*I - sigma)");
  params.addParam<std::vector<VariableName>>(
      "displacements",
      "Displacement variables [disp_x disp_y disp_z] (required if use_displacement_gradient=true)");

  // Optional crack front specification for filtering
  params.addParam<std::vector<BoundaryName>>(
      "crack_front_boundaries",
      "Boundary names defining the crack front (optional - if not specified, computes for all "
      "nodes)");

  return params;
}

ConfigurationalForceUserObject::ConfigurationalForceUserObject(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _energy_density(getMaterialProperty<Real>("energy_density")),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _use_grad_u(getParam<bool>("use_displacement_gradient")),
    _grad_disp_x(nullptr),
    _grad_disp_y(nullptr),
    _grad_disp_z(nullptr)
{
  // Get crack front boundaries if specified
  if (isParamValid("crack_front_boundaries"))
  {
    std::vector<BoundaryName> bnd_names =
        getParam<std::vector<BoundaryName>>("crack_front_boundaries");
    _crack_front_boundaries = _mesh.getBoundaryIDs(bnd_names, true);
  }

  // Setup displacement gradients if using the more accurate formulation
  if (_use_grad_u)
  {
    if (!isParamValid("displacements"))
      mooseError("ConfigurationalForceUserObject: Must specify 'displacements' when "
                 "use_displacement_gradient=true");

    std::vector<VariableName> disp_names = getParam<std::vector<VariableName>>("displacements");
    if (disp_names.size() != 3)
      mooseError("ConfigurationalForceUserObject: Must provide exactly 3 displacement variables "
                 "for 3D problem");

    _grad_disp_x = &coupledGradient("displacements", 0);
    _grad_disp_y = &coupledGradient("displacements", 1);
    _grad_disp_z = &coupledGradient("displacements", 2);
  }
}

void
ConfigurationalForceUserObject::initialize()
{
  // Clear nodal forces from previous timestep
  _nodal_forces.clear();

  // Identify crack front nodes if boundaries specified
  if (!_crack_front_boundaries.empty())
    identifyCrackFrontNodes();
}

void
ConfigurationalForceUserObject::execute()
{
  // Get shape function gradients from assembly
  const auto & grad_phi = _assembly.gradPhi();

  // Loop over all quadrature points in the element
  for (unsigned int qp = 0; qp < _qrule->n_points(); qp++)
  {
    // Compute Eshelby stress tensor: Σ = ΨI - σ (or with ∇u)
    RankTwoTensor Sigma = computeEshelbyStress(qp);

    // Loop over all nodes in the element
    for (unsigned int i = 0; i < _current_elem->n_nodes(); i++)
    {
      // Get node ID
      dof_id_type node_id = _current_elem->node_id(i);

      // Optional: filter to crack front nodes only
      if (!_crack_front_nodes.empty() &&
          _crack_front_nodes.find(node_id) == _crack_front_nodes.end())
        continue;

      // Get shape function gradient ∇N^i at this quadrature point
      RealVectorValue grad_phi_i = grad_phi[i][qp];

      // Compute configurational force contribution: F = Σ · ∇N
      // F_j = Σ_{jk} * (∂N/∂x_k)
      RealVectorValue force_contrib;
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          force_contrib(j) += Sigma(j, k) * grad_phi_i(k);

      // Weight by Jacobian × quadrature weight
      force_contrib *= _JxW[qp] * _coord[qp];

      // Accumulate to nodal force
      _nodal_forces[node_id] += force_contrib;
    }
  }
}

RankTwoTensor
ConfigurationalForceUserObject::computeEshelbyStress(unsigned int qp)
{
  // Get energy density Ψ
  Real psi = _energy_density[qp];

  // Get stress tensor σ
  RankTwoTensor sigma = _stress[qp];

  // Identity tensor
  RankTwoTensor I;
  I.setToIdentity();

  if (!_use_grad_u)
  {
    // Simple small deformation formulation: Σ = ΨI - σ
    return psi * I - sigma;
  }
  else
  {
    // More accurate formulation with displacement gradient: Σ = ΨI - (∇u)^T · σ
    // Build displacement gradient tensor from components
    RankTwoTensor grad_u;
    for (unsigned int i = 0; i < 3; i++)
    {
      grad_u(0, i) = (*_grad_disp_x)[qp](i);
      grad_u(1, i) = (*_grad_disp_y)[qp](i);
      grad_u(2, i) = (*_grad_disp_z)[qp](i);
    }

    return psi * I - grad_u.transpose() * sigma;
  }
}

void
ConfigurationalForceUserObject::identifyCrackFrontNodes()
{
  _crack_front_nodes.clear();

  // Collect all nodes on specified crack front boundaries
  for (auto boundary_id : _crack_front_boundaries)
  {
    ConstBndNodeRange & bnd_nodes = *_mesh.getBoundaryNodeRange();
    for (const auto & bnode : bnd_nodes)
    {
      if (bnode->_bnd_id == boundary_id)
        _crack_front_nodes.insert(bnode->_node->id());
    }
  }
}

void
ConfigurationalForceUserObject::threadJoin(const UserObject & y)
{
  const auto & uo = static_cast<const ConfigurationalForceUserObject &>(y);

  // Merge nodal forces from different threads
  for (const auto & pair : uo._nodal_forces)
    _nodal_forces[pair.first] += pair.second;
}

void
ConfigurationalForceUserObject::finalize()
{
  // Manual parallel sum for map of vectors
  // Collect all unique node IDs across all processors
  std::set<dof_id_type> local_nodes;
  for (const auto & pair : _nodal_forces)
    local_nodes.insert(pair.first);

  // Get union of all node IDs from all processors
  std::set<dof_id_type> all_nodes = local_nodes;
  _communicator.set_union(all_nodes);

  // Sum forces for each node across processors
  std::map<dof_id_type, RealVectorValue> summed_forces;
  for (auto node_id : all_nodes)
  {
    // Get local force (zero if not on this processor)
    RealVectorValue local_force(0, 0, 0);
    auto it = _nodal_forces.find(node_id);
    if (it != _nodal_forces.end())
      local_force = it->second;

    // Sum each component separately
    Real fx = local_force(0);
    Real fy = local_force(1);
    Real fz = local_force(2);

    _communicator.sum(fx);
    _communicator.sum(fy);
    _communicator.sum(fz);

    summed_forces[node_id] = RealVectorValue(fx, fy, fz);
  }

  // Replace local forces with summed forces
  _nodal_forces = summed_forces;
}

// Public access methods for AuxKernels

Real
ConfigurationalForceUserObject::getForceComponent(dof_id_type node_id, unsigned int component) const
{
  mooseAssert(component < 3, "Component must be 0, 1, or 2 for x, y, z");

  auto it = _nodal_forces.find(node_id);
  if (it != _nodal_forces.end())
    return it->second(component);

  return 0.0;
}

Real
ConfigurationalForceUserObject::getForceMagnitude(dof_id_type node_id) const
{
  auto it = _nodal_forces.find(node_id);
  if (it != _nodal_forces.end())
    return it->second.norm();

  return 0.0;
}

RealVectorValue
ConfigurationalForceUserObject::getForceVector(dof_id_type node_id) const
{
  auto it = _nodal_forces.find(node_id);
  if (it != _nodal_forces.end())
    return it->second;

  return RealVectorValue(0, 0, 0);
}
