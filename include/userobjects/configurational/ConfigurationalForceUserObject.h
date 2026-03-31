//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"
#include <map>

/**
 * Computes configurational forces at nodes using the Eshelby stress tensor.
 *
 * For small deformation: Σ = ΨI - σ
 * Nodal forces: F_CNF = Σ_elements ∫ Σ · ∇N dV
 *
 * Based on: Santarossa et al. (2025) "Configurational forces explain echelon cracks in soft materials"
 */
class ConfigurationalForceUserObject : public ElementUserObject
{
public:
  static InputParameters validParams();

  ConfigurationalForceUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  /**
   * Get a specific component of configurational force at a node
   * @param node_id The node ID
   * @param component The component (0=x, 1=y, 2=z)
   * @return Force component value
   */
  Real getForceComponent(dof_id_type node_id, unsigned int component) const;

  /**
   * Get magnitude of configurational force at a node
   * @param node_id The node ID
   * @return Force magnitude
   */
  Real getForceMagnitude(dof_id_type node_id) const;

  /**
   * Get configurational force vector at a node
   * @param node_id The node ID
   * @return Force vector
   */
  RealVectorValue getForceVector(dof_id_type node_id) const;

protected:
  /// Material property: total energy density Ψ
  const MaterialProperty<Real> & _energy_density;

  /// Material property: stress tensor σ
  const MaterialProperty<RankTwoTensor> & _stress;

  /// Flag to use displacement gradient formulation
  const bool _use_grad_u;

  /// Displacement gradient components (optional)
  const VariableGradient * _grad_disp_x;
  const VariableGradient * _grad_disp_y;
  const VariableGradient * _grad_disp_z;

  /// Optional crack front boundaries for filtering
  std::vector<BoundaryID> _crack_front_boundaries;

  /// Set of crack front node IDs
  std::set<dof_id_type> _crack_front_nodes;

  /// Storage for nodal configurational forces
  std::map<dof_id_type, RealVectorValue> _nodal_forces;

  /**
   * Compute Eshelby stress tensor at quadrature point
   * Small deformation: Σ = ΨI - σ
   * With gradient: Σ = ΨI - (∇u)^T · σ
   */
  RankTwoTensor computeEshelbyStress(unsigned int qp);

  /**
   * Identify crack front nodes from boundary IDs
   */
  void identifyCrackFrontNodes();
};
