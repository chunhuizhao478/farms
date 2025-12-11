//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*
Material Description of Slip Weakening Friction 3d
*/

#include "PoroSlipWeakeningFrictionczm3dCDBM.h"
#include "InterfaceKernel.h"
#include "FEProblemBase.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"


registerMooseObject("farmsApp", PoroSlipWeakeningFrictionczm3dCDBM);

InputParameters
PoroSlipWeakeningFrictionczm3dCDBM::validParams()
{
  InputParameters params = CZMComputeLocalTractionTotalBase::validParams();
  params.addClassDescription("linear slip weakening traction separation law.");
  params.addRequiredParam<Real>("mu_s", "value of static friction parameter");
  params.addRequiredParam<Real>("mu_d", "value of dynamic friction parameter");
  params.addRequiredParam<Real>("Dc", "value of characteristic length");
  params.addRequiredParam<Real>("len", "element edge length");
  params.addRequiredCoupledVar("disp_slipweakening_x", "displacement in x dir");
  params.addRequiredCoupledVar("disp_slipweakening_y", "displacement in y dir");
  params.addRequiredCoupledVar("disp_slipweakening_z", "displacement in z dir");
  params.addRequiredCoupledVar("vel_slipweakening_x","velocity in x dir");
  params.addRequiredCoupledVar("vel_slipweakening_y","velocity in y dir");
  params.addRequiredCoupledVar("vel_slipweakening_z","velocity in z dir");
  params.addRequiredCoupledVar("reaction_slipweakening_x", "reaction in x dir");
  params.addRequiredCoupledVar("reaction_slipweakening_y", "reaction in y dir");
  params.addRequiredCoupledVar("reaction_slipweakening_z", "reaction in z dir");
  params.addRequiredCoupledVar("reaction_pressure_x","reaction in x dir");
  params.addRequiredCoupledVar("reaction_pressure_y","reaction in y dir");
  params.addRequiredCoupledVar("reaction_pressure_z","reaction in z dir");
  params.addParam<bool>("use_forced_rupture", false,
                        "use forced rupture or not, default is false");
  params.addParam<Real>("t0", 0.1,
                        "time at which the forced rupture starts, default is 0.1");
  params.addCoupledVar("cohesion_aux", "auxiliary variable for cohesion");
  params.addCoupledVar("forced_rupture_aux", "auxiliary variable for forced rupture");
  params.addRequiredCoupledVar("fault_pressure", "variable for fluid pressure");

  return params;
}

PoroSlipWeakeningFrictionczm3dCDBM::PoroSlipWeakeningFrictionczm3dCDBM(const InputParameters & parameters)
  : CZMComputeLocalTractionTotalBase(parameters),
    _mu_s(getParam<Real>("mu_s")),
    _mu_d(getParam<Real>("mu_d")),
    _Dc(getParam<Real>("Dc")),
    _len(getParam<Real>("len")),
    _density(getMaterialPropertyByName<Real>(_base_name + "density")),
    _rot(getMaterialPropertyByName<RankTwoTensor>(_base_name + "czm_total_rotation")),
    _disp_slipweakening_x(coupledValue("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x(coupledNeighborValue("disp_slipweakening_x")),
    _disp_slipweakening_y(coupledValue("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y(coupledNeighborValue("disp_slipweakening_y")),
    _disp_slipweakening_z(coupledValue("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z(coupledNeighborValue("disp_slipweakening_z")),
    _vel_slipweakening_x(coupledValue("vel_slipweakening_x")),
    _vel_slipweakening_neighbor_x(coupledNeighborValue("vel_slipweakening_x")),
    _vel_slipweakening_y(coupledValue("vel_slipweakening_y")),
    _vel_slipweakening_neighbor_y(coupledNeighborValue("vel_slipweakening_y")),
    _vel_slipweakening_z(coupledValue("vel_slipweakening_z")),
    _vel_slipweakening_neighbor_z(coupledNeighborValue("vel_slipweakening_z")),
    _reaction_slipweakening_x(coupledValue("reaction_slipweakening_x")),
    _reaction_slipweakening_neighbor_x(coupledNeighborValue("reaction_slipweakening_x")),
    _reaction_slipweakening_y(coupledValue("reaction_slipweakening_y")),
    _reaction_slipweakening_neighbor_y(coupledNeighborValue("reaction_slipweakening_y")),
    _reaction_slipweakening_z(coupledValue("reaction_slipweakening_z")),
    _reaction_slipweakening_neighbor_z(coupledNeighborValue("reaction_slipweakening_z")),
    _reaction_pressure_x(coupledValue("reaction_pressure_x")),
    _reaction_pressure_neighbor_x(coupledNeighborValue("reaction_pressure_x")),
    _reaction_pressure_y(coupledValue("reaction_pressure_y")),
    _reaction_pressure_neighbor_y(coupledNeighborValue("reaction_pressure_y")),
    _reaction_pressure_z(coupledValue("reaction_pressure_z")),
    _reaction_pressure_neighbor_z(coupledNeighborValue("reaction_pressure_z")),
    _disp_slipweakening_x_old(coupledValueOld("disp_slipweakening_x")),
    _disp_slipweakening_neighbor_x_old(coupledNeighborValueOld("disp_slipweakening_x")),
    _disp_slipweakening_y_old(coupledValueOld("disp_slipweakening_y")),
    _disp_slipweakening_neighbor_y_old(coupledNeighborValueOld("disp_slipweakening_y")),
    _disp_slipweakening_z_old(coupledValueOld("disp_slipweakening_z")),
    _disp_slipweakening_neighbor_z_old(coupledNeighborValueOld("disp_slipweakening_z")),
    _displacement_jump_strike(declareProperty<Real>("displacement_jump_strike")),
    _displacement_jump_dip(declareProperty<Real>("displacement_jump_dip")),
    _displacement_jump_normal(declareProperty<Real>("displacement_jump_normal")),
    _displacement_jump_rate_strike(declareProperty<Real>("displacement_jump_rate_strike")),
    _displacement_jump_rate_dip(declareProperty<Real>("displacement_jump_rate_dip")),
    _displacement_jump_rate_normal(declareProperty<Real>("displacement_jump_rate_normal")),
    _traction_strike(declareProperty<Real>("traction_strike")),
    _traction_dip(declareProperty<Real>("traction_dip")),
    _traction_normal(declareProperty<Real>("traction_normal")),
    _static_initial_stress_tensor(getMaterialPropertyByName<RankTwoTensor>(_base_name + "sts_initial_tensor")),
    _use_forced_rupture(getParam<bool>("use_forced_rupture")),
    _t0(getParam<Real>("t0")),
    _cohesion_aux(coupledValue("cohesion_aux")),
    _forced_rupture_aux(coupledValue("forced_rupture_aux")),
    _fault_pressure(coupledValue("fault_pressure")),
    _fault_pressure_neighbor(coupledNeighborValue("fault_pressure")),
    _initial_porepressure(getMaterialProperty<Real>("initial_porepressure"))

{

  // only works for small strain
  if (hasBlockMaterialProperty<RankTwoTensor>(_base_name + "strain_increment"))
  {
    mooseError("SlipWeakening only works for small strain!");
  }
}


void
PoroSlipWeakeningFrictionczm3dCDBM::initialSetup()
{
  CZMComputeLocalTractionTotalBase::initialSetup();
  computeNodalVolumePatches();
}

void
PoroSlipWeakeningFrictionczm3dCDBM::computeNodalVolumePatches()
{
  _nodal_volume_patches.clear();
  _nodal_areas.clear();
  
  MeshBase & mesh = _mesh.getMesh();
  
  std::set<dof_id_type> interface_node_ids;
  
  const BoundaryInfo & boundary_info = mesh.get_boundary_info();
  const std::set<BoundaryID> & all_boundary_ids = boundary_info.get_boundary_ids();
  
  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  
  std::map<dof_id_type, Real> interface_face_areas;
  
  for (; el != end_el; ++el)
  {
    const Elem * elem = *el;
    
    for (unsigned int side = 0; side < elem->n_sides(); side++)
    {
      for (auto bid : all_boundary_ids)
      {
        if (boundary_info.has_boundary_id(elem, side, bid))
        {
          const Elem * neighbor = elem->neighbor_ptr(side);
          
          if (neighbor == nullptr)
            continue;
          
          std::unique_ptr<const Elem> side_elem = elem->build_side_ptr(side);
          Real face_area = side_elem->volume();
          unsigned int n_face_nodes = side_elem->n_nodes();
          
          for (unsigned int n = 0; n < n_face_nodes; n++)
          {
            dof_id_type node_id = side_elem->node_id(n);
            interface_node_ids.insert(node_id);
            interface_face_areas[node_id] += face_area / static_cast<Real>(n_face_nodes);
          }
          
          break;
        }
      }
    }
  }
  
  el = mesh.active_elements_begin();
  for (; el != end_el; ++el)
  {
    const Elem * elem = *el;
    
    if (elem->dim() != 3)
      continue;
    
    bool has_interface_node = false;
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
    {
      if (interface_node_ids.count(elem->node_id(n)) > 0)
      {
        has_interface_node = true;
        break;
      }
    }
    
    if (!has_interface_node)
      continue;
    
    Real elem_volume = elem->volume();
    unsigned int n_nodes = elem->n_nodes();
    
    for (unsigned int n = 0; n < n_nodes; n++)
    {
      dof_id_type node_id = elem->node_id(n);
      
      if (interface_node_ids.count(node_id) > 0)
      {
        _nodal_volume_patches[node_id] += elem_volume / static_cast<Real>(n_nodes);
      }
    }
  }
  
  for (const auto & pair : interface_face_areas)
  {
    _nodal_areas[pair.first] = pair.second;
  }

  // // -----------------------------------------------------
  // // (4) PRINT all results
  // // -----------------------------------------------------
  // Moose::out << "\n=== NODAL VOLUME PATCHES & AREAS ===\n";

  // for (const auto & kv : _nodal_volume_patches)
  // {
  //   dof_id_type node_id = kv.first;
  //   Real V = kv.second;
  //   Real A = _nodal_areas.count(node_id) ? _nodal_areas[node_id] : 0.0;

  //   Moose::out << "Node " << node_id
  //              << " | Volume = " << V
  //              << " | Area = " << A
  //              << std::endl;
  // }

  // Moose::out << "=== END PRINT ===\n\n";
}

void
PoroSlipWeakeningFrictionczm3dCDBM::computeInterfaceTractionAndDerivatives()
{
  if(_fe_problem.getCurrentExecuteOnFlag()=="LINEAR")
  {
  // Global Displacement Jump
  RealVectorValue displacement_jump_global(
      _disp_slipweakening_x[_qp] - _disp_slipweakening_neighbor_x[_qp],
      _disp_slipweakening_y[_qp] - _disp_slipweakening_neighbor_y[_qp],
      _disp_slipweakening_z[_qp] - _disp_slipweakening_neighbor_z[_qp]);
  
  // Global Displacement Jump Old
  RealVectorValue displacement_jump_old_global(
    _disp_slipweakening_x_old[_qp] - _disp_slipweakening_neighbor_x_old[_qp],
    _disp_slipweakening_y_old[_qp] - _disp_slipweakening_neighbor_y_old[_qp],
    _disp_slipweakening_z_old[_qp] - _disp_slipweakening_neighbor_z_old[_qp]);

  // Global Displacement Jump Rate
  RealVectorValue displacement_jump_rate_global =
      (displacement_jump_global - displacement_jump_old_global) * (1 / _dt);

  // Local Displacement Jump / Displacement Jump Rate
  RealVectorValue displacement_jump = _rot[_qp].transpose() * displacement_jump_global;
  RealVectorValue displacement_jump_rate = _rot[_qp].transpose() * displacement_jump_rate_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real displacement_jump_n = displacement_jump(0);
  Real displacement_jump_t = displacement_jump(1);
  Real displacement_jump_d = displacement_jump(2);
  Real displacement_jump_rate_n = displacement_jump_rate(0);
  Real displacement_jump_rate_t = displacement_jump_rate(1);
  Real displacement_jump_rate_d = displacement_jump_rate(2);

  // Parameter initialization
  Real tau_f = 0;

  // Reaction force in local coordinate
  RealVectorValue R_plus_global(-_reaction_slipweakening_x[_qp],
                                -_reaction_slipweakening_y[_qp],
                                -_reaction_slipweakening_z[_qp]);
  RealVectorValue R_minus_global(-_reaction_slipweakening_neighbor_x[_qp],
                                 -_reaction_slipweakening_neighbor_y[_qp],
                                 -_reaction_slipweakening_neighbor_z[_qp]);

  RealVectorValue R_plus_local = _rot[_qp].transpose() * R_plus_global;
  RealVectorValue R_minus_local = _rot[_qp].transpose() * R_minus_global;

 // Reaction pressure in local coordinate
  RealVectorValue R_plus_pressure_global(-_reaction_pressure_x[_qp],
                                         -_reaction_pressure_y[_qp],
                                        -_reaction_pressure_z[_qp]);
  RealVectorValue R_minus_pressure_global(-_reaction_pressure_neighbor_x[_qp],
                                          -_reaction_pressure_neighbor_y[_qp],
                                          -_reaction_pressure_neighbor_z[_qp]);

  RealVectorValue R_plus_pressure_local = _rot[_qp].transpose() * R_plus_pressure_global;
  RealVectorValue R_minus_pressure_local = _rot[_qp].transpose() * R_minus_pressure_global;

  // n is along normal direction; t is along tangential direction; d is along dip direction
  Real R_plus_local_n = R_plus_local(0);
  Real R_plus_local_t = R_plus_local(1);
  Real R_plus_local_d = R_plus_local(2);
  Real R_minus_local_n = R_minus_local(0);
  Real R_minus_local_t = R_minus_local(1);
  Real R_minus_local_d = R_minus_local(2);

  Real R_plus_pressure_local_n = R_plus_pressure_local(0);
  Real R_plus_pressure_local_t = R_plus_pressure_local(1);
  Real R_plus_pressure_local_d = R_plus_pressure_local(2);
  Real R_minus_pressure_local_n = R_minus_pressure_local(0);
  Real R_minus_pressure_local_t = R_minus_pressure_local(1);
  Real R_minus_pressure_local_d = R_minus_pressure_local(2);

// // ===== FIXED SECTION: Compute M and A for node nearest to quadrature point =====
//   // Find the node closest to the current quadrature point that has volume patch data
//   Real M = 0.0;
//   Real A = 0.0;
//   bool found = false;
  
//   // Get the physical coordinates of the current quadrature point
//   const Point & qp_coords = _q_point[_qp];
  
//   Real min_distance = std::numeric_limits<Real>::max();
//   dof_id_type nearest_node_id = 0;
  
//   unsigned int n_elem_nodes = _current_elem->n_nodes();
//   for (unsigned int i = 0; i < n_elem_nodes; i++)
//   {
//     const Node * node = _current_elem->node_ptr(i);
//     dof_id_type node_id = node->id();
    
//     // Check if this node has volume patch data
//     auto vol_it = _nodal_volume_patches.find(node_id);
//     auto area_it = _nodal_areas.find(node_id);
    
//     if (vol_it != _nodal_volume_patches.end() && area_it != _nodal_areas.end())
//     {
//       // Calculate distance from quadrature point to this node
//       Real distance = (*node - qp_coords).norm();
      
//       if (distance < min_distance)
//       {
//         min_distance = distance;
//         nearest_node_id = node_id;
//         M = _density[_qp] * vol_it->second;  // M = rho * V
//         A = area_it->second;                  // A
//         found = true;
//       }
//     }
//   }

  Real M = 0.0;
  Real A = 0.0;
  
  M = _density[_qp] * _len * _len * _len / 2;
  A = _len * _len ;
  
  // ===== END FIXED SECTION =====

  // Compute T1_o, T2_o, T3_o for current qp
  //!!! rotation matrix is not applied here !!!
  Real T1_o = _static_initial_stress_tensor[_qp](0, 1); // shear stress in t dir
  Real T2_o = -1.0 * _static_initial_stress_tensor[_qp](1, 1); // normal stress in n dir
  Real T3_o = _static_initial_stress_tensor[_qp](0, 2); // shear stress in d dir

  // Compute sticking stress
  Real T1 = (1 / _dt) * M * displacement_jump_rate_t / (2 * A) +
            (R_plus_local_t + R_plus_pressure_local_t - R_minus_local_t - R_minus_pressure_local_t) / (2 * A) + T1_o;
  Real T3 = (1 / _dt) * M * displacement_jump_rate_d / (2 * A) +
            (R_plus_local_d + R_plus_pressure_local_d- R_minus_local_d - R_minus_pressure_local_d) / (2 * A) + T3_o;
  // Real T2 = -(1 / _dt) * M * (displacement_jump_rate_n + (1 / _dt) * displacement_jump_n) /
  //               (2 * A) +
  //           ((R_minus_local_n + R_minus_pressure_local_n - R_plus_local_n - R_plus_pressure_local_n) / (2 * A)) - T2_o;
  
  Real T2 = ((R_minus_local_n + R_minus_pressure_local_n - R_plus_local_n - R_plus_pressure_local_n) / (2 * A)) - T2_o;

  Real pressure_primary = _fault_pressure[_qp];
  Real pressure_neighbor = _fault_pressure_neighbor[_qp];

  Real Pmax = std::max(pressure_primary, pressure_neighbor);

  Real Pf = _initial_porepressure[_qp] + Pmax; // fault pressure

  //T2: total normal stress acting on the fault, taken to be "positive" in compression: -T2
  //treat tension on the fault the same as if the effective normal stress equals zero.
  Real effective_stress = (-T2) - Pf;


  // Overstress nucleation
  if (!_use_forced_rupture){

    //tau_f
    //T2: total normal stress acting on the fault, taken to be "positive" in compression: -T2
    //treat tension on the fault the same as if the effective normal stress equals zero.


    // Compute fault traction
    if (effective_stress < 0)
    {
    }
    else
    {
      effective_stress = 0;
    }

    // Compute friction strength
    Real slip_total = std::sqrt(displacement_jump_t*displacement_jump_t+displacement_jump_d*displacement_jump_d);
    if (slip_total < _Dc)
    {
      tau_f = (_mu_s - (_mu_s - _mu_d) * slip_total / _Dc) *
              (-effective_stress); // square for shear component
    }
    else
    {
      tau_f = _mu_d * (-effective_stress);
    }
  
  }
  //Forced rupture nucleation
  else{

    //parameter f1
    Real f1 = 0.0;
    Real slip_total = std::sqrt(displacement_jump_t*displacement_jump_t+displacement_jump_d*displacement_jump_d);
    if ( slip_total < _Dc ){
      f1 = ( 1.0 * slip_total ) / ( 1.0 * _Dc );
    }
    else{
      f1 = 1;
    }

    //parameter f2
    //here we close the gradual reduction on mud, replace it by overstress
    Real f2 = 0.0;
    Real t0 = _t0;
    Real T = _forced_rupture_aux[_qp];
    if ( _t < T ){
      f2 = 0.0;
    }
    else if ( _t > T && _t < T + t0 ){
      f2 = ( _t - T ) / t0;
    }
    else{
      f2 = 1;
    }

    Real mu = _mu_s + ( _mu_d - _mu_s ) * std::max(f1,f2);

    tau_f = _cohesion_aux[_qp] + mu * std::max(effective_stress,0.0);

  }

  // Compute fault traction
  if (std::sqrt(T1 * T1 + T3 * T3) < tau_f)
  {
  }
  else
  {
    T1 = tau_f * T1 / std::sqrt(T1 * T1 + T3 * T3);
    T3 = tau_f * T3 / std::sqrt(T1 * T1 + T3 * T3);
  }


  // Save displacement jump in local coordinate
  _displacement_jump_strike[_qp] = displacement_jump_t; // strike direction
  _displacement_jump_dip[_qp] = displacement_jump_d; // dip direction
  _displacement_jump_normal[_qp] = displacement_jump_n; // normal direction

  // Save displacement jump rate in local coordinate
  _displacement_jump_rate_strike[_qp] = displacement_jump_rate_t; // strike direction
  _displacement_jump_rate_dip[_qp] = displacement_jump_rate_d; // dip direction
  _displacement_jump_rate_normal[_qp] = displacement_jump_rate_n; // normal direction

  // Save traction in local coordinate
  _traction_strike[_qp] = T1; // strike direction
  _traction_normal[_qp] = T2; // normal direction
  _traction_dip[_qp] = T3;  // dip direction  

  // Assign back traction in CZM
  RealVectorValue traction(T2 + T2_o, -T1 + T1_o, -T3 + T3_o);
  //RealVectorValue traction(0.0,0.0,0.0);
  _interface_traction[_qp] = traction;
  _dinterface_traction_djump[_qp] = 0;
}
}