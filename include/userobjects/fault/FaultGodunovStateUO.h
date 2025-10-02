#pragma once
#include "InterfaceUserObject.h"
#include "ElasticDGUtils3D.h"
#include <unordered_map>
#include <vector>
#include <limits>

/**
 * FaultGodunovStateUO
 * Computes and stores Godunov elastic interface states and imposed frictional tractions
 * on a specified boundary (fault) each time step so that DGKernels and AuxKernels can reuse
 * them without recomputation.
 */
class FaultGodunovStateUO : public InterfaceUserObject
{
public:
  static InputParameters validParams();
  FaultGodunovStateUO(const InputParameters & params);

  void initialize() override;
  void execute() override;   // called for each interface face
  void finalize() override {}
  void threadJoin(const UserObject & y) override;

  struct State
  {
    Real snnG, snt1G, snt2G;
    Real unG, ut1_star, ut2_star;
    Real tauimp_t1, tauimp_t2;
  Real sigmaN; // snnG + sigma0n
  Real sr_t1, sr_t2; // slip-rate components (minus - plus) after radiation damping adjustment
  };

  // Access by (elem id, side, qp). Throws if missing.
  const State & get(dof_id_type elem_id, unsigned short side, unsigned short qp) const;
  bool has(dof_id_type elem_id, unsigned short side, unsigned short qp) const;

protected:
  // Parameters / names
  const MaterialProperty<Real> & _lambda;  const MaterialProperty<Real> & _lambda_n;
  const MaterialProperty<Real> & _mu;      const MaterialProperty<Real> & _mu_n;
  const MaterialProperty<Real> & _rho;     const MaterialProperty<Real> & _rho_n;

  // Coupled fields
  const VariableValue & _ux; const VariableValue & _ux_n;
  const VariableValue & _uy; const VariableValue & _uy_n;
  const VariableValue & _uz; const VariableValue & _uz_n;
  const VariableValue & _sxx; const VariableValue & _sxx_n;
  const VariableValue & _syy; const VariableValue & _syy_n;
  const VariableValue & _szz; const VariableValue & _szz_n;
  const VariableValue & _sxy; const VariableValue & _sxy_n;
  const VariableValue & _sxz; const VariableValue & _sxz_n;
  const VariableValue & _syz; const VariableValue & _syz_n;
  const VariableValue & _slip_old;

  // Friction params
  const Real _mu_s, _mu_d, _Dc; const Real _tau0_t1, _tau0_t2, _sigma0n;

  inline Real mu_slipweak(Real d) const
  { if (d < _Dc) return std::max(_mu_d, _mu_s - (_mu_s - _mu_d) * (d/_Dc)); return _mu_d; }

  // Storage
  std::vector<State> _states;
  std::unordered_map<unsigned long long, std::size_t> _lookup;

  static unsigned long long pack(dof_id_type e, unsigned short side, unsigned short qp)
  { return ( (unsigned long long)e << 24 ) | ( (unsigned long long)side << 12 ) | (unsigned long long)qp; }
};
