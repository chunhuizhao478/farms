//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/**
 * Created by Chunhui Zhao, Jul 2nd, 2024
 * Implement div v^s term
 * v^s: velocity of solid skeleton
 * Newmark method, see section 2.4.2.1 equation (52) (53) 
 */
#include "ElkADMassConservationNewmark.h"

registerMooseObject("farmsApp", ElkADMassConservationNewmark);

InputParameters ElkADMassConservationNewmark::validParams()
{
    InputParameters params = ADKernel::validParams();
    params.set<bool>("use_displaced_mesh") = false;
    params.addRequiredCoupledVar("displacement_x", "String of displacement component along x direction");
    params.addCoupledVar("displacement_y", 0, "String of displacement component along y direction");
    params.addCoupledVar("displacement_z", 0, "String of displacement component along z direction");
    params.addRequiredCoupledVar("velocity_x", "String of velocity component along x direction");
    params.addCoupledVar("velocity_y", 0, "String of velocity component along x direction");
    params.addCoupledVar("velocity_z", 0, "String of velocity component along x direction");
    params.addRequiredCoupledVar("acceleration_x", "String of acceleration component along x direction");
    params.addCoupledVar("acceleration_y", 0, "String of acceleration component along y direction");
    params.addCoupledVar("acceleration_z", 0, "String of acceleration component along z direction");
    params.addRequiredParam<Real>("beta","beta parameter");
    params.addRequiredParam<Real>("gamma","gamma parameter");
    params.addRequiredParam<bool>("multiply_biot_coefficient","whether or not multiply biot coefficient with the weak form");
    // Plane strain correction (2D simulations): include out-of-plane strain contribution
    params.addParam<bool>("plane_strain_correction", false, "If true (and mesh is 2D), apply epsilon_zz = -nu/(1-nu)*(epsilon_xx+epsilon_yy) to volumetric strain");
    params.addParam<Real>("poissons_ratio", -1.0, "Poisson's ratio used for plane strain correction (required if plane_strain_correction=true)");
    return params;
}

ElkADMassConservationNewmark::ElkADMassConservationNewmark(const InputParameters & parameters):ADKernel(parameters),
    _grad_ux(adCoupledGradient("displacement_x")),
    _grad_uy(adCoupledGradient("displacement_y")),
    _grad_uz(adCoupledGradient("displacement_z")),
    _grad_ux_old(coupledGradientOld("displacement_x")),
    _grad_uy_old(coupledGradientOld("displacement_y")),
    _grad_uz_old(coupledGradientOld("displacement_z")),
    _grad_vx_old(coupledGradientOld("velocity_x")),
    _grad_vy_old(coupledGradientOld("velocity_y")),
    _grad_vz_old(coupledGradientOld("velocity_z")),
    _grad_ax_old(coupledGradientOld("acceleration_x")),
    _grad_ay_old(coupledGradientOld("acceleration_y")),
    _grad_az_old(coupledGradientOld("acceleration_z")),
    _biot_alpha(getADMaterialProperty<Real>("biot_coefficient")),
    _beta(getParam<Real>("beta")),
    _gamma(getParam<Real>("gamma")),
    _multiply_biot_coefficient(getParam<bool>("multiply_biot_coefficient")),
    _biot_modulus(getADMaterialProperty<Real>("biot_modulus")),
    _plane_strain_correction(getParam<bool>("plane_strain_correction")),
    _nu(getParam<Real>("poissons_ratio"))
{}

ADReal
ElkADMassConservationNewmark::computeQpResidual()
{

    //compute divergence of current solid displacement
    ADReal div_u = _grad_ux[_qp](0) + _grad_uy[_qp](1) + _grad_uz[_qp](2);

    //compute divergence of old solid displacement, velocity, acceleration
    ADReal div_u_old = _grad_ux_old[_qp](0) + _grad_uy_old[_qp](1) + _grad_uz_old[_qp](2);
    ADReal div_v_old = _grad_vx_old[_qp](0) + _grad_vy_old[_qp](1) + _grad_vz_old[_qp](2);
    ADReal div_a_old = _grad_ax_old[_qp](0) + _grad_ay_old[_qp](1) + _grad_az_old[_qp](2);

    // If running in 2D with plane strain correction requested, augment divergence
    if (_plane_strain_correction && _mesh.dimension() == 2)
    {
        if (_nu < 0.0 || _nu >= 0.5)
            mooseError("ElkADMassConservationNewmark: provide a valid poissons_ratio (0<=nu<0.5) when plane_strain_correction=true");
        // epsilon_zz = -nu/(1-nu)*(epsilon_xx + epsilon_yy)
        ADReal epszz = -_nu/(1.0 - _nu) * (_grad_ux[_qp](0) + _grad_uy[_qp](1));
        ADReal epszz_old = -_nu/(1.0 - _nu) * (_grad_ux_old[_qp](0) + _grad_uy_old[_qp](1));
        div_u += epszz;
        div_u_old += epszz_old;
        // The old velocity/acceleration divergences already represent time derivatives of divergence; for consistency we
        // recompute their out-of-plane contributions assuming same proportionality (small-strain assumption)
        ADReal div_v_zz_old = -_nu/(1.0 - _nu) * (_grad_vx_old[_qp](0) + _grad_vy_old[_qp](1));
        ADReal div_a_zz_old = -_nu/(1.0 - _nu) * (_grad_ax_old[_qp](0) + _grad_ay_old[_qp](1));
        div_v_old += div_v_zz_old;
        div_a_old += div_a_zz_old;
    }
    
    //compute divergence of current solid acceleration
    //take divergence of both sides of the solid skeleton definition a^s_{n+1} 
    ADReal div_a = 1.0/(_beta*_dt*_dt) * (div_u - div_u_old - div_v_old * _dt - (0.5-_beta) * _dt * _dt * div_a_old);
    
    //compute divergence of current solid velocity
    ADReal div_v = div_v_old + (1-_gamma)*_dt*div_a_old + _gamma*_dt*div_a;
    
    if (_multiply_biot_coefficient){
        return _test[_i][_qp]*_biot_alpha[_qp]*div_v*_biot_modulus[_qp];
    }
    else{
        return _test[_i][_qp]*div_v*_biot_modulus[_qp];
    }
}