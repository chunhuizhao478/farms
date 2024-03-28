#include "InitialAlphaAD.h"

registerMooseObject("farmsApp", InitialAlphaAD);

InputParameters
InitialAlphaAD::validParams()
{
  InputParameters params = Function::validParams();
  return params;
}

InitialAlphaAD::InitialAlphaAD(const InputParameters & parameters)
  : Function(parameters)
{
}

Real
InitialAlphaAD::value(Real /*t*/, const Point & p) const
{


  Real x_coord = p(0); //along the strike direction
  Real y_coord = p(1); //along the normal direction
  Real z_coord = p(2); //along the dip direction  

  Real alpha_o = 0;

  if (z_coord >= 0-0.01 and z_coord <= 0){
    if (x_coord >= -0.6-2*0.01 and x_coord <= -0.6+2*0.01 and y_coord >= 0-2*0.01 and y_coord <= 0+2*0.01){
      alpha_o = 0.8;
    }    
    else if (y_coord >= -0.10 and y_coord <= 0.10 and x_coord >= -0.70 and x_coord <= 0.70){
      alpha_o = 0.7;
    }
    else if (y_coord >= -0.11 and y_coord <= 0.11 and x_coord >= -0.71 and x_coord <= 0.71){
      alpha_o = 0.6;
    }
    else if (y_coord >= -0.12 and y_coord <= 0.12 and x_coord >= -0.72 and x_coord <= 0.72){
      alpha_o = 0.5;
    }
    else if (y_coord >= -0.13 and y_coord <= 0.13 and x_coord >= -0.73 and x_coord <= 0.73){
      alpha_o = 0.4;
    }
    else if (y_coord >= -0.14 and y_coord <= 0.14 and x_coord >= -0.74 and x_coord <= 0.74){
      alpha_o = 0.3;
    }     
    else if (y_coord >= -0.15 and y_coord <= 0.15 and x_coord >= -0.75 and x_coord <= 0.75){
      alpha_o = 0.2;
    } 
    else if (y_coord >= -0.16 and y_coord <= 0.16 and x_coord >= -0.76 and x_coord <= 0.76){
      alpha_o = 0.1;
    }
  }
  else{
    alpha_o = 0.0;
  }
  
  return alpha_o;

}