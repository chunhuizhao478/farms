//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "farmsTestApp.h"
#include "farmsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
farmsTestApp::validParams()
{
  InputParameters params = farmsApp::validParams();
  return params;
}

farmsTestApp::farmsTestApp(InputParameters parameters) : MooseApp(parameters)
{
  farmsTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

farmsTestApp::~farmsTestApp() {}

void
farmsTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  farmsApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"farmsTestApp"});
    Registry::registerActionsTo(af, {"farmsTestApp"});
  }
}

void
farmsTestApp::registerApps()
{
  registerApp(farmsApp);
  registerApp(farmsTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
farmsTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  farmsTestApp::registerAll(f, af, s);
}
extern "C" void
farmsTestApp__registerApps()
{
  farmsTestApp::registerApps();
}
