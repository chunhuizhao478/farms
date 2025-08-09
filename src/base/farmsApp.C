#include "farmsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
farmsApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

farmsApp::farmsApp(InputParameters parameters) : MooseApp(parameters)
{
  farmsApp::registerAll(_factory, _action_factory, _syntax);
}

farmsApp::~farmsApp() {}

void
farmsApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<farmsApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"farmsApp"});
  Registry::registerActionsTo(af, {"farmsApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
farmsApp::registerApps()
{
  registerApp(farmsApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
farmsApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  farmsApp::registerAll(f, af, s);
}
extern "C" void
farmsApp__registerApps()
{
  farmsApp::registerApps();
}
