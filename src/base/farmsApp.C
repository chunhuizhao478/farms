#include "farmsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
farmsApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

registerKnownLabel("farmsApp");

farmsApp::farmsApp(InputParameters parameters) : MooseApp(parameters)
{
  farmsApp::registerAll(_factory, _action_factory, _syntax);
}

farmsApp::~farmsApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
registerSyntax("PoroCohesiveZoneAction", "Actions/PoroCohesiveZoneAction/*");
}


void
farmsApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<farmsApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"farmsApp"});
  Registry::registerActionsTo(af, {"farmsApp"});
  associateSyntaxInner(syntax, af);
  registerDataFilePath();

  /* register custom execute flags, action syntax, etc. here */
}


void
farmsApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"farmsApp"});
  associateSyntaxInner(syntax, action_factory);
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
