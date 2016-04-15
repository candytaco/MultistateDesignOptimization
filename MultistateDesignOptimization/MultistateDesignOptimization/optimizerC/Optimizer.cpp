#include "Optimizer.h"

using namespace OPTIMIZER;

Optimizer::Optimizer()
{
	Optimizer::Optimizer(0, NULL, false);
}

Optimizer::Optimizer(const Optimizer &existing)
{
	// TODO: implement this
}

Optimizer::Optimizer(int nMacrostates, string *macrostates)
{
	Optimizer::Optimizer(nMacrostates, macrostates, false);
}

Optimizer::Optimizer(int nMacrostates, bool continuousBoltzmann)
{
	Optimizer::Optimizer(nMacrostates, NULL, continuousBoltzmann);
}

Optimizer::Optimizer(int nMacrostates, string *macrostates, bool continuousBoltzmann)
{
	// TODO: implement this
}