#include "cobylaOptimizerClass.h"

nlopt_opt CobylaOptimizerFunction::GetOptimizer(const unsigned int numberOfParameters)
{
    return nlopt_create(NLOPT_LN_COBYLA, numberOfParameters);
}
