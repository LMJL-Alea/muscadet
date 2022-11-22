#include "neldermeadOptimizerClass.h"

nlopt_opt NeldermeadOptimizerFunction::GetOptimizer(const unsigned int numberOfParameters)
{
    return nlopt_create(NLOPT_LN_NELDERMEAD, numberOfParameters);
}
