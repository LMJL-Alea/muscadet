#ifndef NELDERMEADOPTIMIZERCLASS_H
#define NELDERMEADOPTIMIZERCLASS_H

#include "baseOptimizerClass.h"

class NeldermeadOptimizerFunction : public BaseOptimizerFunction
{
public:
  nlopt_opt GetOptimizer(const unsigned int numberOfParameters);
};

#endif /* NELDERMEADOPTIMIZERCLASS_H */
