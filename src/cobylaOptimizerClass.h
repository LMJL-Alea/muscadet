#ifndef COBYLAOPTIMIZERCLASS_H
#define COBYLAOPTIMIZERCLASS_H

#include "baseOptimizerClass.h"

class CobylaOptimizerFunction : public BaseOptimizerFunction
{
public:
  nlopt_opt GetOptimizer(const unsigned int numberOfParameters);
};

#endif /* COBYLAOPTIMIZERCLASS_H */
