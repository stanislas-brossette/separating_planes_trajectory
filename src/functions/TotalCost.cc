#include <iostream>

#include <cube-stacks/functions/TotalCost.hh>

namespace cubestacks
{
TotalCost::TotalCost(const IndexManager& indexMngr)
    : DiffFunction_t(indexMngr.totalDim(), 1, "TotalCost"), nCubes_(indexMngr.nCubes())
{
}

TotalCost::~TotalCost() {}

void TotalCost::impl_compute(result_ref res, const_argument_ref x) const
{
  res[0] = 0;
  for (int i = 0; i < nCubes_; i++)
  {
    res[0] += x[7 * i + 2];
  }
}

void TotalCost::impl_gradient(gradient_ref grad, const_argument_ref,
                              size_type) const
{
  for (int i = 0; i < nCubes_; i++)
  {
    grad[7 * i + 0] = 0;
    grad[7 * i + 1] = 0;
    grad[7 * i + 2] = 1;
    grad[7 * i + 4] = 0;
    grad[7 * i + 5] = 0;
    grad[7 * i + 6] = 0;
  }
}
} /* cubestacks */
