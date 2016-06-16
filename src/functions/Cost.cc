#include <iostream>

#include <cube-stacks/functions/Cost.hh>

namespace cubestacks
{
  Cost::Cost()
    : DiffFunction_t(7, 1, "Cost")
  {
  }

  Cost::~Cost()
  {
  }

  void Cost::impl_compute(result_ref res, const_argument_ref x) const
  {
    res[0] = x[2]*x[2];
  }

  void Cost::impl_gradient(gradient_ref grad, const_argument_ref x, size_type) const
  {
    grad[0] = 0;
    grad[1] = 0;
    grad[2] = 2*x[2];
    grad[4] = 0;
    grad[5] = 0;
    grad[6] = 0;
  }
} /* cubestacks */ 
