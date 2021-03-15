//
// Created by tomev on 10/03/2021.
//

#include "somkeNormalKernel.h"

double SOMKENormalKernel::GetValue(const Point &pt) {
  auto x = pt;
  return normalKernel::getValue(&x);
}
