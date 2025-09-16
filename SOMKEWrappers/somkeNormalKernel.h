//
// Created by tomev on 10/03/2021.
//

#ifndef DEDSTA_SOMKENORMALKERNEL_H
#define DEDSTA_SOMKENORMALKERNEL_H

#include "Kernel.h"
#include "Functions/Kernels/normalkernel.h"

class SOMKENormalKernel : public Kernel, normalKernel {
  public:
    double GetValue(const Point &pt) override;
};

#endif //DEDSTA_SOMKENORMALKERNEL_H
