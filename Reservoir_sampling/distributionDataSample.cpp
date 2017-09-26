#include "distributionDataSample.h"
#include <iostream>
#include <string>

void distributionDataSample::print()
{
  std::cout << "Sample: " << std::endl;

  for(auto nameValue : attributesValues)
    std::cout << nameValue.first << " = " << nameValue.second << std::endl;
}
