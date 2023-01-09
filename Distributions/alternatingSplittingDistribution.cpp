//
// Created by tomev on 1/5/2023.
//

#include "alternatingSplittingDistribution.h"

alternatingSplittingDistribution::alternatingSplittingDistribution(int seed,
                                                                   vector<std::shared_ptr<distribution>> *elementalDistributions,
                                                                   vector<double> *contributions, double addendMultiplier) : complexDistribution(seed,
                                                                                                              elementalDistributions,
                                                                                                              contributions), _addendMultiplier(addendMultiplier)
{}

void alternatingSplittingDistribution::increaseMeans(double addend, int index) {

    for(int i = 0; i < elementalDistributions.size(); ++i){
      if(i % 2 ==0){
        elementalDistributions[i]->increaseMeans(addend * _addendMultiplier, index);
      }
      else {
        elementalDistributions[i]->increaseMeans(addend, index);
      }
    }

}

int alternatingSplittingDistribution::randomizeDistributionIndex() {
  return complexDistribution::randomizeDistributionIndex();
  // return _generatedSamples++ % elementalDistributions.size();  // For alternation
}
