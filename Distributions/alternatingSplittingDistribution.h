//
// Created by tomev on 1/5/2023.
//
// This class multiplies added from even numbered distributions in increaseMeans method by a given number.
// It also ignores the contributions and samples from the elementary distributions in the order they were given.
//

#ifndef DEDSTA_ALTERNATINGSPLITTINGDISTRIBUTION_H
#define DEDSTA_ALTERNATINGSPLITTINGDISTRIBUTION_H

#include "complexdistribution.h"

class alternatingSplittingDistribution : public complexDistribution{

  public:

    alternatingSplittingDistribution(int seed, vector<std::shared_ptr<distribution>>* elementalDistributions, vector<double>* contributions, double addendMultiplier = 1);
    virtual void increaseMeans(double addend, int index=-1);

  protected:
    double _addendMultiplier;
    int _generatedSamples = 0;

    int randomizeDistributionIndex();

};

#endif //DEDSTA_ALTERNATINGSPLITTINGDISTRIBUTION_H
