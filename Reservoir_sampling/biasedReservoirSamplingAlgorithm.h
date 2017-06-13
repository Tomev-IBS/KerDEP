//
// Created by Tomev on 03.06.2017.
//

#ifndef RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H

#include "dataReader.h"
#include "dataParser.h"
#include "reservoirSamplingAlgorithm.h"

#include <random>

/*
 *      A class for biased reservoir sampling as proposed in:
 *      On biased reservoir sampling in the presence of stream evolution by Aggarwal C..
 */

class biasedReservoirSamplingAlgorithm : public reservoirSamplingAlgorithm
{
  public:

  biasedReservoirSamplingAlgorithm(dataReader *reader, dataParser *parser, int reservoirMaxSize, int stepsNumber);

  void fillReservoir(std::vector<sample*> *reservoir);
  void performSingleStep(std::vector<sample*> *reservoir, int stepNumber);

  private:

  // Bias rate. According to Aggarwal should be in [0,1]
  double  biasRate = 0.001;
  // Fraction of reservoir filled is also a success rate of deleting a reservoir member.
  double fractionOfReservoirFilled = 0.0;

  void updateFractionOfReservoirFilled(int currentReservoirSize);
};


#endif //RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H
