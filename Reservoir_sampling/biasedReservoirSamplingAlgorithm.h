//
// Created by Tomev on 03.06.2017.
//

#ifndef RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H
#define RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H

#include "dataReader.h"
#include "dataParser.h"

#include <random>

/*
 *      A class for biased reservoir sampling as proposed in:
 *      On biased reservoir sampling in the presence of stream evolution by Aggarwal C..
 */

class biasedReservoirSamplingAlgorithm {

    public:

        biasedReservoirSamplingAlgorithm(dataParser *parser, dataReader *reader);
        void fillReservoir(void *reservoir);

    private:

    // Bias rate. According to Aggarwal should be in [0,1]
    double  biasRate = 0.001;
    // Fraction of reservoir filled is also a success rate of deleting a reservoir member.
    double fractionOfReservoirFilled = 0.0;

    dataReader *reader;
    dataParser *parser;

    const int STEPS_NUMBER = 20000;
    //const int RESERVOIR_SIZE = 1000;

    void updateFractionOfReservoirFilled(int currentReservoirSize);
};


#endif //RESERVOIRALGORITHM_BIASEDRESERVOIRSAMPLINGALGORITHM_H
