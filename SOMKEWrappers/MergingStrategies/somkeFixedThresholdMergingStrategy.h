//
// Created by tomev on 12/03/2021.
//

#ifndef KERDEP_SOMKEFIXEDTHRESHOLDMERGINGSTRATEGY_H
#define KERDEP_SOMKEFIXEDTHRESHOLDMERGINGSTRATEGY_H

#include "SOMKE/include/SOMKE/SOMKEMergingStrategy.h"

class SOMKEFixedThresholdMergingStrategy : public SOMKEMergingStrategy {
  public:
    explicit SOMKEFixedThresholdMergingStrategy(const double &threshold = 1,
                                                const double &beta = 0);
    virtual bool ShouldMergeBePerformed(const vector<SOMSequenceEntry> &entries) const override;
    virtual size_t FindEntryToMerge(const vector<SOMSequenceEntry> &entries) const override;
    virtual void SetDataWindowIterator(int *iterator) override;
  protected:
    double threshold_;
    double beta_;
    int *data_window_iterator_ = nullptr;

    double ComputeModifiedDivergenceOfSOMSequenceEntry(const SOMSequenceEntry &entry) const;
};

#endif //KERDEP_SOMKEFIXEDTHRESHOLDMERGINGSTRATEGY_H
