//
// Created by tomev on 12/03/2021.
//

#ifndef DEDSTA_SOMKEFIXEDMEMORYMERGINGSTRATEGY_H
#define DEDSTA_SOMKEFIXEDMEMORYMERGINGSTRATEGY_H

#include "SOMKE/include/SOMKE/SOMKEMergingStrategy.h"

class SOMKEFixedMemoryMergingStrategy : public SOMKEMergingStrategy {
  public:
    explicit SOMKEFixedMemoryMergingStrategy(const size_t &maximal_number_of_entries_in_sequence = 1,
                                             const double &beta = 0);
    virtual bool ShouldMergeBePerformed(const vector<SOMSequenceEntry> &entries) const override;
    virtual size_t FindEntryToMerge(const vector<SOMSequenceEntry> &entries) const override;
    virtual void SetDataWindowIterator(int *iterator) override;

  protected:

    size_t maximal_number_of_entries_in_sequence_;
    double beta_ = 0; // Rate of increased probability of older entries merge
    int *data_window_iterator_ = nullptr;

    double ComputeModifiedDivergenceOfSOMSequenceEntry(const SOMSequenceEntry &entry) const;


};

#endif //DEDSTA_SOMKEFIXEDMEMORYMERGINGSTRATEGY_H
