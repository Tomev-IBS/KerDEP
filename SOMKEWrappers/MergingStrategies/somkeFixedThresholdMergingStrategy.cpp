//
// Created by tomev on 12/03/2021.
//

#include "somkeFixedThresholdMergingStrategy.h"

SOMKEFixedThresholdMergingStrategy::SOMKEFixedThresholdMergingStrategy(const double &threshold, const double &beta) :
    beta_(beta), threshold_(threshold) {}

bool SOMKEFixedThresholdMergingStrategy::ShouldMergeBePerformed(const vector<SOMSequenceEntry> &entries) const {

  for(unsigned int i = 1; i < entries.size(); ++i) {
    double current_modified_divergence = ComputeModifiedDivergenceOfSOMSequenceEntry(entries[i]);
    if(current_modified_divergence < threshold_) {
      return true;
    }
  }

  return false;
}

size_t SOMKEFixedThresholdMergingStrategy::FindEntryToMerge(const vector<SOMSequenceEntry> &entries) const {
  for(size_t i = 1; i < entries.size(); ++i) {
    double current_modified_divergence = ComputeModifiedDivergenceOfSOMSequenceEntry(entries[i]);
    if(current_modified_divergence < threshold_) {
      return i;
    }
  }
  return 0;
}

void SOMKEFixedThresholdMergingStrategy::SetDataWindowIterator(int *iterator) {
  data_window_iterator_ = iterator;
}

double SOMKEFixedThresholdMergingStrategy::ComputeModifiedDivergenceOfSOMSequenceEntry(
    const SOMSequenceEntry &entry) const {
  double modified_divergence = entry.kl_divergence;

  modified_divergence *= exp(-beta_ * (*data_window_iterator_ - (entry.range.first + entry.range.second) / 2));

  return modified_divergence;
}
