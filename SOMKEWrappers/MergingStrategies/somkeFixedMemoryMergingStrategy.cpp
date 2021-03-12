//
// Created by tomev on 12/03/2021.
//

#include "somkeFixedMemoryMergingStrategy.h"

SOMKEFixedMemoryMergingStrategy::SOMKEFixedMemoryMergingStrategy(const size_t &maximal_number_of_entries_in_sequence,
                                                                 const double &beta) :
  maximal_number_of_entries_in_sequence_(maximal_number_of_entries_in_sequence), beta_(beta)
{ }

bool SOMKEFixedMemoryMergingStrategy::ShouldMergeBePerformed(const vector<SOMSequenceEntry> &entries) const {
  return entries.size() > maximal_number_of_entries_in_sequence_;
}

size_t SOMKEFixedMemoryMergingStrategy::FindEntryToMerge(const vector<SOMSequenceEntry> &entries) const {
  size_t smallest_divergence_som_sequence_entry_index = 1;
  double smallest_divergence = ComputeModifiedDivergenceOfSOMSequenceEntry(entries[1]);

  for(unsigned int i = 2; i < entries.size(); ++i){
    double current_modified_divergence = ComputeModifiedDivergenceOfSOMSequenceEntry(entries[i]);
    if(smallest_divergence > current_modified_divergence){
      smallest_divergence_som_sequence_entry_index = i;
      smallest_divergence = current_modified_divergence;
    }
  }

  return smallest_divergence_som_sequence_entry_index;
}

double SOMKEFixedMemoryMergingStrategy::ComputeModifiedDivergenceOfSOMSequenceEntry(const SOMSequenceEntry &entry) const {
  double modified_divergence = entry.kl_divergence;

  modified_divergence *= exp(- beta_ * (*data_window_iterator_ - (entry.range.first + entry.range.second) / 2));

  return modified_divergence;
}

void SOMKEFixedMemoryMergingStrategy::SetDataWindowIterator(int *iterator) {
  data_window_iterator_ = iterator;
}
