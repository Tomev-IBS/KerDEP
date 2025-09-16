//
// Created by tomev on 1/9/2023.
// A class that implements multivariate linear binning algorithm, as presented by Wand.
//

#ifndef DEDSTA_LINEARBINNING_H
#define DEDSTA_LINEARBINNING_H

#include <armadillo>
#include <map>
#include <vector>

using namespace arma;

class LinearBinning {

  public:
    LinearBinning(const arma::mat &data, const int &gridSize);
    [[nodiscard]] std::map<std::vector<double>, double> compute_grid_counts() const;

  protected:
    mat _data;
    mat _gridInfo;
    int _gridSize;

    void compute_grid_info();
    [[nodiscard]] vec compute_standard_deviations() const;

    // Helper functions
    [[nodiscard]] double l(vec x, int dimension) const;
    [[nodiscard]] std::vector<std::vector<int>> get_points_indices(std::vector<int> found_indices) const; // Recursive
    [[nodiscard]] std::vector<std::vector<double>> get_relevant_points(std::vector<std::vector<int>> indices_sets) const;
    [[nodiscard]] std::vector<double> get_respective_weights(const std::vector<double> &found_weights) const;
    [[nodiscard]] std::vector<std::vector<double>> get_partial_weights_sets(std::vector<double> weights) const; // Recursive

};

#endif //DEDSTA_LINEARBINNING_H
