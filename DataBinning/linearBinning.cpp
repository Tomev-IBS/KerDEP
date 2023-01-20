//
// Created by tomev on 1/9/2023.
//

#include "linearBinning.h"
#include <numeric>

using namespace arma;

LinearBinning::LinearBinning(const arma::mat &data, const int &gridSize) {

  _data = data;
  _gridSize = gridSize; // Assume same size (number of points) for every dimension.

  compute_grid_info();
}

void LinearBinning::compute_grid_info() {
  /*
   * A grid info should contain the data about all the grid values in every dimension. This way we can create the points
   * dynamically.
   */

  arma::vec stDevs = compute_standard_deviations();

  _gridInfo = mat(_gridSize, stDevs.size());

  for(int col_id = 0; col_id < stDevs.size(); ++col_id){

    double col_min = min(_data.col(col_id)) - stDevs[col_id];
    double col_max = max(_data.col(col_id)) + stDevs[col_id];
    double grid_width = (col_max - col_min) / (_gridSize - 1);

    for(int row_id = 0; row_id < _gridSize; ++row_id){
      _gridInfo(row_id, col_id) = col_min + row_id * grid_width;
    }
  }

}

vec LinearBinning::compute_standard_deviations() const{

  arma::vec stDevs = arma::vec(arma::size(_data)[1]);

  for(int i = 0; i < arma::size(_data)[1]; ++i){
    stDevs[i] = stddev(_data.col(i));
  }

  return stDevs;
}



std::map<std::vector<double>, double> LinearBinning::compute_grid_counts() const {
  // The idea is taken from Wands (1994) Fast Computation of Multivariate Kernel Estimators.
  std::map<std::vector<double>, double> counts;


  for(int row_id = 0; row_id < size(_data)[0]; ++row_id){

    vec datum = _data.row(row_id); // It represents multidimensional data point.

    std::vector<double> weights = {};
    std::vector<int> grid_indices = {};

    for(int dimension = 0; dimension < datum.size(); ++dimension){
      double l_val = l(datum, dimension);
      int j = int(l_val);

      grid_indices.push_back(j);
      weights.push_back(1 - (l_val - j));
    }

    auto relevant_points_indices = get_points_indices(grid_indices);
    auto relevant_points = get_relevant_points(relevant_points_indices);
    auto respective_weights = get_respective_weights(weights);

    for(int i = 0; i < relevant_points.size(); ++i) {

      if(counts.count(relevant_points[i]) == 0) {
        counts[relevant_points[i]] = 0;
      }

      counts[relevant_points[i]] += weights[i];

    }
  }

  return counts;
}

double LinearBinning::l(vec x, int dimension) const{
  // A helper function from Fan & Marron 1994.
  vec grid = _gridInfo.col(dimension);
  double grid_width = grid[1] - grid[0];
  return (1 / grid_width) * (x[dimension] - min(grid));
}

std::vector<std::vector<int>> LinearBinning::get_points_indices(std::vector<int> found_indices) const {

  std::vector<std::vector<int>> points_indices = {};

  int j = found_indices.back();
  found_indices.pop_back();

  if(found_indices.empty()){
    points_indices.push_back({j});
    points_indices.push_back({j + 1});
    return points_indices;
  }

  std::vector<std::vector<int>> subsets = get_points_indices(found_indices);

  for(int i = 0; i < 2; ++i){ // We're only adding j and j + 1 for every vector in subsets. Hence, 2 in the loop.
    points_indices.insert(points_indices.begin(), subsets.begin(), subsets.end());
    points_indices[points_indices.size() - 2].push_back(j + i);
    points_indices[points_indices.size() - 1].push_back(j + i);
  }

  return points_indices;
}

std::vector<std::vector<double>> LinearBinning::get_relevant_points(std::vector<std::vector<int>> indices_sets) const{
  std::vector<std::vector<double>> relevant_points = {};

  for(auto set : indices_sets){
    relevant_points.push_back({});

    for(int i = 0; i < set.size(); ++i){
      relevant_points.back().push_back(_gridInfo.col(i)[set[i]]);
    }

  }

  return relevant_points;
}

std::vector<double> LinearBinning::get_respective_weights(const std::vector<double> &found_weights) const {

  auto partial_weights_sets = get_partial_weights_sets(found_weights);

  std::vector<double> weights = {};

  for(auto set : partial_weights_sets){
    weights.push_back(std::accumulate(set.begin(), set.end(), 1, std::multiplies<double>{}));
  }

  return weights;
}

std::vector<std::vector<double>> LinearBinning::get_partial_weights_sets(std::vector<double> weights) const {
  std::vector<std::vector<double>> partial_weights_sets = {};

  double w = weights.back();
  weights.pop_back();

  if(weights.empty()){
    partial_weights_sets.push_back({w});
    partial_weights_sets.push_back({1 - w});
    return partial_weights_sets;
  }

  std::vector<std::vector<double>> subsets = get_partial_weights_sets(weights);

  for(int i = 0; i < 2; ++i){ // We're only adding w and w + 1 for every vector in subsets. Hence, 2 in the loop.
    partial_weights_sets.insert(partial_weights_sets.begin(), subsets.begin(), subsets.end());

    // Update w for second batch.
    if(i == 1){
      w = 1 - w;
    }

    partial_weights_sets[partial_weights_sets.size() - 2].push_back(w);
    partial_weights_sets[partial_weights_sets.size() - 1].push_back(w);
  }

  return partial_weights_sets;
}



