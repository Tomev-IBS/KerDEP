#ifndef KERDEP_ENHANCEDCLUSTERKERNELALGORITHM_H
#define KERDEP_ENHANCEDCLUSTERKERNELALGORITHM_H

#include "ClusterKernelsAlgorithm.h"

class EnhancedClusterKernelAlgorithm : public ClusterKernelsAlgorithm{
  /*
   * This class provides the same functionality as the original Cluster Kernel algorithm with the only
   * difference being that it has some utilities.
   */
  public:
    EnhancedClusterKernelAlgorithm(const int &m,
                                   ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
    std::vector<Point> GetErrorDomain(const int &dimension=0);
    Point GetKDEValuesOnDomain(std::vector<Point> domain);
    unsigned int FindMinimalValueOnDimension(const int &dimension);
    unsigned int FindMaximalValueOnDimension(const int &dimension);
};

#endif //KERDEP_ENHANCEDCLUSTERKERNELALGORITHM_H
