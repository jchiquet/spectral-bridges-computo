#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_affinity(NumericMatrix kmeans_centers, 
                                    IntegerVector kmeans_labels, 
                                    NumericMatrix X_centered, 
                                    IntegerVector kmeans_size, 
                                    int n_nodes) {
  NumericMatrix norm_diff(n_nodes, n_nodes);
  NumericMatrix affinity(n_nodes, n_nodes);
  
  for (int k = 0; k < n_nodes; ++k) {
    for (int l = k; l < n_nodes; ++l) {
      if (k != l) {
        // Compute diff_kl
        NumericVector diff_kl = kmeans_centers(k, _) - kmeans_centers(l, _);
        
        // Compute norm_diff
        double norm_sum = std::inner_product(diff_kl.begin(), diff_kl.end(), diff_kl.begin(), 0.0);
        norm_diff(k, l) = norm_sum;
        norm_diff(l, k) = norm_sum; // Since norm_diff is symmetric
        
        // Calculate term1
        double term1 = 0;
        for (int i = 0; i < X_centered.nrow(); ++i) {
          if (kmeans_labels[i] == k+1) {
            double dot_product = std::inner_product(diff_kl.begin(), diff_kl.end(), X_centered(i, _).begin(), 0.0);
            term1 += std::max(0.0, dot_product);
          }
        }
        
        // Calculate term2
        double term2 = 0;
        for (int i = 0; i < X_centered.nrow(); ++i) {
          if (kmeans_labels[i] == l+1) {
            double dot_product = std::inner_product(diff_kl.begin(), diff_kl.end(), X_centered(i, _).begin(), 0.0);
            term2 += std::max(0.0, -dot_product);
          }
        }
        
        // Calculate affinity
        affinity(k, l) = (term1 + term2) / (norm_diff(k, l) * (kmeans_size[k] + kmeans_size[l]));
        affinity(l, k) = affinity(k, l); // Since affinity is symmetric
      }
    }
  }
  
  return affinity;
}
