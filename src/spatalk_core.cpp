// SpaTalk Core C++ Functions
// Copyright (c) 2022-2024 Zaoqu Liu

#include <Rcpp.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <random>
#include <string>

using namespace Rcpp;

// ============================================================================
// Co-expression calculation
// ============================================================================

//' Calculate co-expression ratios
//' @param ligand_expr Ligand expression matrix (n_genes x n_cells)
//' @param receptor_expr Receptor expression matrix (n_genes x n_cells)
//' @return Co-expression ratios
//' @export
// [[Rcpp::export]]
NumericVector cpp_coexp_fast(NumericMatrix ligand_expr, NumericMatrix receptor_expr) {
    int n_genes = ligand_expr.nrow();
    int n_cells = ligand_expr.ncol();
    
    // Input validation
    if(n_genes == 0 || n_cells == 0) {
        return NumericVector(0);
    }
    if(receptor_expr.nrow() != n_genes || receptor_expr.ncol() != n_cells) {
        stop("Ligand and receptor expression matrices must have the same dimensions");
    }
    
    NumericVector result(n_genes);
    
    for(int g = 0; g < n_genes; g++) {
        int count = 0;
        for(int c = 0; c < n_cells; c++) {
            if(ligand_expr(g, c) > 0 && receptor_expr(g, c) > 0) {
                count++;
            }
        }
        result[g] = static_cast<double>(count) / n_cells;
    }
    
    return result;
}

// ============================================================================
// Permutation test for ligand-receptor pairs
// ============================================================================

//' Permutation test for LR co-expression significance
//' @param st_data_mat Expression matrix
//' @param ligand_genes Ligand gene indices (1-based R indexing)
//' @param receptor_genes Receptor gene indices (1-based R indexing)
//' @param sender_cells Sender cell indices (1-based)
//' @param receiver_cells Receiver cell indices (1-based)
//' @param per_num Number of permutations
//' @param seed Random seed for reproducibility
//' @return List with real_ratios and pvalues
//' @export
// [[Rcpp::export]]
List cpp_permutation_test(NumericMatrix st_data_mat,
                          IntegerVector ligand_genes,
                          IntegerVector receptor_genes,
                          IntegerVector sender_cells,
                          IntegerVector receiver_cells,
                          int per_num = 1000,
                          int seed = 123) {
    
    int n_lrs = ligand_genes.size();
    int n_pairs = sender_cells.size();
    int n_total_cells = st_data_mat.ncol();
    int n_genes = st_data_mat.nrow();
    
    // Input validation
    if(n_lrs == 0 || n_pairs == 0 || n_total_cells == 0) {
        return List::create(
            Named("real_ratios") = NumericVector(0),
            Named("pvalues") = NumericVector(0)
        );
    }
    
    // Validate indices are within bounds
    for(int i = 0; i < n_lrs; i++) {
        if(ligand_genes[i] < 1 || ligand_genes[i] > n_genes ||
           receptor_genes[i] < 1 || receptor_genes[i] > n_genes) {
            stop("Gene indices out of bounds");
        }
    }
    for(int i = 0; i < n_pairs; i++) {
        if(sender_cells[i] < 1 || sender_cells[i] > n_total_cells ||
           receiver_cells[i] < 1 || receiver_cells[i] > n_total_cells) {
            stop("Cell indices out of bounds");
        }
    }
    
    // Calculate observed co-expression ratios
    NumericVector real_ratios(n_lrs);
    for(int lr = 0; lr < n_lrs; lr++) {
        int lig_idx = ligand_genes[lr] - 1;
        int rec_idx = receptor_genes[lr] - 1;
        
        int count = 0;
        for(int p = 0; p < n_pairs; p++) {
            int sender_idx = sender_cells[p] - 1;
            int receiver_idx = receiver_cells[p] - 1;
            
            if(st_data_mat(lig_idx, sender_idx) > 0 && 
               st_data_mat(rec_idx, receiver_idx) > 0) {
                count++;
            }
        }
        real_ratios[lr] = static_cast<double>(count) / n_pairs;
    }
    
    // Initialize random number generator (Mersenne Twister for reproducibility)
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> dist(0, n_total_cells - 1);
    
    // Generate permuted indices
    std::vector<std::vector<int>> perm_sender(per_num, std::vector<int>(n_pairs));
    std::vector<std::vector<int>> perm_receiver(per_num, std::vector<int>(n_pairs));
    
    for(int p = 0; p < per_num; p++) {
        for(int i = 0; i < n_pairs; i++) {
            perm_sender[p][i] = dist(gen);
            perm_receiver[p][i] = dist(gen);
        }
    }
    
    // Calculate p-values
    NumericVector pvalues(n_lrs);
    
    for(int lr = 0; lr < n_lrs; lr++) {
        int lig_idx = ligand_genes[lr] - 1;
        int rec_idx = receptor_genes[lr] - 1;
        
        int exceed_count = 0;
        
        for(int p = 0; p < per_num; p++) {
            int count = 0;
            for(int i = 0; i < n_pairs; i++) {
                if(st_data_mat(lig_idx, perm_sender[p][i]) > 0 && 
                   st_data_mat(rec_idx, perm_receiver[p][i]) > 0) {
                    count++;
                }
            }
            double perm_ratio = static_cast<double>(count) / n_pairs;
            
            if(perm_ratio >= real_ratios[lr]) {
                exceed_count++;
            }
        }
        
        pvalues[lr] = static_cast<double>(exceed_count) / per_num;
    }
    
    return List::create(
        Named("real_ratios") = real_ratios,
        Named("pvalues") = pvalues
    );
}

// ============================================================================
// Random walk algorithm for transcription factor scoring
// ============================================================================

//' Random walk on gene-gene interaction network
//' @param ggi_src Source gene names
//' @param ggi_dest Destination gene names  
//' @param receptor_name Starting receptor name
//' @param tf_names TF names to score
//' @param n_walks Number of random walks
//' @param max_hop Maximum number of hops per walk
//' @param seed Random seed for reproducibility
//' @return TF visit frequency scores
//' @export
// [[Rcpp::export]]
NumericVector cpp_random_walk(CharacterVector ggi_src,
                               CharacterVector ggi_dest,
                               std::string receptor_name,
                               CharacterVector tf_names,
                               int n_walks = 10000,
                               int max_hop = 10,
                               int seed = 123) {
    
    int n_edges = ggi_src.size();
    int n_tfs = tf_names.size();
    
    // Build adjacency list
    std::map<std::string, std::vector<int>> adj_list;
    for(int i = 0; i < n_edges; i++) {
        std::string src = as<std::string>(ggi_src[i]);
        adj_list[src].push_back(i);
    }
    
    // Initialize TF tracking
    std::set<std::string> tf_set;
    std::map<std::string, int> tf_scores;
    for(int i = 0; i < n_tfs; i++) {
        std::string tf = as<std::string>(tf_names[i]);
        tf_set.insert(tf);
        tf_scores[tf] = 0;
    }
    
    // Perform random walks
    std::mt19937 gen(seed);
    
    for(int walk = 0; walk < n_walks; walk++) {
        std::string current = receptor_name;
        
        for(int hop = 0; hop < max_hop; hop++) {
            auto it = adj_list.find(current);
            if(it == adj_list.end()) break;
            
            const std::vector<int>& neighbors = it->second;
            if(neighbors.empty()) break;
            
            std::uniform_int_distribution<size_t> neighbor_dist(0, neighbors.size() - 1);
            int rand_edge = neighbors[neighbor_dist(gen)];
            std::string next_node = as<std::string>(ggi_dest[rand_edge]);
            
            if(tf_set.count(next_node)) {
                tf_scores[next_node]++;
            }
            
            current = next_node;
        }
    }
    
    // Convert to result vector
    NumericVector result(n_tfs);
    for(int i = 0; i < n_tfs; i++) {
        std::string tf = as<std::string>(tf_names[i]);
        result[i] = static_cast<double>(tf_scores[tf]) / n_walks;
    }
    
    return result;
}

// ============================================================================
// Cell sampling for deconvolution
// ============================================================================

//' Sample cells to reconstruct spot expression
//' @param spot_ndata Spot expression vector
//' @param sc_ndata_mat Single-cell expression matrix
//' @param cell_indices_by_type List of cell indices for each cell type
//' @param spot_celltypes Cell types to sample for this spot
//' @param iter_num Maximum number of iterations
//' @param tolerance Convergence tolerance (unused, kept for API compatibility)
//' @param seed Random seed for reproducibility
//' @return Best cell combination and correlation
//' @export
// [[Rcpp::export]]
List cpp_fast_sampling(NumericVector spot_ndata,
                       NumericMatrix sc_ndata_mat,
                       List cell_indices_by_type,
                       CharacterVector spot_celltypes,
                       int iter_num = 200,
                       double tolerance = 0.001,
                       int seed = 123) {
    
    int n_genes = spot_ndata.size();
    int n_cells_in_spot = spot_celltypes.size();
    
    // Build cell type to indices mapping
    std::map<std::string, std::vector<int>> celltype_map;
    CharacterVector celltype_names = cell_indices_by_type.names();
    
    for(int i = 0; i < celltype_names.size(); i++) {
        std::string ct_name = as<std::string>(celltype_names[i]);
        IntegerVector cells = cell_indices_by_type[i];
        for(int j = 0; j < cells.size(); j++) {
            celltype_map[ct_name].push_back(cells[j] - 1);
        }
    }
    
    std::mt19937 gen(seed);
    
    double best_cor = -1.0;
    IntegerVector best_cells(n_cells_in_spot);
    int no_improve = 0;
    
    for(int iter = 0; iter < iter_num; iter++) {
        // Sample cells for each cell type
        IntegerVector sampled_cells(n_cells_in_spot);
        bool valid_sample = true;
        for(int j = 0; j < n_cells_in_spot; j++) {
            std::string ct = as<std::string>(spot_celltypes[j]);
            auto it = celltype_map.find(ct);
            if(it == celltype_map.end() || it->second.empty()) {
                valid_sample = false;
                break;
            }
            const std::vector<int>& available = it->second;
            std::uniform_int_distribution<size_t> cell_dist(0, available.size() - 1);
            sampled_cells[j] = available[cell_dist(gen)];
        }
        if(!valid_sample) continue;
        
        // Predict spot expression
        NumericVector pred(n_genes);
        for(int g = 0; g < n_genes; g++) {
            double sum = 0.0;
            for(int j = 0; j < n_cells_in_spot; j++) {
                sum += sc_ndata_mat(g, sampled_cells[j]);
            }
            pred[g] = sum;
        }
        
        // Calculate Pearson correlation
        double mean_spot = 0, mean_pred = 0;
        for(int g = 0; g < n_genes; g++) {
            mean_spot += spot_ndata[g];
            mean_pred += pred[g];
        }
        mean_spot /= n_genes;
        mean_pred /= n_genes;
        
        double num = 0, denom_spot = 0, denom_pred = 0;
        for(int g = 0; g < n_genes; g++) {
            double diff_spot = spot_ndata[g] - mean_spot;
            double diff_pred = pred[g] - mean_pred;
            num += diff_spot * diff_pred;
            denom_spot += diff_spot * diff_spot;
            denom_pred += diff_pred * diff_pred;
        }
        
        double cor_val = 0.0;
        double denom_product = std::sqrt(denom_spot) * std::sqrt(denom_pred);
        if(denom_product > 0) {
            cor_val = num / denom_product;
        }
        
        // Track best result
        if(cor_val > best_cor) {
            best_cor = cor_val;
            best_cells = clone(sampled_cells);
            no_improve = 0;
        } else {
            no_improve++;
        }
        
        // Early stopping when converged
        if(iter >= 50 && (no_improve >= 30 || best_cor > 0.99)) {
            break;
        }
    }
    
    // Convert to 1-based R indexing
    for(int i = 0; i < n_cells_in_spot; i++) {
        best_cells[i]++;
    }
    
    return List::create(
        Named("cell_indices") = best_cells,
        Named("correlation") = best_cor
    );
}

// ============================================================================
// Batch co-expression calculation
// ============================================================================

//' Calculate co-expression for multiple gene pairs
//' @param st_data_mat Expression matrix
//' @param src_genes Source gene indices (1-based)
//' @param dest_genes Destination gene indices (1-based)
//' @param cell_indices Cell indices to use (1-based)
//' @return Co-expression ratios
//' @export
// [[Rcpp::export]]
NumericVector cpp_batch_coexp(NumericMatrix st_data_mat,
                               IntegerVector src_genes,
                               IntegerVector dest_genes,
                               IntegerVector cell_indices) {
    
    int n_pairs = src_genes.size();
    int n_cells = cell_indices.size();
    int n_total_genes = st_data_mat.nrow();
    int n_total_cells = st_data_mat.ncol();
    
    // Input validation
    if(n_pairs == 0 || n_cells == 0) {
        return NumericVector(0);
    }
    if(dest_genes.size() != n_pairs) {
        stop("src_genes and dest_genes must have the same length");
    }
    
    // Validate indices
    for(int i = 0; i < n_pairs; i++) {
        if(src_genes[i] < 1 || src_genes[i] > n_total_genes ||
           dest_genes[i] < 1 || dest_genes[i] > n_total_genes) {
            stop("Gene indices out of bounds");
        }
    }
    for(int i = 0; i < n_cells; i++) {
        if(cell_indices[i] < 1 || cell_indices[i] > n_total_cells) {
            stop("Cell indices out of bounds");
        }
    }
    
    NumericVector result(n_pairs);
    
    for(int i = 0; i < n_pairs; i++) {
        int src_idx = src_genes[i] - 1;
        int dest_idx = dest_genes[i] - 1;
        
        int count = 0;
        for(int c = 0; c < n_cells; c++) {
            int cell_idx = cell_indices[c] - 1;
            if(st_data_mat(src_idx, cell_idx) > 0 && 
               st_data_mat(dest_idx, cell_idx) > 0) {
                count++;
            }
        }
        result[i] = static_cast<double>(count) / n_cells;
    }
    
    return result;
}

// ============================================================================
// Distance matrix computation
// ============================================================================

//' Calculate Euclidean distance matrix
//' @param x X coordinates
//' @param y Y coordinates
//' @return Distance matrix
//' @export
// [[Rcpp::export]]
NumericMatrix cpp_fast_dist(NumericVector x, NumericVector y) {
    int n = x.size();
    NumericMatrix dist_mat(n, n);
    
    for(int i = 0; i < n; i++) {
        dist_mat(i, i) = 0.0;
        for(int j = i + 1; j < n; j++) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double d = std::sqrt(dx * dx + dy * dy);
            dist_mat(i, j) = d;
            dist_mat(j, i) = d;
        }
    }
    
    return dist_mat;
}

// ============================================================================
// K-nearest neighbors
// ============================================================================

//' Find K nearest neighbors
//' @param dist_mat Distance matrix
//' @param query_idx Query index (0-based)
//' @param k Number of neighbors
//' @return Neighbor indices (0-based)
//' @export
// [[Rcpp::export]]
IntegerVector cpp_knn(NumericMatrix dist_mat, int query_idx, int k) {
    int n = dist_mat.ncol();
    
    std::vector<std::pair<double, int>> dist_idx;
    dist_idx.reserve(n);
    
    for(int i = 0; i < n; i++) {
        if(i != query_idx && dist_mat(query_idx, i) > 0) {
            dist_idx.push_back(std::make_pair(dist_mat(query_idx, i), i));
        }
    }
    
    int actual_k = std::min(k, static_cast<int>(dist_idx.size()));
    if(actual_k > 0) {
        std::partial_sort(dist_idx.begin(), 
                         dist_idx.begin() + actual_k,
                         dist_idx.end());
    }
    
    IntegerVector result(actual_k);
    for(int i = 0; i < actual_k; i++) {
        result[i] = dist_idx[i].second;
    }
    
    return result;
}

// ============================================================================
// Batch correlation computation
// ============================================================================

//' Calculate correlation between a vector and matrix columns
//' @param vec1 Reference vector
//' @param mat Matrix where each column is compared to vec1
//' @return Pearson correlation coefficients
//' @export
// [[Rcpp::export]]
NumericVector cpp_batch_cor(NumericVector vec1, NumericMatrix mat) {
    int n = vec1.size();
    int n_cols = mat.ncol();
    NumericVector result(n_cols);
    
    // Pre-compute vec1 statistics
    double mean1 = 0;
    for(int i = 0; i < n; i++) {
        mean1 += vec1[i];
    }
    mean1 /= n;
    
    std::vector<double> diff1(n);
    double denom1 = 0;
    for(int i = 0; i < n; i++) {
        diff1[i] = vec1[i] - mean1;
        denom1 += diff1[i] * diff1[i];
    }
    double sqrt_denom1 = std::sqrt(denom1);
    
    for(int col = 0; col < n_cols; col++) {
        double mean2 = 0;
        for(int i = 0; i < n; i++) {
            mean2 += mat(i, col);
        }
        mean2 /= n;
        
        double num = 0, denom2 = 0;
        for(int i = 0; i < n; i++) {
            double diff2 = mat(i, col) - mean2;
            num += diff1[i] * diff2;
            denom2 += diff2 * diff2;
        }
        
        double denom_product = sqrt_denom1 * std::sqrt(denom2);
        if(denom_product > 0) {
            result[col] = num / denom_product;
        } else {
            result[col] = 0.0;
        }
    }
    
    return result;
}
