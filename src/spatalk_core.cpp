#include <Rcpp.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// ============================================================================
// 核心1: 超快速共表达计算（CCI最关键）
// ============================================================================

//' Fast co-expression calculation (C++ vectorized)
//' @param ligand_expr Ligand expression matrix (n_genes x n_cells)
//' @param receptor_expr Receptor expression matrix (n_genes x n_cells)
//' @return Co-expression ratios
//' @export
// [[Rcpp::export]]
NumericVector cpp_coexp_fast(NumericMatrix ligand_expr, NumericMatrix receptor_expr) {
    int n_genes = ligand_expr.nrow();
    int n_cells = ligand_expr.ncol();
    NumericVector result(n_genes);
    
    for(int g = 0; g < n_genes; g++) {
        int count = 0;
        for(int c = 0; c < n_cells; c++) {
            if(ligand_expr(g, c) > 0 && receptor_expr(g, c) > 0) {
                count++;
            }
        }
        result[g] = (double)count / n_cells;
    }
    
    return result;
}

// ============================================================================
// 核心2: 超快速置换检验（CCI瓶颈）
// ============================================================================

//' Fast permutation test - THE KEY OPTIMIZATION
//' @param st_data_mat Expression matrix
//' @param ligand_genes Ligand gene indices (1-based R indexing)
//' @param receptor_genes Receptor gene indices (1-based R indexing)
//' @param sender_cells Sender cell indices (1-based)
//' @param receiver_cells Receiver cell indices (1-based)
//' @param per_num Number of permutations
//' @param seed Random seed
//' @return List with real_ratio and pvalues
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
    
    // Calculate real co-expression ratios
    NumericVector real_ratios(n_lrs);
    for(int lr = 0; lr < n_lrs; lr++) {
        int lig_idx = ligand_genes[lr] - 1;  // Convert to 0-based
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
        real_ratios[lr] = (double)count / n_pairs;
    }
    
    // Pre-generate ALL random indices at once (HUGE speedup!)
    std::srand(seed);
    std::vector<std::vector<int>> perm_sender(per_num, std::vector<int>(n_pairs));
    std::vector<std::vector<int>> perm_receiver(per_num, std::vector<int>(n_pairs));
    
    for(int p = 0; p < per_num; p++) {
        for(int i = 0; i < n_pairs; i++) {
            perm_sender[p][i] = std::rand() % n_total_cells;
            perm_receiver[p][i] = std::rand() % n_total_cells;
        }
    }
    
    // Permutation test
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
            double perm_ratio = (double)count / n_pairs;
            
            if(perm_ratio >= real_ratios[lr]) {
                exceed_count++;
            }
        }
        
        pvalues[lr] = (double)exceed_count / per_num;
    }
    
    return List::create(
        Named("real_ratios") = real_ratios,
        Named("pvalues") = pvalues
    );
}

// ============================================================================
// 核心3: 超快速随机游走（CCI关键）
// ============================================================================

//' Fast random walk for TF scoring
//' @param ggi_src Source gene names
//' @param ggi_dest Destination gene names  
//' @param receptor_name Starting receptor name
//' @param tf_names TF names to score
//' @param n_walks Number of walks
//' @param max_hop Maximum hops
//' @param seed Random seed
//' @return TF scores
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
    
    // Build adjacency list (O(1) lookup!)
    std::map<std::string, std::vector<int>> adj_list;
    for(int i = 0; i < n_edges; i++) {
        std::string src = as<std::string>(ggi_src[i]);
        adj_list[src].push_back(i);
    }
    
    // Build TF set and score map
    std::set<std::string> tf_set;
    std::map<std::string, int> tf_scores;
    for(int i = 0; i < n_tfs; i++) {
        std::string tf = as<std::string>(tf_names[i]);
        tf_set.insert(tf);
        tf_scores[tf] = 0;
    }
    
    // Random walks
    std::srand(seed);
    
    for(int walk = 0; walk < n_walks; walk++) {
        std::string current = receptor_name;
        
        for(int hop = 0; hop < max_hop; hop++) {
            // Fast adjacency lookup
            auto it = adj_list.find(current);
            if(it == adj_list.end()) break;
            
            const std::vector<int>& neighbors = it->second;
            if(neighbors.empty()) break;
            
            // Random next node
            int rand_edge = neighbors[std::rand() % neighbors.size()];
            std::string next_node = as<std::string>(ggi_dest[rand_edge]);
            
            // Update TF score
            if(tf_set.count(next_node)) {
                tf_scores[next_node]++;
            }
            
            current = next_node;
        }
    }
    
    // Convert to result
    NumericVector result(n_tfs);
    for(int i = 0; i < n_tfs; i++) {
        std::string tf = as<std::string>(tf_names[i]);
        result[i] = (double)tf_scores[tf] / n_walks;
    }
    
    return result;
}

// ============================================================================
// 核心4: 超快速细胞采样（解卷积优化）
// ============================================================================

//' Fast cell sampling with early stopping
//' @param spot_ndata Spot expression vector
//' @param sc_ndata_mat SC expression matrix
//' @param cell_indices_by_type List of cell indices for each celltype
//' @param spot_celltypes Cell types needed for this spot
//' @param iter_num Maximum iterations
//' @param tolerance Convergence tolerance
//' @param seed Random seed
//' @return Best cell combination
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
    
    // Pre-process cell indices
    std::map<std::string, std::vector<int>> celltype_map;
    CharacterVector celltype_names = cell_indices_by_type.names();
    
    for(int i = 0; i < celltype_names.size(); i++) {
        std::string ct_name = as<std::string>(celltype_names[i]);
        IntegerVector cells = cell_indices_by_type[i];
        for(int j = 0; j < cells.size(); j++) {
            celltype_map[ct_name].push_back(cells[j] - 1);  // 0-based
        }
    }
    
    std::srand(seed);
    
    double best_cor = -1.0;
    IntegerVector best_cells(n_cells_in_spot);
    int no_improve = 0;
    
    for(int iter = 0; iter < iter_num; iter++) {
        // Sample cells
        IntegerVector sampled_cells(n_cells_in_spot);
        for(int j = 0; j < n_cells_in_spot; j++) {
            std::string ct = as<std::string>(spot_celltypes[j]);
            const std::vector<int>& available = celltype_map[ct];
            sampled_cells[j] = available[std::rand() % available.size()];
        }
        
        // Predict expression
        NumericVector pred(n_genes);
        for(int g = 0; g < n_genes; g++) {
            double sum = 0.0;
            for(int j = 0; j < n_cells_in_spot; j++) {
                sum += sc_ndata_mat(g, sampled_cells[j]);
            }
            pred[g] = sum;
        }
        
        // Calculate correlation
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
        
        double cor = num / (std::sqrt(denom_spot) * std::sqrt(denom_pred));
        
        // Track best
        if(cor > best_cor) {
            best_cor = cor;
            best_cells = sampled_cells;
            no_improve = 0;
        } else {
            no_improve++;
        }
        
        // AGGRESSIVE early stopping
        if(iter >= 50 && (no_improve >= 30 || best_cor > 0.99)) {
            break;
        }
    }
    
    // Convert to 1-based indexing for R
    for(int i = 0; i < n_cells_in_spot; i++) {
        best_cells[i]++;
    }
    
    return List::create(
        Named("cell_indices") = best_cells,
        Named("correlation") = best_cor
    );
}

// ============================================================================
// 核心5: 批量共表达计算（向量化）
// ============================================================================

//' Batch co-expression for gene-gene interactions
//' @param st_data_mat Expression matrix
//' @param src_genes Source gene indices (1-based)
//' @param dest_genes Dest gene indices (1-based)
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
        result[i] = (double)count / n_cells;
    }
    
    return result;
}

// ============================================================================
// 核心6: 距离矩阵快速计算
// ============================================================================

//' Fast Euclidean distance matrix
//' @param x X coordinates
//' @param y Y coordinates
//' @return Distance matrix
//' @export
// [[Rcpp::export]]
NumericMatrix cpp_fast_dist(NumericVector x, NumericVector y) {
    int n = x.size();
    NumericMatrix dist(n, n);
    
    for(int i = 0; i < n; i++) {
        dist(i, i) = 0.0;
        for(int j = i + 1; j < n; j++) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double d = std::sqrt(dx * dx + dy * dy);
            dist(i, j) = d;
            dist(j, i) = d;
        }
    }
    
    return dist;
}

// ============================================================================
// 核心7: 快速找到K近邻
// ============================================================================

//' Fast K-nearest neighbors
//' @param dist_mat Distance matrix
//' @param query_idx Query cell index (0-based)
//' @param k Number of neighbors
//' @return Neighbor indices (0-based)
//' @export
// [[Rcpp::export]]
IntegerVector cpp_knn(NumericMatrix dist_mat, int query_idx, int k) {
    int n = dist_mat.ncol();
    
    // Get distances and indices
    std::vector<std::pair<double, int>> dist_idx;
    for(int i = 0; i < n; i++) {
        if(i != query_idx && dist_mat(query_idx, i) > 0) {
            dist_idx.push_back(std::make_pair(dist_mat(query_idx, i), i));
        }
    }
    
    // Partial sort to get k smallest
    int actual_k = std::min(k, (int)dist_idx.size());
    std::partial_sort(dist_idx.begin(), 
                     dist_idx.begin() + actual_k,
                     dist_idx.end());
    
    IntegerVector result(actual_k);
    for(int i = 0; i < actual_k; i++) {
        result[i] = dist_idx[i].second;
    }
    
    return result;
}

// ============================================================================
// 核心8: 批量相关系数计算
// ============================================================================

//' Fast batch correlation
//' @param vec1 First vector
//' @param mat Matrix where each column is a vector to correlate with vec1
//' @return Correlations
//' @export
// [[Rcpp::export]]
NumericVector cpp_batch_cor(NumericVector vec1, NumericMatrix mat) {
    int n = vec1.size();
    int n_cols = mat.ncol();
    NumericVector result(n_cols);
    
    // Calculate mean of vec1
    double mean1 = 0;
    for(int i = 0; i < n; i++) {
        mean1 += vec1[i];
    }
    mean1 /= n;
    
    for(int col = 0; col < n_cols; col++) {
        // Calculate mean of this column
        double mean2 = 0;
        for(int i = 0; i < n; i++) {
            mean2 += mat(i, col);
        }
        mean2 /= n;
        
        // Calculate correlation
        double num = 0, denom1 = 0, denom2 = 0;
        for(int i = 0; i < n; i++) {
            double diff1 = vec1[i] - mean1;
            double diff2 = mat(i, col) - mean2;
            num += diff1 * diff2;
            denom1 += diff1 * diff1;
            denom2 += diff2 * diff2;
        }
        
        if(denom1 > 0 && denom2 > 0) {
            result[col] = num / (std::sqrt(denom1) * std::sqrt(denom2));
        } else {
            result[col] = 0.0;
        }
    }
    
    return result;
}

