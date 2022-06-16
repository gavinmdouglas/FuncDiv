// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

//  Function assumes that the func_tab column names are equal (and identically ordered) to the abun_tab rownames
// [[Rcpp::export]]
List prep_all_sample_func_vec(NumericMatrix abun_tab, NumericMatrix func_tab) {
  
  int num_samples = abun_tab.ncol();
  int num_func = func_tab.nrow();
  
  int num_instances = num_func * num_samples;
  
  List all_func_sample_taxa_abun(num_instances);
  
  int instance_num = 0;
  
  for (int i = 0; i < num_func; ++i) {
    
    NumericVector func_row = func_tab(i, _); 
    LogicalVector func_row_nonzero = func_row > 0;
    
    for (int j = 0; j < num_samples; ++j) {
      
      NumericVector sample_taxa_abun = abun_tab(_, j);
      sample_taxa_abun = sample_taxa_abun[func_row_nonzero];
      LogicalVector sample_taxa_abun_nonzero_flag = sample_taxa_abun > 0;
      all_func_sample_taxa_abun[instance_num] = sample_taxa_abun[sample_taxa_abun_nonzero_flag];
      
      instance_num += 1;
      
    }
    
  }
  
  return all_func_sample_taxa_abun;
  
}


// [[Rcpp::export]]
List prep_all_sample_func_taxa_vec(NumericMatrix abun_tab, NumericMatrix func_tab) {
  
  int num_samples = abun_tab.ncol();
  int num_func = func_tab.nrow();
  
  int num_instances = num_func * num_samples;
  
  CharacterVector abun_taxa = rownames(abun_tab);
  
  List all_func_sample_taxa_present(num_instances);
  
  int instance_num = 0;
  
  for (int i = 0; i < num_func; ++i) {
    
    NumericVector func_row = func_tab(i, _); 
    LogicalVector func_row_nonzero = func_row > 0;
    
    for (int j = 0; j < num_samples; ++j) {
      
      NumericVector sample_taxa_abun = abun_tab(_, j);
        
      sample_taxa_abun = sample_taxa_abun[func_row_nonzero];

      CharacterVector func_sample_abun_taxa = abun_taxa[func_row_nonzero];
      
      LogicalVector sample_taxa_abun_nonzero_flag = sample_taxa_abun > 0;
  
      all_func_sample_taxa_present[instance_num] = func_sample_abun_taxa[sample_taxa_abun_nonzero_flag];

      instance_num += 1;
      
    }
    
  }
  
  return all_func_sample_taxa_present;
  
}


// [[Rcpp::export]]
List prep_func_contributor_dimnames(arma::mat abun_tab, 
                                    arma::mat func_tab) {
  
  int num_func = func_tab.n_rows;
  
  List all_func_contrib_abun(num_func);

  for (int i = 0; i < num_func; ++i) {

    arma::rowvec func_row = func_tab.row(i);

    arma::uvec taxa_w_func_i = arma::find(func_row > 0);
    arma::mat func_contributor_abun = abun_tab.rows(taxa_w_func_i);

    arma::rowvec func_contributor_abun_sample_sums = arma::sum(func_contributor_abun, 0);
    arma::uvec nonzero_samples_i = arma::find(func_contributor_abun_sample_sums > 0);

    List func_info(2);
 
    NumericVector taxa_w_func_i_vec = Rcpp::NumericVector(taxa_w_func_i.begin(), taxa_w_func_i.end());
    NumericVector nonzero_samples_i_vec = Rcpp::NumericVector(nonzero_samples_i.begin(), nonzero_samples_i.end());
    
    func_info[0] = taxa_w_func_i_vec;
    func_info[1] = nonzero_samples_i_vec;
    
    all_func_contrib_abun[i] = func_info;

  }

  return all_func_contrib_abun;

}
