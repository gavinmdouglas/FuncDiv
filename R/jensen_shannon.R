jensen_shannon_divergence_FuncPtr <- cppXPtr(
  "double customDist(const arma::mat &A, const arma::mat &B) {
    arma::mat p = A / arma::accu(A);
    arma::mat q = B / arma::accu(B);
    arma::mat m = (p + q) * 0.5;
    double result = 0.5 * arma::accu(p * arma::log(p / m).t()) + 0.5 * arma::accu(q * arma::log(q / m).t());
    return std::isinf(result) ? std::numeric_limits<double>::quiet_NaN()
    : result;
}", depends = c("RcppArmadillo"))