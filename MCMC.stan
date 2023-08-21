functions{
// BM kernel
real kern_bm(
  real x_i,
  real y_i){
    real result=0.5 *(abs(x_i)+abs(y_i)-abs(x_i-y_i));
    return result;
  }
// gram matrix
matrix sq_bm_K_gram(int N,vector x){
    matrix[N,N] K_gram;
    matrix[N,N] A;
    matrix[N,N] result;
    matrix[N,N] C = diag_matrix(rep_vector(1.0, N)) - (1.0 / N) * rep_matrix(1.0, N, N);
    for (i in 1:N) {
      for (j in 1:N) {
        K_gram[i, j] = kern_bm(x[i], x[j]);
        }
     }
    A = C * K_gram * C;
    result = A*A;
    return result;
}
// kernel vector
matrix sq_bm_k_vec(int N, vector x, vector x_new){
    matrix[N,N] K_gram;
    matrix[N,N] C = diag_matrix(rep_vector(1.0, N)) - (1.0 / N) * rep_matrix(1, N, N);
    matrix[N,N] K_cent;
    int M = num_elements(x_new);
    matrix[N,M] k_vec;
    matrix[N,M] centered;
    matrix[N,M] result;
    for (i in 1:N) {
      for (j in 1:N) {
        K_gram[i, j] = kern_bm(x[i], x[j]);
        }
     }
    K_cent = C' * K_gram * C;
    for (i in 1:N) {
      for (j in 1:M) {
        k_vec[i, j] = kern_bm(x[i], x_new[j]);
        }
     }
    {
      vector[N] first_element = K_gram * rep_vector(1.0,N);
      row_vector[M] second_element = rep_row_vector(1.0,N)*k_vec;
      real third_element = rep_row_vector(1.0,N) * K_gram * rep_vector(1.0,N);
      centered = k_vec - (1.0/N * first_element * rep_row_vector(1.0,M))-  (1.0/N * rep_vector(1.0,N)*second_element)+ (1.0/(square(N)) * third_element * rep_matrix(1, N, M));
      result = K_cent * centered;
    }
    return result;
}
matrix sq_bm_k_star(int N, vector x, vector x_new){
    matrix[N,N] K_gram;
    matrix[N,N] C = diag_matrix(rep_vector(1.0, N)) - (1.0 / N) * rep_matrix(1, N, N);
    matrix[N,N] K_cent;
    int M = num_elements(x_new);
    matrix[N,M] k_vec;
    matrix[N,M] centered;
    matrix[M,M] result;
    for (i in 1:N) {
      for (j in 1:N) {
        K_gram[i, j] = kern_bm(x[i], x[j]);
        }
     }
    K_cent = C' * K_gram * C;
    for (i in 1:N) {
      for (j in 1:M) {
        k_vec[i, j] = kern_bm(x[i], x_new[j]);
        }
     }
    {
      vector[N] first_element = K_gram * rep_vector(1.0,N);
      row_vector[M] second_element = rep_row_vector(1.0,N)*k_vec;
      real third_element = rep_row_vector(1.0,N) * K_gram * rep_vector(1.0,N);
      centered = k_vec - (1.0/N * first_element * rep_row_vector(1.0,M))-  (1.0/N * rep_vector(1.0,N)*second_element)+ (1.0/(square(N)) * third_element * rep_matrix(1, N, M));
      result = centered' * centered;
    }
    return result;
}
}
data {
  int<lower=1> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
  vector[N] x4;
  vector[N] x5;
  vector[N] x6;
  vector[N] y;
  int<lower=1> M;
  matrix[M,1] x_new1;
  matrix[M,1] x_new2;
  matrix[M,1] x_new3;
  matrix[M,1] x_new4;
  matrix[M,1] x_new5;
  matrix[M,1] x_new6;
}
transformed data {
  vector[N] mu = rep_vector(0, N);
  matrix[N,N] K_gram1 = sq_bm_K_gram(N,x1);
  matrix[N,N] K_gram2 = sq_bm_K_gram(N,x2);
  matrix[N,N] K_gram3 = sq_bm_K_gram(N,x3);
  matrix[N,N] K_gram4 = sq_bm_K_gram(N,x4);
  matrix[N,N] K_gram5 = sq_bm_K_gram(N,x5);
  matrix[N,N] K_gram6 = sq_bm_K_gram(N,x6);
  matrix[N,M] k_vec1 = sq_bm_k_vec(N, x1, x_new1[,1]);
  matrix[N,M] k_vec2 = sq_bm_k_vec(N, x2, x_new2[,1]);
  matrix[N,M] k_vec3 = sq_bm_k_vec(N, x3, x_new3[,1]);
  matrix[N,M] k_vec4 = sq_bm_k_vec(N, x4, x_new4[,1]);
  matrix[N,M] k_vec5 = sq_bm_k_vec(N, x5, x_new5[,1]);
  matrix[N,M] k_vec6 = sq_bm_k_vec(N, x6, x_new6[,1]);
  matrix[M,M] k_star1 = sq_bm_k_star(N,x1, x_new1[,1]);
  matrix[M,M] k_star2 = sq_bm_k_star(N,x2, x_new2[,1]);
  matrix[M,M] k_star3 = sq_bm_k_star(N,x3, x_new3[,1]);
  matrix[M,M] k_star4 = sq_bm_k_star(N,x4, x_new4[,1]);
  matrix[M,M] k_star5 = sq_bm_k_star(N,x5, x_new5[,1]);
  matrix[M,M] k_star6 = sq_bm_k_star(N,x6, x_new6[,1]);
}
parameters {
  real<lower=0> sigma;
  real<lower=0> lambda_1;
  real<lower=0> lambda_2;
  //real<lower=0> lambda_3;
  //real<lower=0> lambda_4;
  //real<lower=0> lambda_5;
  //real<lower=0> lambda_6;
  real<lower=0> alpha0;
}
model{
  matrix[N,N] K = square(alpha0)*(rep_matrix(1,N,N)+lambda_1*K_gram1+lambda_2*K_gram2);//+lambda_3*K_gram3+lambda_4*K_gram4+lambda_5*K_gram5);
  // diagonal elements
  matrix[N,N] L;
  for (n in 1:N)
    K[n, n] = K[n, n] + square(sigma);
  L=cholesky_decompose(K);
  
  //prior
  // sigma ~ normal(0, 1);
  // lambda_1 ~ normal(0, 1);
  // lambda_2 ~ normal(0, 1);
  // alpha0 ~ normal(0, 1);
  target += std_normal_lpdf(alpha0);
  target += std_normal_lpdf(lambda_1);
  target += std_normal_lpdf(lambda_2);
  //target += std_normal_lpdf(lambda_3);
  //target += std_normal_lpdf(lambda_4);
  //target += std_normal_lpdf(lambda_5);
  //target += std_normal_lpdf(lambda_6);
  target += std_normal_lpdf(sigma);
  //likelihood
  //target += -0.5 * y'*inverse(K)*y -sum(log(diagonal(L)))-N/2.0*log(2.0*pi());//y'*mdivide_left_tri(capital_sigma,y)
  //target += multi_normal_lpdf(y | mu, capital_sigma); // they should be the same
  //y ~ multi_normal_cholesky(mu, L);
  target += multi_normal_cholesky_lpdf(y | mu, L);
}
generated quantities{
  real mloglik;
  vector[M] mu_predicted;
  matrix[M,M] var_predicted;
  // vector[N] f_zero;
  // vector[M] f_one;
  // vector[M] f_two;
  // vector[M] f_three;
  // vector[N] f_four;
  // vector[N] f_five;
  // vector[N] f;
  {
  matrix[N,N] K_eval;
  matrix[N,N] L_eval;
  matrix[N,M] k_vec_eval;
  matrix[M,M] k_star_eval;
  vector[N] alpha_eval;
  matrix[N,M] v;
  K_eval = square(alpha0)*(rep_matrix(1,N,N)+lambda_1*K_gram1+lambda_2*K_gram2);//+lambda_3*K_gram3+lambda_4*K_gram4+lambda_5*K_gram5);
  for (n in 1:N)
    K_eval[n, n] = K_eval[n, n] + square(sigma);
  //matrix[N,N] capital_sigma_eval = (K_eval + square(sigma) * diag_matrix(rep_vector(1, N)));
  L_eval=cholesky_decompose(K_eval);
  k_vec_eval = square(alpha0)*(rep_matrix(1,N,M)+lambda_1*k_vec1+lambda_2*k_vec2);//+lambda_3*k_vec3+lambda_4*k_vec4+lambda_5*k_vec5);
  k_star_eval = square(alpha0)*(rep_matrix(1,M,M)+lambda_1*k_star1+lambda_2*k_star2);//+lambda_3*k_star3+lambda_4*k_star4+lambda_5*k_star5);
  //Rasmussen
  alpha_eval=L_eval'\(L_eval\y);
  v = L_eval\k_vec_eval;
  mu_predicted = k_vec_eval'*alpha_eval;
  var_predicted = k_star_eval - v'*v;
  mloglik = -0.5 * y'*alpha_eval - sum(log(diagonal(L_eval)))- N/2.0*log(2.0*pi());
  // f_zero = square(alpha0)*sum(alpha_eval)*rep_vector(1,N);
  // f_one =square(alpha0)*lambda_1*K_gram1*alpha_eval;//K_gram1*alpha_eval
  // f_two =square(alpha0)*lambda_2*K_gram2*alpha_eval;
  // f_three =square(alpha0)*lambda_3*K_gram3*alpha_eval;
  // f_four =square(alpha0)*lambda_4*K_gram4*alpha_eval;
  // f_five =square(alpha0)*lambda_5*K_gram5*alpha_eval;
  // f = square(alpha0)*(rep_matrix(1,N,N)+lambda_1*K_gram1+lambda_2*K_gram2+lambda_3*K_gram3+lambda_4*K_gram4+lambda_5*K_gram5)*alpha_eval;
  }
}






