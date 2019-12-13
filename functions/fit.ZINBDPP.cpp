// #include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

static double norm_rs(double a, double b);
static double half_norm_rs(double a, double b);
static double unif_rs(double a, double b);
static double exp_rs(double a, double b);
static double rnorm_trunc(double mu, double sigma, double lower, double upper);

// [[Rcpp::export]]
Rcpp::List fitZINBDPP(IntegerMatrix Y, 
                      IntegerVector z, 
                      NumericVector s, 
                      int iter, 
                      bool DPP, 
                      IntegerMatrix S, 
                      bool aggregate, 
                      double b, 
                      double h, 
                      bool MRF, 
                      IntegerMatrix G,
                      double phi_low = 1,
                      int mincount = 2) {
  // Read data dimensionality
  int n = Y.nrow();
  int p = Y.ncol();
  int K = max(z) + 1;
  int p_ext = S.nrow();
  int pp = S.ncol();
  if(!aggregate) 
  {
    p_ext = 0;
  }
  
  // Extend data matrix Y if phylogenetic tree provided
  IntegerMatrix Y_ext(n, p + p_ext);
  int i, j, jj;
  int p_b = p;
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < p; j++)
    {
      Y_ext(i, j) = Y(i, j);
    }
  }
  if(aggregate)
  {
    for(j = 0; j < p_ext; j++)
    {
      for(i = 0; i < n; i++)
      {
        Y_ext(i, p + j) = 0;
        for(jj = 0; jj < pp; jj++)
        {
          if(S(j, jj) == 0)
          {
            break;
          }
          else
          {
            Y_ext(i, p + j) = Y_ext(i, p + j) + Y_ext(i, S(j, jj) - 1);
          }
        }
      }
    }
    p = p + p_ext;
  }
  
  // Set hyperparameters
  double a_pi = 1;
  double b_pi = 1;
  double a_phi = 0.001; 
  double b_phi = 0.001; 
  double a_omega = 1.0;
  double b_omega = 1.0;
  double dd = -2.2;
  double ff = 0.5;
  
  
  double a_0 = 2.0;
  double a = 2.0;
  double b_0 = b;
  // double b = 1.0;
  double h_0 = h;
  //double h = 1000;
  
  // Set algorithm settings
  int burn = iter*0.5;
  int E = 20;
  int M = n*0.5;
  double tau_phi = 1.0;
  double tau_alpha = 1.0;
  double tau_s = 1.0;
  
  // Set temporary variables
  int it, ii, k, e, gamma_temp, count_2, m, mm, count_temp, count_temp_2;
  int count = 0;
  int gamma_sum = 0;
  int H_sum = 0;
  int p_valid = 0;
  double pi, alpha_temp, max_temp, sum_temp, hastings, phi_temp, temp, temp_2, temp_4, temp_5, s_temp, loga_lwr, loga_upp, min_Y;
  double accept_gamma = 0.0;
  double accept_phi = 0.0;
  double accept_alpha = 0.0;
  
  IntegerVector gamma(p);
  IntegerVector flag(p);
  IntegerVector gamma_sum_store(iter);
  IntegerVector H_sum_store(iter);
  IntegerVector nk(K);
  NumericVector phi(p);
  NumericVector phi_mean(p);
  NumericVector gamma_ppi(p);
  NumericVector pi_store(iter);
  NumericVector prob_temp(2);
  IntegerMatrix H(n, p);
  NumericMatrix H_ppi(n, p);
  NumericMatrix A(n, p);
  NumericMatrix A_mean(n, p);
  NumericVector temp_3(n);
  NumericMatrix logafc(iter - burn, p);
  NumericMatrix s_store(iter - burn, n);
  
  // Load settings for Dirichlet process prior
  double c_s = 0.0;
  double a_m = 1.0;
  double b_m = 1.0;
  double a_t = 1.0;
  double b_t = 1.0;
  double tau_eta = 1.0;
  double sigma_s = 1.0;
  double logs_infty;
  IntegerVector nu_s(n);
  IntegerVector epsilon_s(n);
  IntegerVector state_s(M);
  NumericVector prob_temp_s(M);
  NumericVector psi_s(M);
  NumericVector t_s(M);
  NumericVector eta_s(M);
  NumericVector s_mean(n);
  
  // Initialization
  for(j = 0; j < p; j++)
  {
    for(k = 0; k < K; k++)
    {
      count_2 = 0;
      for(i = 0; i < n; i++)
      {
        if(z(i) == k && Y_ext(i, j) != 0)
        {
          count_2++;
        }
      }
      if(count_2 < mincount)
      {
        flag(j) = 1;
        break;
      }
    }
    if(flag(j) == 0)
    {
      p_valid++;
      gamma(j) = rbinom(1, 1, 0.05)(0);
      phi(j) = 10;
      phi_mean(j) = 0;
      gamma_ppi(j) = 0;
      if(gamma(j) == 1)
      {
        gamma_sum++;
      }
      for(i = 0; i < n; i++)
      {
        A(i, j) = 1.0;
        A_mean(i, j) = 0;
        if(Y_ext(i, j) == 0)
        {
          H(i, j) = rbinom(1, 1, 0.5)(0);
        }
        if (H(i, j) == 1) {
          H_sum++;
        }
      }
    }
    else
    {
      gamma(j) = -1;
      phi(j) = -1;
      phi_mean(j) = -1;
      gamma_ppi(j) = 0;
      for(i = 0; i < n; i++)
      {
        A(i, j) = -1;
        A_mean(i, j) = -1;
        H(i, j) = -1;
      }
    }
  }
  
  for(m = 0; m < M; m++)
  {
    psi_s(m) = 1.0/M;
    t_s(m) = 0.5;
    eta_s(m) = 0;
    prob_temp_s(m) = 1.0/M;
    state_s(m) = m;
  }
  min_Y = max(Y_ext);
  for(i = 0; i < n; i++)
  {
    temp = RcppArmadillo::sample(state_s, 1, true, prob_temp_s)(0);
    nu_s(i) = temp;
    epsilon_s(i) = rbinom(1, 1, 0.5)(0);
    s_mean(i) = 0;
    for(j = 0; j < p_b; j++)
    {
      temp_3(i) = temp_3(i) + Y_ext(i, j);
      if(Y_ext(i, j) < min_Y && Y_ext(i, j) != 0)
      {
        min_Y = Y_ext(i, j);
      }
    }
  }
  logs_infty = log(max(temp_3)/min(temp_3));
  loga_upp = log(max(Y_ext)/exp(c_s - logs_infty));
  loga_lwr = log(min_Y/exp(c_s + logs_infty));
  
  // MCMC
  for(it = 0; it < iter; it++)
  {
    // Update pi
    pi = rbeta(1, a_pi + H_sum, b_pi + n*p - H_sum)(0);
    
    // Update H
    H_sum = 0;
    for(j = 0; j < p; j++)
    {
      if(flag(j) == 0)
      {
        for(i = 0; i < n; i++)
        {
          if(Y_ext(i, j) == 0)
          {
            prob_temp(0) = phi(j)*(log(phi(j)) - log(s(i)*A(i, j) + phi(j))) + log(1 - pi);
            prob_temp(1) = log(pi);
            max_temp = max(prob_temp);
            prob_temp(0) = prob_temp(0) - max_temp;
            prob_temp(1) = prob_temp(1) - max_temp;
            prob_temp(0) = exp(prob_temp(0));
            prob_temp(1) = exp(prob_temp(1));
            sum_temp = prob_temp(0) + prob_temp(1);
            prob_temp(0) = prob_temp(0)/sum_temp;
            prob_temp(1) = prob_temp(1)/sum_temp;
            H(i, j) = rbinom(1, 1, prob_temp(1))(0);
            if (H(i, j) == 1) {
              H_sum++;
            }
          }
        }
      }
    }
    
    // Update phi
    for(j = 0; j < p; j++)
    {
      if(flag(j) == 0)
      {
        count_2 = 0;
        do {
          phi_temp = rgamma(1, phi(j)*phi(j)/tau_phi, tau_phi/phi(j))(0);
          count_2++;
        } while (phi_temp < phi_low && count_2 < 1000);
        if (count_2 == 1000)
        {
          phi_temp = 10;
        }
        hastings = 0;
        for(i = 0; i < n; i++)
        {
          if(H(i, j) == 0) {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + Y_ext(i, j)) - (phi_temp + Y_ext(i, j))*log(phi_temp + s(i)*A(i, j));
            hastings = hastings - (phi(j)*log(phi(j)) - lgamma(phi(j)) + lgamma(phi(j) + Y_ext(i, j)) - (phi(j) + Y_ext(i, j))*log(phi(j) + s(i)*A(i, j)));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi(j)) - b_phi*phi(j));
        if(hastings >= log(double(rand()%10001)/10000))
        {
          phi(j) = phi_temp;
          if (it > burn) {
            accept_phi++;
          }
        }
      }
    }
    
    // Update A
    for(j = 0; j < p; j++)
    {
      if(flag(j) == 0)
      {
        for(i = 0; i < n; i++)
        {
          alpha_temp = exp(rnorm_trunc(log(A(i, j)), tau_alpha, loga_lwr, loga_upp));
          hastings = 0;
          if(H(i, j) == 0) {
            hastings = hastings + Y_ext(i, j)*log(s(i)*alpha_temp) - (phi(j) + Y_ext(i, j))*log(phi(j) + s(i)*alpha_temp);
            hastings = hastings - (Y_ext(i, j)*log(s(i)*A(i, j)) - (phi(j) + Y_ext(i, j))*log(phi(j) + s(i)*A(i, j)));
          }
          if(gamma(j) == 0)
          {
            temp = 0;
            temp_2 = 0;
            count_temp = 1;
            for(ii = 0; ii < n; ii++)
            {
              if(ii != i && H(ii, j) == 0)
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
                count_temp++;
              }
            }
            temp_4 = temp + log(alpha_temp)*log(alpha_temp);
            temp_5 = temp_2 + log(alpha_temp);
            temp_4 = temp_4 - temp_5*temp_5/(count_temp + 1/h_0);
            hastings = hastings - (a_0 + count_temp*0.5)*log(b_0 + temp_4*0.5);
            temp_4 = temp + log(A(i, j))*log(A(i, j));
            temp_5 = temp_2 + log(A(i, j));
            temp_4 = temp_4 - temp_5*temp_5/(count_temp + 1/h_0);
            hastings = hastings + (a_0 + count_temp*0.5)*log(b_0 + temp_4*0.5);
          }
          else
          {
            temp = 0;
            temp_2 = 0;
            count_temp = 1;
            for(ii = 0; ii < n; ii++)
            {
              if(ii != i && H(ii, j) == 0 && z[ii] == z[i])
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
                count_temp++;
              }
            }
            temp_4 = temp + log(alpha_temp)*log(alpha_temp);
            temp_5 = temp_2 + log(alpha_temp);
            temp_4 = temp_4 - temp_5*temp_5/(count_temp + 1/h);
            hastings = hastings - (a + count_temp*0.5)*log(b + temp_4*0.5);
            temp_4 = temp + log(A(i, j))*log(A(i, j));
            temp_5 = temp_2 + log(A(i, j));
            temp_4 = temp_4 - temp_5*temp_5/(count_temp + 1/h);
            hastings = hastings + (a + count_temp*0.5)*log(b + temp_4*0.5);
          }
          if (hastings >= log(double(rand()%10001)/10000))
          {
            A(i, j) = alpha_temp;
            if (it > burn) {
              accept_alpha++;
            }
          }
        }
      }
    }
    
    // Update gamma
    if (K > 1)
    {
      for(e = 0; e < E; e++)
      {
        j = rand()%p;
        if(flag(j) == 0)
        {
          gamma_temp = 1 - gamma(j);
          if(gamma_temp == 0) // Delete
          {
            if(MRF)
            {
              hastings = -dd;
              for(jj = 0; jj < p; jj++)
              {
                if(G(j, jj) == 1)
                {
                  hastings = hastings - ff;
                }
              }
            }
            else 
            {
              hastings = log(b_omega) - log(a_omega);
            }
            temp = 0;
            temp_2 = 0;
            count_temp = 0;
            for(ii = 0; ii < n; ii++)
            {
              if(H(ii, j) == 0)
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
                count_temp++;
              }
            }
            temp = temp - temp_2*temp_2/(count_temp + 1/h_0);
            //hastings = hastings + (-log(count_temp*h_0 + 1)*0.5 + lgamma(a_0 + count_temp*0.5) - lgamma(a_0) + a_0*log(b_0(j)) - (a_0 + count_temp*0.5)*log(b_0(j) + temp*0.5));
            hastings = hastings + (-log(count_temp*h_0 + 1)*0.5 + lgamma(a_0 + count_temp*0.5) - lgamma(a_0) + a_0*log(b_0) - (a_0 + count_temp*0.5)*log(b_0 + temp*0.5));
            for(k = 0; k < K; k++)
            {
              temp = 0;
              temp_2 = 0;
              count_temp = 0;
              for(ii = 0; ii < n; ii++)
              {
                if(z(ii) == k && H(ii, j) == 0)
                {
                  temp = temp + log(A(ii, j))*log(A(ii, j));
                  temp_2 = temp_2 + log(A(ii, j));
                  count_temp++;
                }
              }
              temp = temp - temp_2*temp_2/(count_temp + 1/h);
              //hastings = hastings - (-log(count_temp*h + 1)*0.5 + lgamma(a + count_temp*0.5) - lgamma(a) + a*log(b(k, j)) - (a + count_temp*0.5)*log(b(k, j) + temp*0.5));
              hastings = hastings - (-log(count_temp*h + 1)*0.5 + lgamma(a + count_temp*0.5) - lgamma(a) + a*log(b) - (a + count_temp*0.5)*log(b + temp*0.5));
            }
          }
          else // Add
          {
            if(MRF)
            {
              hastings = dd;
              for(jj = 0; jj < p; jj++)
              {
                if(G(j, jj) == 1)
                {
                  hastings = hastings + ff;
                }
              }
            }
            else 
            {
              hastings = log(a_omega) - log(b_omega);
            }
            for(k = 0; k < K; k++)
            {
              temp = 0;
              temp_2 = 0;
              count_temp = 0;
              for(ii = 0; ii < n; ii++)
              {
                if(z(ii) == k && H(ii, j) == 0)
                {
                  temp = temp + log(A(ii, j))*log(A(ii, j));
                  temp_2 = temp_2 + log(A(ii, j));
                  count_temp++;
                }
              }
              temp = temp - temp_2*temp_2/(count_temp + 1/h);
              //hastings = hastings + (-log(count_temp*h + 1)*0.5 + lgamma(a + count_temp*0.5) - lgamma(a) + a*log(b(k, j)) - (a + count_temp*0.5)*log(b(k, j) + temp*0.5));
              hastings = hastings + (-log(count_temp*h + 1)*0.5 + lgamma(a + count_temp*0.5) - lgamma(a) + a*log(b) - (a + count_temp*0.5)*log(b + temp*0.5));
            }
            temp = 0;
            temp_2 = 0;
            count_temp = 0;
            for(ii = 0; ii < n; ii++)
            {
              if(H(ii, j) == 0)
              {
                temp = temp + log(A(ii, j))*log(A(ii, j));
                temp_2 = temp_2 + log(A(ii, j));
                count_temp++;
              }
            }
            temp = temp - temp_2*temp_2/(count_temp + 1/h_0);
            //hastings = hastings - (-log(count_temp*h_0 + 1)*0.5 + lgamma(a_0 + count_temp*0.5) - lgamma(a_0) + a_0*log(b_0(j)) - (a_0 + count_temp*0.5)*log(b_0(j) + temp*0.5));
            hastings = hastings - (-log(count_temp*h_0 + 1)*0.5 + lgamma(a_0 + count_temp*0.5) - lgamma(a_0) + a_0*log(b_0) - (a_0 + count_temp*0.5)*log(b_0 + temp*0.5));
          }
          if (hastings >= log(double(rand()%10001)/10000))
          {
            gamma(j) = gamma_temp;
            if(gamma_temp == 1) // Add
            {
              gamma_sum++;
            }
            else // Delete
            {
              gamma_sum--;
            }
            if(it > burn) {
              accept_gamma++;
            }
          }
        }
      }
    }
    
    // Update s
    if(DPP)
    {
      // Update s
      for(i = 0; i < n; i++)
      {
        s_temp = exp(rnorm_trunc(log(s(i)), tau_s, c_s - logs_infty, c_s + logs_infty));
        hastings = 0;
        for(j = 0; j < p_b; j++)
        {
          if(flag(j) == 0 && H(i, j) == 0)
          //if(H(i, j) == 0)
          {
            hastings = hastings + Y_ext(i, j)*log(s_temp) - (Y_ext(i, j) + phi(j))*log(s_temp*A(i, j) + phi(j));
            hastings = hastings - (Y_ext(i, j)*log(s(i)) - (Y_ext(i, j) + phi(j))*log(s(i)*A(i, j) + phi(j)));
          }
        }
        if (hastings >= log(double(rand()%10001)/10000))
        {
          s(i) = s_temp;
        }
      }
      // Update nu_s
      for(i = 0; i < n; i++)
      {
        for(m = 0; m < M; m++)
        {
          prob_temp_s(m) = log(psi_s(m)) - (log(s(i)) - epsilon_s(i)*eta_s(m) - (1 - epsilon_s(i))*(c_s - t_s(m)*eta_s(m))/(1 - t_s(m)))*(log(s(i)) - epsilon_s(i)*eta_s(m) - (1 - epsilon_s(i))*(c_s - t_s(m)*eta_s(m))/(1 - t_s(m)))*0.5/sigma_s/sigma_s;
        }
        sum_temp = 0;
        max_temp = max(prob_temp_s);
        for(m = 0; m < M; m++)
        {
          prob_temp_s(m) = prob_temp_s(m) - max_temp;
          prob_temp_s(m) = exp(prob_temp_s(m));
          sum_temp = sum_temp + prob_temp_s(m);
        }
        for(m = 0; m < M; m++)
        {
          prob_temp_s(m) = prob_temp_s(m)/sum_temp;
        }
        temp = RcppArmadillo::sample(state_s, 1, true, prob_temp_s)(0);
        nu_s(i) = temp;
      }
      // Update psi_s
      for(m = 0; m < M; m++)
      {
        temp = 0;
        temp_2 = 0;
        for(i = 0; i < n; i++)
        {
          if(nu_s(i) == m)
          {
            temp++;
          }
          if(nu_s(i) > m)
          {
            temp_2++;
          }
        }
        prob_temp_s(m) = rbeta(1, a_m + temp, b_m + temp_2)(0);
      }
      psi_s(0) = prob_temp_s(0);
      for(m = 1; m < M; m++)
      {
        psi_s(m) = prob_temp_s(m);
        for(mm = 0; mm < m; mm++)
        {
          psi_s(m) = psi_s(m)*(1 - prob_temp_s(mm));
        }
      }
      // Update epsilon_s
      for(i = 0; i < n; i++)
      {
        prob_temp(0) = log(1 - t_s(nu_s(i))) - (log(s(i)) - (c_s - t_s(nu_s(i))*eta_s(nu_s(i)))/(1 - t_s(nu_s(i))))*(log(s(i)) - (c_s - t_s(nu_s(i))*eta_s(nu_s(i)))/(1 - t_s(nu_s(i))))*0.5/sigma_s/sigma_s;
        prob_temp(1) = log(t_s(nu_s(i))) - (log(s(i)) - eta_s(nu_s(i)))*(log(s(i)) - eta_s(nu_s(i)))*0.5/sigma_s/sigma_s;
        max_temp = max(prob_temp);
        prob_temp(0) = prob_temp(0) - max_temp;
        prob_temp(1) = prob_temp(1) - max_temp;
        prob_temp(0) = exp(prob_temp(0));
        prob_temp(1) = exp(prob_temp(1));
        sum_temp = prob_temp(0) + prob_temp(1);
        prob_temp(0) = prob_temp(0)/sum_temp;
        prob_temp(1) = prob_temp(1)/sum_temp;
        epsilon_s(i) = rbinom(1, 1, prob_temp(1))(0);
      }
      // Update t_s
      for(m = 0; m < M; m++)
      {
        temp = 0;
        temp_2 = 0;
        for(i = 0; i < n; i++)
        {
          if(nu_s(i) == m)
          {
            if(epsilon_s(i) == 1)
            {
              temp++;
            }
            else{
              temp_2++;
            }
          }
        }
        t_s(m) = rbeta(1, a_t + temp, b_t + temp_2)(0);
      }
      // Update eta_s
      for(m = 0; m < M; m++)
      {
        temp = 0;
        temp_2 = 0;
        for(i = 0; i < n; i++)
        {
          if(nu_s(i) == m)
          {
            if(epsilon_s(i) == 1)
            {
              temp = temp + log(s(i));
              temp_2++;
            }
            else
            {
              temp = temp - t_s(nu_s(i))/(1 - t_s(nu_s(i)))*(log(s(i)) - c_s/(1 - t_s(nu_s(i))));
              temp_2 = temp_2 + t_s(nu_s(i))*t_s(nu_s(i))/(1 - t_s(nu_s(i)))/(1 - t_s(nu_s(i)));
            }
          }
        }
        eta_s(m) = rnorm(1, temp/sigma_s/sigma_s/(temp_2/sigma_s/sigma_s + 1/tau_eta/tau_eta), sqrt(1/(temp_2/sigma_s/sigma_s + 1/tau_eta/tau_eta)))(0);
      }
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    gamma_sum_store(it) = gamma_sum;
    H_sum_store(it) = H_sum;
    pi_store(it) = pi;
    if(it >= burn) {
      for(j = 0; j < p; j++)
      {
        if(flag(j) == 0)
        {
          gamma_ppi(j) = gamma_ppi(j) + gamma(j);
          phi_mean(j) = phi_mean(j) + phi(j);
          temp = 0;
          temp_2 = 0;
          count_temp = 0;
          count_temp_2 = 0;
          for(i = 0; i < n; i++)
          {
            H_ppi(i, j) = H_ppi(i, j) + H(i, j);
            A_mean(i, j) = A_mean(i, j) + A(i, j);
            if(H(i, j) == 0)
            {
              if(z(i) == 0)
              {
                temp = temp + A(i, j);
                count_temp++;
              }
              else if(z(i) == 1)
              {
                temp_2 = temp_2 + A(i, j);
                count_temp_2++;
              }
            }
          }
          logafc(it - burn, j) = log(temp_2) - log(count_temp_2) - log(temp) + log(count_temp);
        }
        else
        {
          logafc(it - burn, j) = 0;
        }
      }
      for(i = 0; i < n; i++)
      {
        s_mean(i) = s_mean(i) + s(i);
        s_store(it - burn, i) = s(i);
      }
    }
  }
  
  accept_gamma = accept_gamma/E/(iter - burn);
  accept_phi = accept_phi/p_valid/(iter - burn);
  accept_alpha = accept_alpha/n/p_valid/(iter - burn);
  for(j = 0; j < p; j++)
  {
    if(flag(j) == 0)
    {
      gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
      phi_mean(j) = phi_mean(j)/(iter - burn);
      for(i = 0; i < n; i++)
      {
        H_ppi(i, j) = H_ppi(i, j)/(iter - burn);
        A_mean(i, j) = A_mean(i, j)/(iter - burn);
      }
    }
  }
  for(i = 0; i < n; i++)
  {
    s_mean(i) = s_mean(i)/(iter - burn);
  }
  
  return Rcpp::List::create(Rcpp::Named("fold.change") = logafc, 
                            Rcpp::Named("remove.idx") = flag, 
                            Rcpp::Named("size.factor") = s_store, 
                            Rcpp::Named("alpha.matrix") = A_mean, 
                            Rcpp::Named("phi") = phi_mean,
                            Rcpp::Named("pi") = pi_store,
                            Rcpp::Named("H.ppi") = H_ppi, 
                            Rcpp::Named("gamma.ppi") = gamma_ppi, 
                            Rcpp::Named("gamma.accept.rate") = accept_gamma, 
                            Rcpp::Named("phi.accept.rate") = accept_phi, 
                            Rcpp::Named("alpha.accept.rate") = accept_alpha);
}

double rnorm_trunc(double mu, double sigma, double lower, double upper)
{
  int change;
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;
  
  change = 0;
  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;
  
  // First scenario
  if( (a == R_NegInf)||(b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      change = 1;
      a = -b;
      b = R_PosInf;
    }
    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    if(change) z = -z;
  }
  
  // Second scenario
  else if((a*b) <= 0.0)
  {
    // The two possibilities for this scenario
    if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
    {
      z = norm_rs(a, b);
    }
    else z = unif_rs(a,b);
  }
  
  // Third scenario
  else
  {
    if(b < 0)
    {
      tmp = b; b = -a; a = -tmp; change = 1;
    }
    
    lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
    if(lograt <= logt2)
    {
      z = unif_rs(a,b);
    }
    else if((lograt > logt1)&&(a < t3))
    {
      z = half_norm_rs(a,b);
    }
    else
    {
      z = exp_rs(a,b);
    }
    if(change)
    {
      z = -z;
    }
  }
  double output;
  output = sigma*z + mu;
  return (output);
}

double exp_rs(double a, double b)
{
  double  z, u, rate;
  rate = 1/a;
  
  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a))
  {
    z = R::rexp(rate);
  }
  u = R::runif(0.0, 1.0);
  
  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a))
    {
      z = R::rexp(rate);
    }
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}

double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;
  
  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) 
  {
    xstar = 0.0;
  }
  else 
  {
    xstar = a;
  }
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
  
  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while(logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

double half_norm_rs(double a, double b)
{
  double x;
  x = fabs(norm_rand());
  while((x<a)||(x>b)) 
  {
    x = fabs(norm_rand());
  }
  return x;
}

double norm_rs(double a, double b)
{
  double x;
  x = Rf_rnorm(0.0, 1.0);
  while((x < a)||(x > b)) 
  {
    x = norm_rand();
  }
  return x;
}
