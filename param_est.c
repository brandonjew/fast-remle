#include <math.h>
#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include "param_est.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

struct REMLEstimator* NewREMLE(int num_i, int num_t, int num_c,
                               double* Y, double* covariate_matrix){
  // TODO: malloc to use heap instead of stack space for large sample sizes.
  //       No need for memset, initialize through first loop
  double covariate_kronecker[num_c*num_c];
  memset(covariate_kronecker, 0.0, sizeof covariate_kronecker);
  double sum_Y_ind[num_i];
  memset(sum_Y_ind, 0.0, sizeof sum_Y_ind);
  double sum_Y_covariate_t[num_c*num_t];
  memset(sum_Y_covariate_t, 0.0, sizeof sum_Y_covariate_t);
  double sum_covariate_effects[num_c];
  memset(sum_covariate_effects, 0.0, sizeof sum_covariate_effects);
  double kappa, mu, coefficient, constant, opt_delta, vg, ve;
  struct REMLEstimator* reml_func = malloc(sizeof *reml_func);
  reml_func->num_i = num_i;
  reml_func->num_t = num_t;
  reml_func->num_c = num_c;
  reml_func->Y = Y;
  reml_func->covariate_matrix = covariate_matrix;
  reml_func->covariate_kronecker_matrix = covariate_kronecker;
  reml_func->sum_Y_ind = sum_Y_ind;
  reml_func->sum_Y_covariate_t = sum_Y_covariate_t;
  reml_func->sum_covariate_effects = sum_covariate_effects;
  InvertCovariateKroneckerMatrix(reml_func);
  CalculateYValues(reml_func);
  CalculateConstants(reml_func);
  return reml_func;
}

void LoadValues(char* file_path, double* val_array){
  int max_val_length = 100;
  FILE* in_file = fopen(file_path, "r");
  char val[max_val_length];
  int curr_char = 0, curr_val = 0;
  char c;
  while((c = fgetc(in_file)) != EOF){
    if (c == '\n' || curr_char == (max_val_length-1)){
      val[curr_char] = '\0';
      val_array[curr_val] = strtod(val, NULL);
      curr_char = 0;
      curr_val++;
    }
    else{
      val[curr_char] = c;
      curr_char++;
    }
  }
  fclose(in_file);
}

// Calculates first Kronecker product of X^T*X and inverts it.
int InvertCovariateKroneckerMatrix(struct REMLEstimator* reml_func){
  int c = reml_func->num_c;
  int n = reml_func->num_i, t = reml_func->num_t;
  int cc = c*c;
  // Kronecker
  double* covariate_kronecker = reml_func->covariate_kronecker_matrix;
  for (int i = 0; i < c; i++){
    for (int j = 0; j < c; j++){
      if (i > j)
        covariate_kronecker[(c*i)+j] = covariate_kronecker[(c*j)+i];
      else{
        for (int k = 0; k < reml_func->num_i; k++){
          covariate_kronecker[(c*i)+j] += 
            reml_func->covariate_matrix[(n*i)+k] *
            reml_func->covariate_matrix[(n*j)+k];
        }
      }
    }
  }
  // Inversion
  int error_handler_f, error_handler_i;
  int pivot_array[c];
  error_handler_f = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, c, c,
                                 reml_func->covariate_kronecker_matrix, c,
                                 pivot_array);
  error_handler_i = LAPACKE_dgetri(LAPACK_ROW_MAJOR, c,
                                 reml_func->covariate_kronecker_matrix, c,
                                 pivot_array);
  if (error_handler_f || error_handler_i)
    perror("Covariates NOT linearly independent");
  return error_handler_i;
}

// Sums over products of covariate values and Y measurments
void CalculateYValues(struct REMLEstimator* reml_func){
  int n = reml_func->num_i;
  int t = reml_func->num_t;
  int c = reml_func->num_c;
  // TODO: iterate through all measurements and calculate sums from paper
}

// Determines a number of assorted constants
void CalculateConstants(struct REMLEstimator* reml_func){
  double mse = 0, kappa = 0;
  int n = reml_func->num_i;
  int t = reml_func->num_t;
  int c = reml_func->num_c;
  for (int i = 0; i < (n*t); i++){
    int t = i / n;
    int ind = i % n;
    double mse_term_1 = 0;
    for (int j = 0; j < c; j++){
      double mse_term_2 = 0;
      for (int k = 0; k < c; k++){
        mse_term_2 += reml_func->covariate_matrix[(n*k)+ind] *
                      reml_func->covariate_kronecker_matrix[(c*k)+j];
      }
      mse_term_1 +=
        reml_func->sum_Y_covariate_t[(j*t)+t] * mse_term_2;
    }
    mse_term_1 = reml_func->Y[i] - mse_term_1;

    mse += (mse_term_1*mse_term_1);
    if (i < n){
      double kappa_term = reml_func->sum_Y_ind[ind];
      for (int j = 0; j < c; j++){
        //TODO: Solve for Kappa 
        kappa_term -= 1;
      }
      kappa += (kappa_term * kappa_term);
    }
  }
  reml_func->kappa = kappa/t;
  reml_func->mu = 0; //TODO: Add formulation for mu from paper
  mse = mse/(n*t);
  reml_func->coefficient = (n-c)/2;
  double p = reml_func->coefficient;
  reml_func->constant = -p*t*(log(M_PI/(p*t)) + 1);
  double k = reml_func->kappa;
  double m = reml_func->mu;
  // TODO: Solve for d (from paper)
  double d = 1;
  reml_func->opt_delta = d;
  if (reml_func->opt_delta < 0){
    printf("Root %f not valid, assuming vg is 0\n", reml_func->opt_delta);
    reml_func->vg = 0;
    reml_func->ve = mse;
  }
  else{
    double R = (k/(t+d)) + (m/d);
    reml_func->vg = R/((n*t) - (c*t));
    reml_func->ve = d * reml_func->vg;
  }
}

int main(int argc, char* argv[]){
  //TODO: Usable interface
  int num_i = atoi(argv[1]), num_t = atoi(argv[2]),
      num_c = argc - 3;
  double Y[num_i*num_t];
  LoadValues(argv[3], Y);
  double covariate_matrix[num_c*num_i];
  for (int i = 0; i < num_i; i++)
    covariate_matrix[i] = 1;
  int curr_cov = 1;
  for (int i = 4; i < argc; i++){
    LoadValues(argv[i], &covariate_matrix[num_i*curr_cov]);
    curr_cov++;
  }
  struct REMLEstimator* remle = NewREMLE(num_i, num_t, num_c, Y,
                                 covariate_matrix);
  free(remle);
  return 0;
}
