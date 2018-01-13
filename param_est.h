#ifndef PARAM_EST_H_
#define PARAM_EST_H_

// Stores all data needed to calculate REMLE.
// NewREMLE() provides interface to assign and calculate all members.
struct REMLEstimator{
  int num_i;
  int num_t;
  int num_c;
  double* Y;
  double* covariate_matrix;
  double* covariate_kronecker_matrix;
  double* sum_Y_ind;
  double* sum_Y_covariate_t; 
  double* sum_covariate_effects;
  //Constants in REML function
  double kappa;
  double mu;
  double coefficient;
  double constant;
  double opt_delta;
  double vg;
  double ve;
};

// Calculates REMLE variance components.
//   Returns REMLEstimator struct with params for likelihood
//     Y is an array 
//     covariate matrix is LAPACK_ROW_MAJOR formatted array
struct REMLEstimator* NewREMLE(int num_i, int num_t, int num_c,
                        double* Y, double* covariate_matrix);

// Loads values from file into double array
void LoadValues(char* file_path, double* val_array);

// Calculates inverse of first Kronecker product factor
//   Returns error code from LAPACKE_dgetr[f,i].
int InvertCovariateKroneckerMatrix(struct REMLEstimator* reml_func);

// Calculates Y dependent values for REML function
void CalculateYValues(struct REMLEstimator* reml_func);

// Calculates constants in REML function
void CalculateConstants(struct REMLEstimator* reml_func);

// Solves restricted log-likelihood for given param ratio
double SolveREML(struct REMLEstimator* reml_func, double delta);

// Solves for derivative at given param ratio in likelihood fxn
double DerivativeREML(struct REMLEstimator* reml_func, double delta);

#endif // PARAM_EST_H_
