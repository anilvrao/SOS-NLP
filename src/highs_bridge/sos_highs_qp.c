#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "highs/interfaces/highs_c_api.h"

static double sos_bound(double value, double bigbnd, int lower, double inf) {
  if (fabs(value) >= bigbnd) {
    return lower ? -inf : inf;
  }
  return value;
}

void sos_highs_qp_(int* n, int* nclin, int* lda, int* ldh, int* nhess,
                   double* a, double* bl, double* bu, double* cvec,
                   double* h, int* istate, double* x, int* inform,
                   double* obj, int* msglvl, double* ax, double* clamda,
                   double* bigbnd, double* tolfea, double* tolopt) {
  const int num_col = *n;
  const int num_row = *nclin;
  const int a_nz = num_col * num_row;
  const int quiet = *msglvl <= 0;
  double inf = 1.0e30;
  int q_nz = 0;
  int status = 0;
  void* highs = NULL;

  HighsInt* a_start = NULL;
  HighsInt* a_index = NULL;
  double* a_value = NULL;
  HighsInt* q_start = NULL;
  HighsInt* q_index = NULL;
  double* q_value = NULL;
  double* col_lower = NULL;
  double* col_upper = NULL;
  double* row_lower = NULL;
  double* row_upper = NULL;
  double* col_dual = NULL;
  double* row_value = NULL;
  double* row_dual = NULL;

  *inform = 9;
  *obj = 0.0;

  if (num_col <= 0 || num_row < 0) return;

  highs = Highs_create();
  if (!highs) goto cleanup;
  inf = Highs_getInfinity(highs);

  col_lower = (double*)calloc((size_t)num_col, sizeof(double));
  col_upper = (double*)calloc((size_t)num_col, sizeof(double));
  row_lower = (double*)calloc((size_t)(num_row > 0 ? num_row : 1), sizeof(double));
  row_upper = (double*)calloc((size_t)(num_row > 0 ? num_row : 1), sizeof(double));
  col_dual = (double*)calloc((size_t)num_col, sizeof(double));
  row_value = (double*)calloc((size_t)(num_row > 0 ? num_row : 1), sizeof(double));
  row_dual = (double*)calloc((size_t)(num_row > 0 ? num_row : 1), sizeof(double));
  a_start = (HighsInt*)calloc((size_t)(num_col + 1), sizeof(HighsInt));
  a_index = (HighsInt*)calloc((size_t)(a_nz > 0 ? a_nz : 1), sizeof(HighsInt));
  a_value = (double*)calloc((size_t)(a_nz > 0 ? a_nz : 1), sizeof(double));
  q_start = (HighsInt*)calloc((size_t)(num_col + 1), sizeof(HighsInt));

  if (!col_lower || !col_upper || !row_lower || !row_upper || !col_dual ||
      !row_value || !row_dual || !a_start || !a_index || !a_value ||
      !q_start) {
    goto cleanup;
  }

  for (int j = 0; j < num_col; ++j) {
    col_lower[j] = sos_bound(bl[j], *bigbnd, 1, inf);
    col_upper[j] = sos_bound(bu[j], *bigbnd, 0, inf);
  }
  for (int i = 0; i < num_row; ++i) {
    row_lower[i] = sos_bound(bl[num_col + i], *bigbnd, 1, inf);
    row_upper[i] = sos_bound(bu[num_col + i], *bigbnd, 0, inf);
  }

  for (int j = 0; j < num_col; ++j) {
    a_start[j] = (HighsInt)(j * num_row);
    for (int i = 0; i < num_row; ++i) {
      const int p = j * num_row + i;
      a_index[p] = (HighsInt)i;
      a_value[p] = a[j * (*lda) + i];
    }
  }
  a_start[num_col] = (HighsInt)a_nz;

  if (*ldh > 0 && *nhess > 0) {
    for (int j = 0; j < num_col; ++j) {
      q_start[j] = (HighsInt)q_nz;
      for (int i = j; i < num_col; ++i) {
        const double hij = h[j * (*ldh) + i];
        if (hij != 0.0) ++q_nz;
      }
    }
    q_start[num_col] = (HighsInt)q_nz;
    q_index = (HighsInt*)calloc((size_t)(q_nz > 0 ? q_nz : 1), sizeof(HighsInt));
    q_value = (double*)calloc((size_t)(q_nz > 0 ? q_nz : 1), sizeof(double));
    if (!q_index || !q_value) goto cleanup;
    q_nz = 0;
    for (int j = 0; j < num_col; ++j) {
      for (int i = j; i < num_col; ++i) {
        const double hij = h[j * (*ldh) + i];
        if (hij != 0.0) {
          q_index[q_nz] = (HighsInt)i;
          q_value[q_nz] = hij;
          ++q_nz;
        }
      }
    }
  } else {
    q_start[num_col] = 0;
  }

  Highs_setBoolOptionValue(highs, "output_flag", quiet ? 0 : 1);
  if (*tolfea > 0.0) {
    Highs_setDoubleOptionValue(highs, "primal_feasibility_tolerance", *tolfea);
    Highs_setDoubleOptionValue(highs, "dual_feasibility_tolerance", *tolfea);
  }
  if (*tolopt > 0.0) {
    Highs_setDoubleOptionValue(highs, "optimality_tolerance", *tolopt);
  }

  status = Highs_passModel(highs, (HighsInt)num_col, (HighsInt)num_row,
                           (HighsInt)a_nz, (HighsInt)q_nz,
                           kHighsMatrixFormatColwise,
                           kHighsHessianFormatTriangular,
                           kHighsObjSenseMinimize, 0.0, cvec, col_lower,
                           col_upper, row_lower, row_upper, a_start, a_index,
                           a_value, q_start, q_index, q_value, NULL);
  if (status != kHighsStatusOk) goto cleanup;

  status = Highs_run(highs);
  if (status != kHighsStatusOk && status != kHighsStatusWarning) goto cleanup;

  status = Highs_getSolution(highs, x, col_dual, row_value, row_dual);
  if (status != kHighsStatusOk) goto cleanup;

  *obj = Highs_getObjectiveValue(highs);
  for (int i = 0; i < num_row; ++i) ax[i] = row_value[i];
  for (int j = 0; j < num_col; ++j) clamda[j] = col_dual[j];
  for (int i = 0; i < num_row; ++i) clamda[num_col + i] = row_dual[i];
  for (int i = 0; i < num_col + num_row; ++i) istate[i] = 0;

  switch (Highs_getModelStatus(highs)) {
    case kHighsModelStatusOptimal:
      *inform = 0;
      break;
    case kHighsModelStatusInfeasible:
      *inform = 2;
      break;
    case kHighsModelStatusUnbounded:
    case kHighsModelStatusUnboundedOrInfeasible:
      *inform = 3;
      break;
    default:
      *inform = 9;
      break;
  }

cleanup:
  if (highs) Highs_destroy(highs);
  free(a_start);
  free(a_index);
  free(a_value);
  free(q_start);
  free(q_index);
  free(q_value);
  free(col_lower);
  free(col_upper);
  free(row_lower);
  free(row_upper);
  free(col_dual);
  free(row_value);
  free(row_dual);
}
