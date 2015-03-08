#include <math.h>
#include "logsumexp.h"

#define LOG_SUM_EXP_LOOKUP_MAX 10
#define LOG_SUM_EXP_LOOKUP_PRECISION .0001

#define LOG_SUM_EXP_LOOKUP_ENTRIES (((int) (LOG_SUM_EXP_LOOKUP_MAX / LOG_SUM_EXP_LOOKUP_PRECISION)) + 1)

long double log_sum_exp_lookup[LOG_SUM_EXP_LOOKUP_ENTRIES];

void init_log_sum_exp_lookup() {
  int n;
  long double x;
  for (n = 0; n < LOG_SUM_EXP_LOOKUP_ENTRIES; ++n) {
    x = n * LOG_SUM_EXP_LOOKUP_PRECISION;
    log_sum_exp_lookup[n] = log_sum_exp_unary_slow(x);
  }
}

long double log_sum_exp (long double a, long double b) {
  double min, max, diff;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  return max + log_sum_exp_unary (diff);
}

long double log_sum_exp_unary (long double x) {
#ifdef LOGSUMEXP_DEBUG
  return log_sum_exp_unary_slow(x);
#else /* LOGSUMEXP_DEBUG */
  int n;
  long double dx, f0, f1, df;
  if (x >= LOG_SUM_EXP_LOOKUP_MAX || isnan(x) || isinf(x))
    return 0;
  if (x < 0)
    return -x;
  n = (int) (x / LOG_SUM_EXP_LOOKUP_PRECISION);
  dx = x - (n * LOG_SUM_EXP_LOOKUP_PRECISION);
  f0 = log_sum_exp_lookup[n];
  f1 = log_sum_exp_lookup[n+1];
  df = f1 - f0;
  return f0 + df * (dx / LOG_SUM_EXP_LOOKUP_PRECISION);
#endif /* LOGSUMEXP_DEBUG */
}

long double log_sum_exp_unary_slow (long double x) {
  return log (1. + exp(-x));
}
