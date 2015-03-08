#ifndef LOGSUMEXP_INCLUDED
#define LOGSUMEXP_INCLUDED

/* uncomment to disable lookup table */
/*
#define LOGSUMEXP_DEBUG
*/

void init_log_sum_exp_lookup();  /* call this first, to initialize lookup table */

long double log_sum_exp (long double a, long double b);  /* returns log(exp(a) + exp(b)) */
long double log_sum_exp_unary (long double x);  /* returns log(1 + exp(-x)) for nonnegative x */
long double log_sum_exp_unary_slow (long double x);  /* does not use lookup table */

#endif /* LOGSUMEXP_INCLUDED */
