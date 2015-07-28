#include <math.h>
#include "gamma.h"

#define nemes_log_gamma_stepdown_threshold 4
#define half_log_2pi 0.918938533204673

double log_gamma (double z) {
  if (z <= 0)
    return INFINITY;

  // if below accuracy threshold, step down, using log(Gamma(z)) = log(Gamma(z+1)) - log(z)
  double offset = 0;
  while (z < nemes_log_gamma_stepdown_threshold) {
    offset -= log(z);
    z += 1;
  }

  // Nemes approximation
  // Nemes, G. "New asymptotic expansion for the Gamma function", Archiv der Mathematik 95 (2): 161–169
  double f = 1 / (15*z*z);
  double lg = offset + (z - 0.5)*log(z) - z + half_log_2pi + (1.25*z)*log(1+f);

  return lg;
}

double log_beta_dist (double x, double a, double b) {
  double log_beta = log_gamma(a) + log_gamma(b) - log_gamma(a+b);
  return (a - 1) * log(x) + (b - 1) * log(1-x) - log_beta;
}

double log_gamma_dist (double x, double a, double b) {
  return a * log(b) + (a - 1) * log(x) - x * b - log_gamma(a);
}
