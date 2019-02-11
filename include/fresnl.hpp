#ifndef ASTRO_ACCELERATE_FRESNL_HPP
#define ASTRO_ACCELERATE_FRESNL_HPP

namespace astroaccelerate {

  int fresnl(double xxa, double *ssa, double *cca);
  double polevl(double x, void *coef, int N);
  double p1evl(double x, void *coef, int N);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_FRESNL_HPP
