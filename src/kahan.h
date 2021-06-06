/***
 Kahan summation algorithm.
 */
inline void kahan(double &sum, double &comp, double const new_val) noexcept {
double const y = new_val - comp,
             t = sum + y;
  comp = (t - sum) - y;
  sum = t;
}
