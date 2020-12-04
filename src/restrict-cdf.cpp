#include "restrict-cdf.h"

namespace restrictcdf {

cache_mem<double> likelihood::dmen;
cache_mem<double> deriv::dmem;
cache_mem<double> imputation::dmem;

} // namespace restrictcdf
