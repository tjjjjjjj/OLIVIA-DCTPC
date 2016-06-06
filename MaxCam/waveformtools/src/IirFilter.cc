#include "IirFilter.hh"

//Force compilation of double, float, and long double versions
namespace waveform{
  template class IirFilter<double>;
  template class IirFilter<float>;
  template class IirFilter<long double>;
}
