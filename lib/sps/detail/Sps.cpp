#include "AstroAccelerate/sps/Sps.h"


namespace ska {
namespace astroaccelerate {
namespace sps {

template<typename SpsParameterType>
Sps<SpsParameterType>::Sps()
{
}

template<typename SpsParameterType>
Sps<SpsParameterType>::~Sps()
{
}

template<typename SpsParameterType>
template<typename SpsHandler, typename DmHandler>
Sps<SpsParameterType>::operator()(DeviceId, DataInputType const&, DataOutputType&, DedispersionPlan const&, SpsHandler, DmHandler)
{
}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
