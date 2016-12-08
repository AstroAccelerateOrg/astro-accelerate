#include "../SpsParameters.h"

namespace ska {
namespace astroaccelerate {
namespace sps{

template<typename DerivedType>
SpsParameters<DerivedType>::SpsParameters()
{

}

template<typename DerivedType>
SpsParameters<DerivedType>::~SpsParameters()
{
}

template<typename DerivedType>
constexpr int		 SpsParameters<DerivedType>::get_accmax(){ return DerivedType::_accmax; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_accstep(){ return DerivedType::_accstep; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_unrolls(){ return DerivedType::_unrolls; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_snumreg(){ return DerivedType::_snumreg; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_sdivint(){ return DerivedType::_sdivint; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_sdivindm(){ return DerivedType::_sdivindm; }
template<typename DerivedType>
constexpr float	SpsParameters<DerivedType>::get_sfdivindm(){ return DerivedType::_sfdivindm; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_card(){ return DerivedType::_card; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_nopsshift(){ return DerivedType::_nopsshift; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_nopsloop(){ return DerivedType::_nopsloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_ndataperloop(){ return DerivedType::_ndataperloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_bindivint(){ return DerivedType::_bindivint; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_bindivinf(){ return DerivedType::_bindivinf; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_ct(){ return DerivedType::_ct; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_cf(){ return DerivedType::_cf; }
template<typename DerivedType>
constexpr float	SpsParameters<DerivedType>::get_nops(){ return DerivedType::_nops; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_statst(){ return DerivedType::_statst; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_statsloop(){ return DerivedType::_statsloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_warp(){return DerivedType::_warp;}
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_half_warp(){ return DerivedType::_half_warp; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_msd_elem_per_thread(){ return DerivedType::_msd_elem_per_thread; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_msd_warps_per_block(){ return DerivedType::_msd_warps_per_block; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_thr_elem_per_thread(){ return DerivedType::_thr_elem_per_thread; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_thr_warps_per_block(){ return DerivedType::_thr_warps_per_block; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_nthreads(){ return DerivedType::_pd_nthreads; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_nwindows(){ return DerivedType::_pd_nwindows; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_maxtaps(){ return DerivedType::_pd_maxtaps; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_smem_size(){ return DerivedType::_pd_smem_size; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_fir_active_warps(){ return DerivedType::_pd_fir_active_warps; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_fir_nwindows(){ return DerivedType::_pd_fir_nwindows; }

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
