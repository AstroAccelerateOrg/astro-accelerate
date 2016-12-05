#include "../SpsParameters.h"

namespace ska {
namespace astroaccelerate {

template<typename DerivedType>
SpsParameters<DerivedType>::SpsParameters()
{

}

template<typename DerivedType>
SpsParameters<DerivedType>::~SpsParameters()
{
}

template<typename DerivedType>
constexpr int		 SpsParameters<DerivedType>::get_accmax(){ return _accmax; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_accstep(){ return _accstep; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_unrolls(){ return _unrolls; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_snumreg(){ return _snumreg; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_sdivint(){ return _sdivint; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_sdivindm(){ return _sdivindm; }
template<typename DerivedType>
constexpr float	SpsParameters<DerivedType>::get_sfdivindm(){ return _sfdivindm; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_card(){ return _card; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_nopsshift(){ return _nopsshift; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_nopsloop(){ return _nopsloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_ndataperloop(){ return _ndataperloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_bindivint(){ return _bindivint; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_bindivinf(){ return _bindivinf; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_ct(){ return _ct; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_cf(){ return _cf; }
template<typename DerivedType>
constexpr float	SpsParameters<DerivedType>::get_nops(){ return _nops; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_statst(){ return _statst; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_statsloop(){ return _statsloop; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_warp(){return DerivedType::_warp;}
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_half_warp(){ return _half_warp; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_msd_elem_per_thread(){ return _msd_elem_per_thread; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_msd_warps_per_block(){ return _msd_warps_per_block; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_thr_elem_per_thread(){ return _thr_elem_per_thread; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_thr_warps_per_block(){ return _thr_warps_per_block; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_nthreads(){ return _pd_nthreads; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_nwindows(){ return _pd_nwindows; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_maxtaps(){ return _pd_maxtaps; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_smem_size(){ return _pd_smem_size; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_pd_fir_active_warps(){ return _pd_fir_active_warps; }
template<typename DerivedType>
constexpr int 	SpsParameters<DerivedType>::get_fir_nwindows(){ return _pd_fir_nwindows; }

} // namespace astroaccelerate
} // namespace ska
