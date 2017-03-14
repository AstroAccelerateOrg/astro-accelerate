#include "../AstroAccelerateParameters.h"

namespace ska {
namespace astroaccelerate {

template<typename DerivedType>
AstroAccelerateParameters<DerivedType>::AstroAccelerateParameters()
{

}

template<typename DerivedType>
AstroAccelerateParameters<DerivedType>::~AstroAccelerateParameters()
{
}

template<typename DerivedType>
constexpr int		 AstroAccelerateParameters<DerivedType>::get_accmax(){ return DerivedType::_accmax; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_accstep(){ return DerivedType::_accstep; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_unrolls(){ return DerivedType::_unrolls; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_snumreg(){ return DerivedType::_snumreg; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_sdivint(){ return DerivedType::_sdivint; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_sdivindm(){ return DerivedType::_sdivindm; }
template<typename DerivedType>
constexpr float	AstroAccelerateParameters<DerivedType>::get_sfdivindm(){ return DerivedType::_sfdivindm; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_card(){ return DerivedType::_card; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_nopsshift(){ return DerivedType::_nopsshift; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_nopsloop(){ return DerivedType::_nopsloop; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_ndataperloop(){ return DerivedType::_ndataperloop; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_bindivint(){ return DerivedType::_bindivint; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_bindivinf(){ return DerivedType::_bindivinf; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_ct(){ return DerivedType::_ct; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_cf(){ return DerivedType::_cf; }
template<typename DerivedType>
constexpr float	AstroAccelerateParameters<DerivedType>::get_nops(){ return DerivedType::_nops; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_statst(){ return DerivedType::_statst; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_statsloop(){ return DerivedType::_statsloop; }
template<typename DerivedType>
// sps parameters
constexpr int 	AstroAccelerateParameters<DerivedType>::get_warp(){return DerivedType::_warp;}
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_half_warp(){ return DerivedType::_half_warp; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_msd_elem_per_thread(){ return DerivedType::_msd_elem_per_thread; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_msd_warps_per_block(){ return DerivedType::_msd_warps_per_block; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_thr_elem_per_thread(){ return DerivedType::_thr_elem_per_thread; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_thr_warps_per_block(){ return DerivedType::_thr_warps_per_block; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_pd_nthreads(){ return DerivedType::_pd_nthreads; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_pd_nwindows(){ return DerivedType::_pd_nwindows; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_pd_maxtaps(){ return DerivedType::_pd_maxtaps; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_pd_smem_size(){ return DerivedType::_pd_smem_size; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_pd_fir_active_warps(){ return DerivedType::_pd_fir_active_warps; }
template<typename DerivedType>
constexpr int 	AstroAccelerateParameters<DerivedType>::get_fir_nwindows(){ return DerivedType::_pd_fir_nwindows; }
// fdas parameters
template<typename DerivedType>
constexpr float AstroAccelerateParameters<DerivedType>::get_tsamp(){ return DerivedType::_tsamp; }
template<typename DerivedType>
constexpr float AstroAccelerateParameters<DerivedType>::get_slight(){ return DerivedType::_slight; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_radix(){ return DerivedType::_radix; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_nexp(){ return DerivedType::_nexp; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_potwo(){ return DerivedType::_potwo; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_kernlen(){ return DerivedType::_kernlen; }
template<typename DerivedType>
constexpr float AstroAccelerateParameters<DerivedType>::get_accel_step(){ return DerivedType::_accel_step; }
template<typename DerivedType>
constexpr float AstroAccelerateParameters<DerivedType>::get_accel_step_r(){ return DerivedType::_accel_step_r; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_zmax(){ return DerivedType::_zmax; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_nkern(){ return DerivedType::_nkern; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_tbsizex(){ return DerivedType::_tbsizex; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_tbsizey(){ return DerivedType::_tbsizey; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_ptbsizex(){ return DerivedType::_ptbsizex; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_ptbsizey(){ return DerivedType::_ptbsizey; }
template<typename DerivedType>
constexpr int AstroAccelerateParameters<DerivedType>::get_taps(){ return DerivedType::_taps; }

} // namespace astroaccelerate
}
