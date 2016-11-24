#include "AstroAccelerate/SpsParameters.h"


namespace ska {
namespace astroaccelerate {

SpsParameters::SpsParameters()
{

}

SpsParameters::~SpsParameters()
{
}

constexpr int		SpsParameters::get_accmax(){ return _accmax; }
constexpr int 	SpsParameters::get_accstep(){ return _accstep; }
constexpr int 	SpsParameters::get_unrolls(){ return _unrolls; }
constexpr int 	SpsParameters::get_snumreg(){ return _snumreg; }
constexpr int 	SpsParameters::get_sdivint(){ return _sdivint; }
constexpr int 	SpsParameters::get_sdivindm(){ return _sdivindm; }
constexpr float	SpsParameters::get_sfdivindm(){ return _sfdivindm; }
constexpr int 	SpsParameters::get_card(){ return _card; }
constexpr int 	SpsParameters::get_nopsshift(){ return _nopsshift; }
constexpr int 	SpsParameters::get_nopsloop(){ return _nopsloop; }
constexpr int 	SpsParameters::get_ndataperloop(){ return _ndataperloop; }
constexpr int 	SpsParameters::get_bindivint(){ return _bindivint; }
constexpr int 	SpsParameters::get_bindivinf(){ return _bindivinf; }
constexpr int 	SpsParameters::get_ct(){ return _ct; }
constexpr int 	SpsParameters::get_cf(){ return _cf; }
constexpr float	SpsParameters::get_nops(){ return _nops; }
constexpr int 	SpsParameters::get_statst(){ return _statst; }
constexpr int 	SpsParameters::get_statsloop(){ return _statsloop; }
constexpr int 	SpsParameters::get_warp(){ return _warp; }
constexpr int 	SpsParameters::get_half_warp(){ return _half_warp; }
constexpr int 	SpsParameters::get_msd_elem_per_thread(){ return _msd_elem_per_thread; }
constexpr int 	SpsParameters::get_msd_warps_per_block(){ return _msd_warps_per_block; }
constexpr int 	SpsParameters::get_thr_elem_per_thread(){ return _thr_elem_per_thread; }
constexpr int 	SpsParameters::get_thr_warps_per_block(){ return _thr_warps_per_block; }
constexpr int 	SpsParameters::get_pd_nthreads(){ return _pd_nthreads; }
constexpr int 	SpsParameters::get_pd_nwindows(){ return _pd_nwindows; }
constexpr int 	SpsParameters::get_pd_maxtaps(){ return _pd_maxtaps; }
constexpr int 	SpsParameters::get_pd_smem_size(){ return _pd_smem_size; }
constexpr int 	SpsParameters::get_pd_fir_active_warps(){ return _pd_fir_active_warps; }
constexpr int 	SpsParameters::get_fir_nwindows(){ return _pd_fir_nwindows; }

} // namespace astroaccelerate
} // namespace ska
