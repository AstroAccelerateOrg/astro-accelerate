#ifndef SKA_ASTROACCELERATE_SPS_SPSPARAMETERS_H
#define SKA_ASTROACCELERATE_SPS_SPSPARAMETERS_H


namespace ska {
namespace astroaccelerate {
namespace sps{
/**
 * @brief		This object carries the parameters needed by the SPS kernels
 * 
 * @details
 * 
 */

template<typename DerivedType>
class SpsParameters
{
    public:
  		/**
  		* @brief Default constructor
  		*/
        SpsParameters();

        /**
         * @brief Destructor
         */
        ~SpsParameters();

        /**
         * @brief No setters needed since it's constant values.
         */

        /**
         * @brief getters
         */
        static constexpr int		get_accmax();
        static constexpr int 		get_accstep();
        static constexpr int 		get_unrolls();
        static constexpr int 		get_snumreg();
        static constexpr int 		get_sdivint();
        static constexpr int 		get_sdivindm();
        static constexpr float		get_sfdivindm();
        static constexpr int 		get_card();
        static constexpr int 		get_nopsshift();
        static constexpr int 		get_nopsloop();
        static constexpr int 		get_ndataperloop();
        static constexpr int 		get_bindivint();
        static constexpr int 		get_bindivinf();
        static constexpr int 		get_ct();
        static constexpr int 		get_cf();
        static constexpr float		get_nops();
        static constexpr int 		get_statst();
        static constexpr int 		get_statsloop();
        static constexpr int 		get_warp();
        static constexpr int 		get_half_warp();
        static constexpr int 		get_msd_elem_per_thread();
        static constexpr int 		get_msd_warps_per_block();
        static constexpr int 		get_thr_elem_per_thread();
        static constexpr int 		get_thr_warps_per_block();
        static constexpr int 		get_pd_nthreads();
        static constexpr int 		get_pd_nwindows();
        static constexpr int 		get_pd_maxtaps();
        static constexpr int 		get_pd_smem_size();
        static constexpr int 		get_pd_fir_active_warps();
        static constexpr int 		get_fir_nwindows();
        
    private:
        /**
         * @brief some parameters
         */
        static constexpr int 	_accmax 	= 350;
        static constexpr int 	_accstep 	= 11;
        static constexpr int 	_unrolls	= 16;
        static constexpr int 	_snumreg	= 10;
        static constexpr int 	_sdivint	= 8;
        static constexpr int 	_sdivindm	= 40;
        static constexpr float _sfdivindm	= 40.0f;
        static constexpr int 	_card		= 0;
        static constexpr int 	_nopsshift	= 5;
        static constexpr int 	_nopsloop		= 3;
        static constexpr int 	_ndataperloop	= 1;
        static constexpr int 	_bindivint		= 6;
        static constexpr int 	_bindivinf		= 32;
        static constexpr int 	_ct				= 256;
        static constexpr int 	_cf				= 2;
        static constexpr float	_nops			= 4.0;
        static constexpr int 	_statst			= 128;
        static constexpr int 	_statsloop		= 8;
        /**
         * @brief some other parameters // added by KA
         */
        static constexpr int _warp					= 32;
        static constexpr int _half_warp				= 16;
        static constexpr int _msd_elem_per_thread	= 8;
        static constexpr int _msd_warps_per_block	= 16;
        static constexpr int _thr_elem_per_thread	= 4;
        static constexpr int _thr_warps_per_block	= 4;
        static constexpr int _pd_nthreads			= 512;
        static constexpr int _pd_nwindows			= 2;
        static constexpr int _pd_maxtaps			= 16;
        static constexpr int _pd_smem_size			= 1280;
        static constexpr int _pd_fir_active_warps	= 2;
        static constexpr int _pd_fir_nwindows		= 2;

};

} // namespace sps
} // namespace astroaccelerate
} // namespace ska


#include "detail/SpsParameters.cpp"

#endif // SKA_ASTROACCELERATE_TOOLS_SPSPARAMETERS_H
