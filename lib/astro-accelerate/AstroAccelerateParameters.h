#ifndef ASTROACCELERATE_ASTROACCELERATEPARAMETERS_H
#define ASTROACCELERATE_ASTROACCELERATEPARAMETERS_H

namespace astroaccelerate {

/**
 * @brief		This object carries the parameters needed by the astroaccelerate kernels
 * 				It is not used yet
 * 
 * @details
 * 
 */

template<typename DerivedType>
class AstroAccelerateParameters
{
    public:
  		/**
  		* @brief Default constructor
  		*/
        AstroAccelerateParameters();

        /**
         * @brief Destructor
         */
        ~AstroAccelerateParameters();

        /**
         * @brief getters
         */
		static constexpr int get_accmax();
		static constexpr int get_accstep();
		static constexpr int get_unrolls();
		static constexpr int get_snumreg();
		static constexpr int get_sdivint();
		static constexpr int get_sdivindm();
		static constexpr float get_sfdivindm();
		static constexpr int get_card();
		static constexpr int get_nopsshift();
		static constexpr int get_nopsloop();
		static constexpr int get_ndataperloop();
		static constexpr int get_bindivint();
		static constexpr int get_bindivinf();
		static constexpr int get_ct();
		static constexpr int get_cf();
		static constexpr float get_nops();
		static constexpr int get_statst();
		static constexpr int get_statsloop();
		static constexpr int get_warp();
		static constexpr int get_half_warp();
		static constexpr int get_msd_elem_per_thread();
		static constexpr int get_msd_warps_per_block();
		static constexpr int get_thr_elem_per_thread();
		static constexpr int get_thr_warps_per_block();
		static constexpr int get_pd_nthreads();
		static constexpr int get_pd_nwindows();
		static constexpr int get_pd_maxtaps();
		static constexpr int get_pd_smem_size();
		static constexpr int get_pd_fir_active_warps();
		static constexpr int get_fir_nwindows();
        static constexpr float get_tsamp();
        static constexpr float get_slight();
        static constexpr int get_radix();
        static constexpr int get_nexp();
        static constexpr int get_potwo();
        static constexpr int get_kernlen();
        static constexpr float get_accel_step();
        static constexpr float get_accel_step_r();
        static constexpr int get_zmax();
        static constexpr int get_nkern();
        static constexpr int get_tbsizex();
        static constexpr int get_tbsizey();
        static constexpr int get_ptbsizex();
        static constexpr int get_ptbsizey();
        static constexpr int get_taps();

	private:
		/**
		 * @brief dedispersion parameters
		 */
		static constexpr int _accmax = 350;
		static constexpr int _accstep = 11;
		static constexpr int _unrolls = 16;
		static constexpr int _snumreg = 12;
		static constexpr int _sdivint = 12;
		static constexpr int _sdivindm = 40;
		static constexpr float _sfdivindm = 40.0f;
		static constexpr int _nopsshift = 5;
		static constexpr int _nopsloop = 3;
		static constexpr int _ndataperloop = 1;
		static constexpr int _bindivint = 8;
		static constexpr int _bindivinf = 32;
		static constexpr int _ct = 32;
		static constexpr int _cf = 8;
		static constexpr float _nops = 4.0;
		static constexpr int _statst = 128;
		static constexpr int _statsloop = 8;
		/**
         * @brief single pulse search parameters
         */
        static constexpr int _warp = 32;
        static constexpr int _half_warp = 16;
        static constexpr int _msd_elem_per_thread = 8;
        static constexpr int _msd_warps_per_block = 16;
        static constexpr int _thr_elem_per_thread = 4;
        static constexpr int _thr_warps_per_block = 4;
        static constexpr int _pd_nthreads = 512;
        static constexpr int _pd_nwindows = 2;
        static constexpr int _pd_maxtaps = 32;
        static constexpr int _pd_smem_size = 1280;
        static constexpr int _pd_fir_active_warps = 2;
        static constexpr int _pd_fir_nwindows = 2;
        /**
         * @brief fdas parameters
         */
        static constexpr float _tsamp = 0.000064;
        static constexpr float _slight = 299792458.0;
        static constexpr int _radix = 1;
        static constexpr int _nexp = 10;
        static constexpr int _potwo = 1<<_nexp;
        static constexpr int _kernlen = _radix*_potwo;
        static constexpr float _accel_step = (float)(2.0);
        static constexpr float _accel_step_r = (float)(1/_accel_step);
        static constexpr int _zmax = 96;
        static constexpr int _nkern = _zmax + 1;
        static constexpr int _tbsizex = 32;
        static constexpr int _tbsizey = 1;
        static constexpr int _ptbsizex = 64;
        static constexpr int _ptbsizey = 1;
        static constexpr int _taps = 8;
};

} // namespace astroaccelerate


#include "detail/AstroAccelerateParameters.cpp"

#endif // ASTROACCELERATE_ASTROACCELERATEPARAMETERS_H
