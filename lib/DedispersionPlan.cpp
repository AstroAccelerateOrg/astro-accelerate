#include "AstroAccelerate/DedispersionPlan.h"


namespace ska {
namespace astroaccelerate {
namespace sps {

	DedispersionPlan::DedispersionPlan()
	{
		_in_bin					= NULL;
		_out_bin				= NULL;
		_maxshift				= 0;
	  _user_dm_low		= NULL;
		_user_dm_high 	= NULL;
		_user_dm_step 	= NULL;
		_dm_low 				= NULL;
		_dm_high 				= NULL;
		_dm_step 				= NULL;
		_dmshifts				= NULL;
	}

	DedispersionPlan::~DedispersionPlan()
	{
	}

	void DedispersionPlan::set_user_dm_low(float * user_dm_low)
	{
		_user_dm_low = user_dm_low;
	}

	void DedispersionPlan::set_user_dm_high(float * user_dm_high)
	{
		_user_dm_high = user_dm_high;
	}

	void DedispersionPlan::set_user_dm_step(float * user_dm_step)
	{
		_user_dm_step = user_dm_step;
	}

	void DedispersionPlan::set_dm_low(float * dm_low)
	{
		_dm_low = dm_low;
	}

	void DedispersionPlan::set_dm_high(float * dm_high)
	{
		_dm_high = dm_high;
	}

	void DedispersionPlan::set_dm_step(float * dm_step)
	{
		_dm_step = dm_step;
	}

	void DedispersionPlan::set_dmshifts(float * dmshifts)
	{
		_dmshifts = dmshifts;
	}


	int* DedispersionPlan::get_in_bin() const
	{
		return _in_bin;
	}

	int* DedispersionPlan::get_out_bin() const
	{
		return _out_bin;
	}

	int DedispersionPlan::get_maxshift() const
	{
		return _maxshift ;
	}

	float* DedispersionPlan::get_user_dm_low() const
	{
		return _user_dm_low ;
	}

	float* DedispersionPlan::get_user_dm_high() const
	{
		return _user_dm_high;
	}

	float* DedispersionPlan::get_user_dm_step() const
	{
		return _user_dm_step;
	}

	float* DedispersionPlan::get_dm_low() const
	{
		return _dm_low;
	}

	float* DedispersionPlan::get_dm_high() const
	{
		return _dm_high;
	}

	float* DedispersionPlan::get_dm_step() const
	{
		return _dm_step;
	}

	float* DedispersionPlan::get_dmshifts() const
	{
		return _dmshifts;
	}

} // namespace sps
} // namespace astroaccelerate
} // namespace ska
