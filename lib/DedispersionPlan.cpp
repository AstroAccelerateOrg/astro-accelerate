#include "AstroAccelerate/DedispersionPlan.h"


namespace ska {
namespace astroaccelerate {

DedispersionPlan::DedispersionPlan()
{
   	_user_dm_low	= NULL;
	_user_dm_high 	= NULL;
	_user_dm_step 	= NULL;
	_dm_low 		= NULL;		
	_dm_high 		= NULL;
	_dm_step 		= NULL;
	_dmshifts		= NULL;
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


float * DedispersionPlan::get_user_dm_low()
{
	return _user_dm_low ;
}

float * DedispersionPlan::get_user_dm_high()
{
	return _user_dm_high;
}

float * DedispersionPlan::get_user_dm_step()
{
	return _user_dm_step;
}

DedispersionPlan::get_dm_low()
{
	return _dm_low;
}

DedispersionPlan::get_dm_high()
{
	return _dm_high;
}

DedispersionPlan::get_dm_step()
{
	return _dm_step;
}

DedispersionPlan::get_dm_step()
{
	return _dmshifts;
}
	
} // namespace astroaccelerate
} // namespace ska
