#include "../DmTime.h"

namespace astroaccelerate {

template<typename ValueType>
DmTime<ValueType>::DmTime(DedispersionStrategy const& dedispersion_strategy)
{
    _nsamples.clear();
	std::size_t number_of_dm_ranges = dedispersion_strategy.get_range();
    _ndms.reserve(number_of_dm_ranges);
    _data = (ValueType ***) malloc(number_of_dm_ranges*sizeof(ValueType **));
    // ensure _range never gets above our actual malloced range as this will be
    for(_range = 0; _range < number_of_dm_ranges; ++_range )
    {
        std::size_t dm_block_size = (dedispersion_strategy.get_ndms())[_range];
    	int total_samps = 0;
    	for (int k = 0; k < dedispersion_strategy.get_num_tchunks(); ++k)
    		total_samps += dedispersion_strategy.get_t_processed()[_range][k];
        _nsamples.push_back(total_samps);
        _data[_range] = (ValueType **) malloc(dm_block_size*sizeof(ValueType *));
        if (_data[_range] == NULL)
        {
        	fprintf(stderr, "\nDmTime memory allocation failed\n");
        	exit(0);
        }
        _ndms[_range]=0;
        for(int j = 0; j < dm_block_size; ++j)
        {
            _data[_range][j] = (ValueType *) malloc((total_samps)*sizeof(ValueType));
            if (_data[_range][j] == NULL)
            {
            	fprintf(stderr, "\nDmTime memory allocation failed\n");
            	exit(0);
            }
            ++_ndms[_range];
        }
        _ouput_size += total_samps * dm_block_size * sizeof(float);
    }
    assert(_range == number_of_dm_ranges);
}

template<typename ValueType>
DmTime<ValueType>::~DmTime()
{
    for(int i = 0; i < _range; ++i)
    {
        for(int j = 0; j < _ndms[i]; ++j)
        {
            free(_data[i][j]);
        }
        free(_data[i]);
    }
    free(_data);
}

template<typename ValueType>
std::size_t DmTime<ValueType>::number_of_dm_ranges() const
{
    return _range;
}

template<typename ValueType>
std::size_t DmTime<ValueType>::output_size() const
{
    return _ouput_size;
}

} // namespace astroaccelerate
