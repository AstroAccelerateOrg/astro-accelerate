#ifndef SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
#define SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H


namespace ska {
namespace astroaccelerate {
namespace sps {

/**
 * @brief
 * 
 * @details
 * 
 */

class DedispersionPlan
{
    public:
        DedispersionPlan();
        ~DedispersionPlan();

        /**
         * @brief return the minimum number of samples required for the
         *        algorithm to operate
         */
        std::size_t minimum_number_of_samples() const;

    private:
};


} // namespace sps
} // namespace astroaccelerate
} // namespace ska

#endif // SKA_ASTROACCELERATE_SPS_DEDISPERSIONPLAN_H
