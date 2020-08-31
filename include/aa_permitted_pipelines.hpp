#ifndef ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_HPP
#define ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_HPP

#include "aa_pipeline.hpp"
#include "aa_permitted_pipelines.hpp"
#include "aa_permitted_pipelines_0.hpp"
#include "aa_permitted_pipelines_1.hpp"

namespace astroaccelerate {

  /**
   * \class aa_permitted_pipelines aa_permitted_pipelines.hpp "include/aa_permitted_pipelines.hpp"
   * \brief Class that is used to check whether a pipeline is valid and permitted.
   * \details This class can also be used to obtain a valid pipeline.
   * \author Cees Carels.
   * \date 22 October 2018.
   */

  class aa_permitted_pipelines {
  public:
    
    //Example valid pipelines
    static const aa_pipeline::pipeline pipeline0;
    static const aa_pipeline::pipeline pipeline1;
    static const aa_pipeline::pipeline pipeline2;
    static const aa_pipeline::pipeline pipeline3;
    static const aa_pipeline::pipeline pipeline3_0;
    static const aa_pipeline::pipeline pipeline4;
    static const aa_pipeline::pipeline pipeline4_0;
    static const aa_pipeline::pipeline pipeline5;
    static const aa_pipeline::pipeline pipeline5_0;

    /**
     * \brief Pass a pipeline object and the function validates the components it contains.
     * \returns A boolean to indicate whether the pipeline is permitted (true) or not (false).
     */
    static bool is_permitted(const aa_pipeline::pipeline &pipeline) {
		// With generic pipeline in mind we need to rethink this.
      if(pipeline == pipeline0) {
	return true;
      }
      else if(pipeline == pipeline1) {
	return true;
      }
      else if(pipeline == pipeline2) {
	return true;
      }
      else if(pipeline == pipeline3) {
	return true;
      }
      else if(pipeline == pipeline3_0) {
	return true;
      }
      else if(pipeline == pipeline4) {
	return true;
      }
      else if(pipeline == pipeline4_0) {
	return true;
      }
      else if(pipeline == pipeline5) {
	return true;
      }
      else if(pipeline == pipeline5_0) {
	return true;
      }
      else {
	return true;
      }
        
      return false;
    }
    
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_PERMITTED_PIPELINES_HPP
