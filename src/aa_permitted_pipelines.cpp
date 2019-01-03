#include "aa_permitted_pipelines.hpp"

namespace astroaccelerate {
  /** \brief Static implementation of permitted pipelines. */  
  const aa_compute::pipeline aa_permitted_pipelines::pipeline0   = {aa_compute::modules::empty};
  const aa_compute::pipeline aa_permitted_pipelines::pipeline1   = {aa_compute::modules::dedispersion};
  const aa_compute::pipeline aa_permitted_pipelines::pipeline2   = {aa_compute::modules::dedispersion, aa_compute::modules::analysis};
  const aa_compute::pipeline aa_permitted_pipelines::pipeline3   = {aa_compute::modules::dedispersion, aa_compute::modules::analysis, aa_compute::modules::periodicity};
  const	aa_compute::pipeline aa_permitted_pipelines::pipeline4   = {aa_compute::modules::dedispersion, aa_compute::modules::analysis, aa_compute::modules::fdas};
  const	aa_compute::pipeline aa_permitted_pipelines::pipeline5   = {aa_compute::modules::dedispersion, aa_compute::modules::analysis, aa_compute::modules::periodicity, aa_compute::modules::fdas};

} //namespace astroaccelerate
