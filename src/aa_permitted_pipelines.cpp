#include "aa_permitted_pipelines.hpp"

namespace astroaccelerate {
  /** \brief Static implementation of permitted pipelines. */  
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline0   = {aa_pipeline::component::empty};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline1   = {aa_pipeline::component::dedispersion};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline2   = {aa_pipeline::component::dedispersion, aa_pipeline::component::analysis};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline3   = {aa_pipeline::component::dedispersion, aa_pipeline::component::analysis, aa_pipeline::component::periodicity};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline3_0 = {aa_pipeline::component::dedispersion, aa_pipeline::component::periodicity};
  const	aa_pipeline::pipeline aa_permitted_pipelines::pipeline4   = {aa_pipeline::component::dedispersion, aa_pipeline::component::analysis, aa_pipeline::component::fdas};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline4_0 = {aa_pipeline::component::dedispersion, aa_pipeline::component::fdas};
  const	aa_pipeline::pipeline aa_permitted_pipelines::pipeline5   = {aa_pipeline::component::dedispersion, aa_pipeline::component::analysis, aa_pipeline::component::periodicity, aa_pipeline::component::fdas};
  const aa_pipeline::pipeline aa_permitted_pipelines::pipeline5_0 = {aa_pipeline::component::dedispersion, aa_pipeline::component::periodicity, aa_pipeline::component::fdas};

} //namespace astroaccelerate
