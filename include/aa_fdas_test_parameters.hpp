#ifndef ASTRO_ACCELERATE_AA_FDAS_TEST_PARAMETERS_HPP
#define ASTRO_ACCELERATE_AA_FDAS_TEST_PARAMETERS_HPP
namespace astroaccelerate {

  // basic test with top-hat filters
  // FDAS code uses series of filters. This parameter gives increment for the top-hat filter in this test. Effective number of filters is ZMAX/2.
#define FDAS_TEST_FILTER_INCREMENT 2
  // length of one tooth in sawtooth wave
#define FDAS_TEST_TOOTH_LENGTH 4096

  // Accelerated sin wave
#define FDAS_TEST_FREQUENCY 105.0
  // acceleration
#define FDAS_TEST_ZVALUE 6
#define FDAS_TEST_HAMONICS 4
#define FDAS_TEST_DUTY_CYCLE 1.0
#define FDAS_TEST_SIGNAL_AMPLITUDE 1.0

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_FDAS_TEST_PARAMETERS_HPP
