examples = [
    'dedispersion',
    'dedispersion_and_analysis',
    'fake_signal_periodic',
    'fake_signal_single',
    'filterbank_dedispersion',
    'periodicity'
]

foreach app : examples
  exe = executable(app, './src'/app+'.cpp',
                include_directories: ['src', '../include'],
                dependencies: [cuda_dep, openmp_dep],
                link_with: astroac_lib,
                install: false
                )
endforeach
