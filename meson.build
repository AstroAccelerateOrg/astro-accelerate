# Meson build script for pacer imager
project('astroaccelerate', 'cpp', 'cuda', version: '1.3.1')

# Basic setup
cc = meson.get_compiler('cpp')
c_std = 'c++14'
buildtype = 'debug'
comp_flags = ['-g', '-std=c++14']
add_project_arguments(comp_flags, language : ['c', 'cpp', 'cuda'])


# Load list of files to compile into shared object
# Note: Meson does not allow * wildcard (deliberately)
# Instead, we use file generated from command: ls src/*.cpp > meson_srclist.txt
fs = import('fs') 
srcs = fs.read('meson_srclist.txt').strip().split('\n')
apps = fs.read('meson_applist.txt').strip().split('\n')

# Dependencies: CUDA
cuda_dep = dependency('cuda', version : '>=11', modules : ['cudart', 'nvrtc', 'cufft', 'curand'])
openmp_dep = dependency('openmp')

astroac_lib = library('astroaccelerate', srcs,              
                     include_directories: ['src', 'include'],
                     dependencies: [cuda_dep, openmp_dep],
                     install: true
                    )

astroac_exe = executable('astro-accelerate', 'src/aa_main.cpp',
                include_directories: ['src', 'include'],
                dependencies: [cuda_dep, openmp_dep],
                link_with: astroac_lib,
                install: true
                )


subdir('examples')
subdir('include')