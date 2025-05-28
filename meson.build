project('dpmjet', 'fortran', version: '19.3.7')

add_project_arguments('-O3', '-Wno-uninitialized', '-fPIC', '-std=legacy',
                      '-fno-second-underscore', '-cpp', language: 'fortran')

mold = find_program('mold', required: false)
if mold.found()
  add_project_link_arguments('-fuse-ld=mold', language: 'fortran')
endif

inc = include_directories('include/pythia',
                          'include/phojet',
                          'include/dpmjet',
                          'include/flinclude')

phojet_src = []
foreach f : run_command('sh', '-c', 'ls src/phojet/*.f').stdout().split()
  phojet_src += files(f)
endforeach

dpmjet_src = []
foreach f : run_command('sh', '-c', 'ls src/dpmjet/*.f').stdout().split()
  dpmjet_src += files(f)
endforeach

pythia_src = []
foreach f : run_command('sh', '-c', 'ls src/pythia/*.f').stdout().split()
  pythia_src += files(f)
endforeach

common_src = []
foreach f : run_command('sh', '-c', 'ls common/*.f').stdout().split()
  common_src += files(f)
endforeach

libdpmjet = static_library('DPMJET',
                           phojet_src + dpmjet_src + pythia_src + common_src,
                           include_directories: inc)

foreach f : run_command('sh', '-c', 'ls src/exe/*.f').stdout().split()
  name = f.split('/')[-1].split('.')[0]
  executable(name, files(f),
             include_directories: inc,
             link_with: libdpmjet,
             install: true)
endforeach

install_data('examples/dpmjet/ppLHC.inp', install_dir: get_option('datadir'))
