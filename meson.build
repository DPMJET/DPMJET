project('dpmjet', 'fortran', version: '19.3.7', default_options: ['optimization=3', 'fortran_std=legacy'])

add_project_arguments('-Wno-uninitialized', '-fPIC',
                      '-fno-second-underscore', '-Wno-unused', '-Wno-conversion',
                      '-Wno-unused-dummy-argument', '-Wno-character-truncation', 
                      '-Wno-integer-division','-Wno-surprising', '-cpp', language: 'fortran')

# mold = find_program('mold', required: false)
# if mold.found()
#   add_project_link_arguments('-fuse-ld=mold', language: 'fortran')
# endif

inc = include_directories('include/pythia',
                          'include/phojet',
                          'include/dpmjet',
                          'include/flinclude')

phojet_src = []
foreach f : run_command('sh', '-c', 'ls src/phojet/*.f', check:true).stdout().split()
  phojet_src += files(f)
endforeach

dpmjet_src = []
foreach f : run_command('sh', '-c', 'ls src/dpmjet/*.f', check:true).stdout().split()
  dpmjet_src += files(f)
endforeach

pythia_src = []
foreach f : run_command('sh', '-c', 'ls src/pythia/*.f', check:true).stdout().split()
  pythia_src += files(f)
endforeach

common_src = []
foreach f : run_command('sh', '-c', 'ls common/*.f', check:true).stdout().split()
  common_src += files(f)
endforeach

libdpmjet = static_library('DPMJET',
                           phojet_src + dpmjet_src + pythia_src + common_src,
                           include_directories: inc)

foreach f : run_command('sh', '-c', 'ls src/exe/*.f', check:true).stdout().split()
  name = f.split('/')[-1].split('.')[0]
  executable(name, files(f),
             include_directories: inc,
             link_with: libdpmjet,
             install: true,
             install_dir: 'bin')
endforeach

# # Custom clean target to remove both build artifacts and installed executables
# run_target('distclean',
#   command: [
#     'sh', '-c', 
#     'rm -rf build/* bin/DPMJET* bin/PHOJET* bin/pho_* bin/photo_* || true'
#   ]
# )
