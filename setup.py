import platform, distutils.core, distutils.extension, Cython.Build


## A hack?  Yes -- Anaconda compiler update broke my linking.
##    and the strict prototypes thing is just an obnoxious holdover from c.
import distutils.sysconfig, re
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    # print(key, cfg_vars[key])
    if type(value) == str:
        value = value.replace("-Wstrict-prototypes", "")
        value = re.sub(r"-B.*compiler_compat", "", value)
        cfg_vars[key] = value

import sys

## Macs require this extra build option.
ouff_mac = []
if sys.platform == "darwin":
  ouff_mac = ['-mmacosx-version-min=10.9']

EXTENSION = distutils.extension.Extension(
    name = 'pyc4', language = 'c++',
    sources = ['pyc4.pyx'],
    extra_compile_args = ['-Wno-unused-function', 
                          '-std=c++11', '-Wall'] + ouff_mac,
    undef_macros       = ["NDEBUG"],
    extra_link_args    = ouff_mac,
    # include_path       = ["/usr/local/include/"],
    # libraries          = ["armadillo"]
    )

EXT_MODULES=Cython.Build.cythonize([EXTENSION],
                                   # include_path = ["/usr/local/include/"],
                                   language='c++')#, gdb_debug=True)

distutils.core.setup(
    name = 'simple',
    ext_modules=EXT_MODULES,
    )

