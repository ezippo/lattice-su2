bin_PROGRAMS = conf_check \
               debug_su2 \
               main_su2

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

gauge_conf = ./include/gauge_conf.h \
             ./lib/gauge_conf_def.c \
             ./lib/gauge_conf_meas.c \
             ./lib/gauge_conf_upd.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

su2 = ./include/su2.h \
      ./include/su2_upd.h \
      ./lib/su2.c \
      ./lib/su2_upd.c

everything =  $(dSFMT)  \
              $(endian) \
              $(gauge_conf) \
              $(geometry) \
              $(gparam) \
              ./include/macro.h \
              $(random) \
              $(su2) 

conf_check_SOURCES=./src/conf_check.c \
                   $(dSFMT) \
                   $(endian) \
                   $(gparam) \
                   $(random) \
                   $(su2) 

debug_su2_SOURCES=./src/debug_su2.c \
                  $(dSFMT) \
                  $(endian) \
                  $(gparam) \
                  ./include/macro.h \
                  $(random) \
                  $(su2) 

main_su2_SOURCES=./src/main_su2.c \
                          $(everything)

AM_CFLAGS = -std=c99 -O2 -Wall -Wextra -Werror -pedantic -Wconversion @OPENMP_CFLAGS@

CLEANFILES = log* dati* conf.dat* template_input* conf5D* 
