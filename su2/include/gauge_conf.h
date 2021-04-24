#ifndef GAUGE_CONF_H
#define GAUGE_CONF_H

#include"macro.h"

#include<complex.h>
#include<openssl/md5.h>
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"
#include"su2.h"


typedef struct Gauge_Conf {
  long update_index;

  Su2 **lattice;       // [volume] [STDIM]
  } Gauge_Conf;


// in gauge_conf_def.c
void init_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void read_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void free_gauge_conf(Gauge_Conf *GC,
                     GParam const * const param);
void write_conf_on_file_with_name(Gauge_Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Gauge_Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Gauge_Conf const * const GC,
                             GParam const * const param);
void init_gauge_conf_from_gauge_conf(Gauge_Conf *GC,
                                     Gauge_Conf const * const GC2,
                                     GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Gauge_Conf const * const GC,
                         GParam const * const param);

// in gauge_conf_meas.c
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j);
double wilson_loopp(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i,
                    int j,
                    int const size_i,
                    int const size_j);
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt);
void plaquette_adj(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *plaqs,
              double *plaqt);
void wilson_loop(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 int const size_i,
                 int const size_j,
                 double *wloop_s,
                 double *wloop_t);
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly);
void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep);

// in gauge_conf_upd.c
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const gparam,
                        long r,
                        int i,
                        Su2 *M);
void calcstaples_wilson_nosum(Gauge_Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const gparam,
                              long r,
                              int i,
                              Su2 *M);
void heatbath(Gauge_Conf * GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i);
void overrelaxation(Gauge_Conf * GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i);
int metropolis(Gauge_Conf *GC,
               Geometry const * const geo,
               GParam const * const param,
               long r,
               int i);
int metropolis_fund_plus_adj(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i);
void update(Gauge_Conf *GC,
            Geometry const * const geo,
            GParam const * const param);
double update_metropolis(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param);

#endif
