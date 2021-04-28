#ifndef GAUGE_CONF_MEAS_C
#define GAUGE_CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/gauge_conf.h"


// computation of the plaquette (1/NCOLOR the trace of) in position r and positive directions i,j
double plaquettep(Gauge_Conf const * const GC,
                  Geometry const * const geo,
                  GParam const * const param,
                  long r,
                  int i,
                  int j)
   {
   Su2 matrix;

   #ifdef DEBUG
   if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #else
   (void) param; // just to avoid warning at compile time
   #endif

//
//       ^ i
//       |   (2)
//       +---<---+
//       |       |
//   (3) V       ^ (1)
//       |       |
//       +--->---+---> j
//       r   (4)
//

   equal(&matrix, &(GC->lattice[nnp(geo, r, j)][i]));
   times_equal_dag(&matrix, &(GC->lattice[nnp(geo, r, i)][j]));
   times_equal_dag(&matrix, &(GC->lattice[r][i]));
   times_equal(&matrix, &(GC->lattice[r][j]));

   return retr(&matrix);
   }


// computation of the rectangulare wilson loop (1/NCOLOR the trace of) size_i*size_j in position r and positive directions i,j
double wilson_loopp(Gauge_Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i,
                    int j,
                    int const size_i,
                    int const size_j)
{
  Su2 matrix_loop, matrix_43;
  long r_new12, r_new43;
  int count;

  #ifdef DEBUG
    if(r >= param->d_volume)
     {
     fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
    if(j >= STDIM || i >= STDIM)
     {
     fprintf(stderr, "i or j too large: (i=%d || j=%d) >= %d (%s, %d)\n", i, j, STDIM, __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
  #else
    (void) param; // just to avoid warning at compile time
  #endif

//
//           ^ i
//           |     (3)
// r+size_i  +------<------+
//           |             |
//       (4) V             ^ (2)
//           |             |
//           +------>------+---> j
//         r      (1)     r+size_j
//

  // degenerate rectangular case: loop = Id
  if(size_i==0 || size_j==0)  return 1.0;

  // (1)
  equal(&matrix_loop, &(GC->lattice[r][j]));              // M=U_j(r)
  r_new12 = nnp(geo,r,j);                                   // r -> r+j
  for(count=1; count<size_j; count++)
  {
    times_equal(&matrix_loop, &(GC->lattice[r_new12][j]));   // M=U_j(r)U_j(r+j)...U_j(r+count*j)
    r_new12 = nnp(geo,r_new12,j);
  }

  // (2)
  for(count=0; count<size_i; count++)
  {
    times_equal(&matrix_loop, &(GC->lattice[r_new12][i]));   // M=(1)*U_i(r+size_j)...U_i(r+size_j+count*i)
    r_new12 = nnp(geo,r_new12,i);
  }

  // (4)
  equal(&matrix_43, &(GC->lattice[r][i]));              // M43= U_i(r)
  r_new43 = nnp(geo,r,i);
  for(count=1; count<size_i; count++)
  {
    times_equal(&matrix_43, &(GC->lattice[r_new43][i]));   // M43=U_i(r)U_i(r+i)...U_i(r+count*i)
    r_new43 = nnp(geo,r_new43,i);
  }

  // (3)
  for(count=0; count<size_j; count++)
  {
    times_equal(&matrix_43, &(GC->lattice[r_new43][j]));   // M43=(4)*U_j(r+size_i)...U_j(r+size_i+count*j)
    r_new43 = nnp(geo,r_new43,j);
  }

  #ifdef DEBUG
    if(r_new12!=r_new43)
    {
      fprintf(stderr, "Error in Wilson loop computation: loop seems not to be rectangular...\n r12=%ld , r43=%ld \n", r_new12, r_new43);
      exit(EXIT_FAILURE);
    }
  #endif

  times_equal_dag(&matrix_loop, &matrix_43);      // M=(1)*(2)*( (4)(3) )^dag
  return retr(&matrix_loop);
}


// compute the mean plaquettes (spatial, temporal)
void plaquette(Gauge_Conf const * const GC,
               Geometry const * const geo,
               GParam const * const param,
               double *plaqs,
               double *plaqt)
   {
   long r;
   double ps, pt;

   ps=0.0;
   pt=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      int i, j;
      i=0;
      for(j=1; j<STDIM; j++)
         {
         pt+=plaquettep(GC, geo, param, r, i, j);
         }

      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            ps+=plaquettep(GC, geo, param, r, i, j);
            }
         }
      }

   if(STDIM>2)
     {
     ps*=param->d_inv_vol;
     ps/=((double) (STDIM-1)*(STDIM-2)/2);
     }
   else
     {
     ps=0.0;
     }

   pt*=param->d_inv_vol;
   pt/=((double) STDIM-1);

   *plaqs=ps;
   *plaqt=pt;
   }


// compute the mean adjoint plaquettes (spatial, temporal)
void plaquette_adj(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *plaqs,
              double *plaqt)
  {
  long r;
  double ps, pt;

  ps=0.0;
  pt=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : pt) reduction(+ : ps)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     int i, j;
     i=0;
     for(j=1; j<STDIM; j++)
        {
        double tmp;
        tmp=plaquettep(GC, geo, param, r, i, j);
        pt+=tmp*tmp;
        }

     for(i=1; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           double tmp;
           tmp=plaquettep(GC, geo, param, r, i, j);
           ps+=tmp*tmp;
           }
        }
     }

  if(STDIM>2)
    {
    ps*=param->d_inv_vol;
    ps/=((double) (STDIM-1)*(STDIM-2)/2);
    }
  else
    {
    ps=0.0;
    }

  pt*=param->d_inv_vol;
  pt/=((double) STDIM-1);

  *plaqs=ps*4.0/3.0;
  *plaqt=pt*4.0/3.0;
  }


// compute the mean wilson loop size_i*size_j (spatial, temporal)
void wilson_loop(Gauge_Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param,
                 int const size_i,
                 int const size_j,
                 double *wloop_s,
                 double *wloop_t)
{
  long r;
  double ws=0.0, wt=0.0;

  #ifdef OPENMP_MODE
  #pragma omp parallel for num_threads(NTHREADS) private(r) reduction(+ : wt) reduction(+ : ws)
  #endif
  for(r=0; r<(param->d_volume); r++)
     {
     int i, j;
     i=0;
     for(j=1; j<STDIM; j++)
        {
        wt+=wilson_loopp(GC, geo, param, r, i, j, size_i, size_j);
        }

     for(i=1; i<STDIM; i++)
        {
        for(j=i+1; j<STDIM; j++)
           {
           ws+=wilson_loopp(GC, geo, param, r, i, j, size_i, size_j);
           }
        }
     }

  if(STDIM>2)
    {
    ws*=param->d_inv_vol;
    ws/=((double) (STDIM-1)*(STDIM-2)/2);
    }
  else
    {
    ws=0.0;
    }

  wt*=param->d_inv_vol;
  wt/=((double) STDIM-1);

  *wloop_s=ws;
  *wloop_t=wt;
}


// compute the mean Polyakov loop (the trace of)
void polyakov(Gauge_Conf const * const GC,
              Geometry const * const geo,
              GParam const * const param,
              double *repoly,
              double *impoly)
   {
   long rsp;
   double rep, imp;

   rep=0.0;
   imp=0.0;

   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(rsp) reduction(+ : rep) reduction(+ : imp)
   #endif
   for(rsp=0; rsp<param->d_space_vol; rsp++)
      {
      long r;
      int i;
      Su2 matrix;

      r=sisp_and_t_to_si(geo, rsp, 0);

      one(&matrix);
      for(i=0; i<param->d_size[0]; i++)
         {
         times_equal(&matrix, &(GC->lattice[r][0]));
         r=nnp(geo, r, 0);
         }

      rep+=retr(&matrix);
      imp+=imtr(&matrix);
      }

   *repoly=rep*param->d_inv_space_vol;
   *impoly=imp*param->d_inv_space_vol;
   }


void perform_measures_localobs(Gauge_Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param,
                               FILE *datafilep)
   {
   double plaqs, plaqt, wloops, wloopt, wloops_11, wloopt_11, wloops_10, wloopt_10, wloops_44, wloopt_44, wloops_42, wloopt_42, polyre, polyim;
   int size1, size2;
   size1 = param->d_loop_size[0];
   size2 = param->d_loop_size[1];
   plaquette(GC, geo, param, &plaqs, &plaqt);
   wilson_loop(GC, geo, param, size1, size2, &wloops, &wloopt);
   wilson_loop(GC, geo, param, size1-1, size2-1, &wloops_11, &wloopt_11);
   wilson_loop(GC, geo, param, size1-1, size2, &wloops_10, &wloopt_10);
   wilson_loop(GC, geo, param, 4, 4, &wloops_44, &wloopt_44);
   wilson_loop(GC, geo, param, 4, 2, &wloops_42, &wloopt_42);
   polyakov(GC, geo, param, &polyre, &polyim);

   if(fabs(param->d_adj_beta)<MIN_VALUE)      // wilson action
      {
      fprintf(datafilep, "%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g ", 0.5*(plaqs+plaqt), 0.5*(wloops+wloopt), 0.5*(wloops_11+wloopt_11), 0.5*(wloops_10+wloopt_10), 0.5*(wloops_44+wloopt_44), 0.5*(wloops_42+wloopt_42), polyre, polyim);
      fprintf(datafilep, "\n");
      }
   else                      // fundamental plus adjoint action
      {
      double plaqs_adj, plaqt_adj;
      plaquette_adj(GC, geo, param, &plaqs_adj, &plaqt_adj);
      fprintf(datafilep, "%.12g %.12g %.12g %.12g %.12g %.12g %.12g ", 0.5*(plaqs+plaqt), 0.5*(wloops+wloopt), 0.5*(wloops_11+wloopt_11), 0.5*(wloops_10+wloopt_10), 0.5*(plaqs_adj+plaqt_adj), polyre, polyim);
      fprintf(datafilep, "\n");
      }
   fflush(datafilep);
   }





#endif
