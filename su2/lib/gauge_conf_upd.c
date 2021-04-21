#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#ifdef OPENMP_MODE
  #include<omp.h>
#endif
#include<stdlib.h>

#include"../include/gauge_conf.h"
#include"../include/gparam.h"
#include"../include/random.h"
#include"../include/su2_upd.h"

// compute the staple in position r, direction i and save it in M
void calcstaples_wilson(Gauge_Conf const * const GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i,
                        Su2 *M)
  {
  int j, l;
  long k;
  Su2 link1, link2, link3, link12, stap;

  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #else
  (void) param; // just to avoid warnings
  #endif

  zero(M); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     plus_equal(M, &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     plus_equal(M, &stap);
     }
   }


// compute the components of the staple in position r, direction i and save it in M[2*(STDIM-1)+1]
// position 0 of M is not used.
void calcstaples_wilson_nosum(Gauge_Conf const * const GC,
                              Geometry const * const geo,
                              GParam const * const param,
                              long r,
                              int i,
                              Su2 *M)
  {
  int j, l, count;
  long k;
  Su2 link1, link2, link3, link12, stap;

  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #else
  (void) param; // just to avoid warnings
  #endif

  count=0;
  zero(&M[count]); // M=0

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     count++;
     zero(&M[count]); // M=0

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

     equal(&link1, &(GC->lattice[nnp(geo, r, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[nnp(geo, r, j)][i]));  // link2 = (2)
     equal(&link3, &(GC->lattice[r][j]));               // link3 = (3)

     times_dag2(&link12, &link1, &link2);  // link12=link1*link2^{dag}
     times_dag2(&stap, &link12, &link3);   // stap=link12*stap^{dag}

     equal(&(M[count]), &stap);

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

     count++;
     zero(&M[count]); // M=0

     k=nnm(geo, r, j);

     equal(&link1, &(GC->lattice[nnp(geo, k, i)][j]));  // link1 = (1)
     equal(&link2, &(GC->lattice[k][i]));               // link2 = (2)
     equal(&link3, &(GC->lattice[k][j]));               // link3 = (3)

     times_dag12(&link12, &link1, &link2); // link12=link1^{dag}*link2^{dag}
     times(&stap, &link12, &link3);        // stap=link12*link3

     equal(&(M[count]), &stap);
     }

   #ifdef DEBUG
   Su2 helper;
   int m;

   calcstaples_wilson(GC, geo, param, r, i, &helper);
   for(m=0; m<2*(STDIM-1)+1; m++)
      {
      minus_equal(&helper, &(M[m]));
      }
   if(norm(&helper)>MIN_VALUE)
     {
     fprintf(stderr, "Problems in calcstaples_wilson_nosum (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }
   #endif
   }


// perform an update with heatbath
void heatbath(Gauge_Conf *GC,
              Geometry const * const geo,
              GParam const * const param,
              long r,
              int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  Su2 stap;

  calcstaples_wilson(GC, geo, param, r, i, &stap);

  times_equal_real(&stap, param->d_beta);
  single_heatbath(&(GC->lattice[r][i]), &stap);
  }


// perform an update with overrelaxation
void overrelaxation(Gauge_Conf *GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    int i)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif
  (void) param; // just to avoid wornings

  Su2 stap;

  calcstaples_wilson(GC, geo, param, r, i, &stap);

  single_overrelaxation(&(GC->lattice[r][i]), &stap);
  }


// perform an update with metropolis
// return 1 if the proposed update is accepted
int metropolis(Gauge_Conf *GC,
             Geometry const * const geo,
             GParam const * const param,
             long r,
             int i,
             int numhits)
{
#ifdef DEBUG
if(r >= param->d_volume)
  {
  fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }
if(i >= STDIM)
  {
  fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
  exit(EXIT_FAILURE);
  }
#endif

Su2 stap, new_link, tmp_matrix, rnd_matrix;
double action_new, action_old;
int acc, hits;

calcstaples_wilson(GC, geo, param, r, i, &stap);

acc=0;

for(hits=0; hits<numhits; hits++)
   {
   // compute old action
   times(&tmp_matrix, &(GC->lattice[r][i]), &stap);
   action_old=param->d_beta*(1.0-retr(&tmp_matrix));

   // compute the new link
   rand_matrix_p0(param->d_epsilon_metro, &rnd_matrix);   // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
   if(casuale()<0.5)
     {
     times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
     }
   else
     {
     times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
     }

   // new action
   times(&tmp_matrix, &new_link, &stap);
   action_new=param->d_beta*(1.0-retr(&tmp_matrix));

   if(casuale()< exp(action_old-action_new))
     {
     equal(&(GC->lattice[r][i]), &new_link);
     acc+=1;
     }
   else
     {
     acc+=0;
     }
   }

return acc;
}


// perform an update with metropolis
// action used: wilson + adjoint
// return 1 if the proposed update is accepted
int metropolis_fund_plus_adj(Gauge_Conf *GC,
                             Geometry const * const geo,
                             GParam const * const param,
                             long r,
                             int i,
                             int numhits)
  {
  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  Su2 stap[2*(STDIM-1)+1], new_link, tmp_matrix, rnd_matrix;
  double action_new, action_old, tmp;
  int acc, hits, i_stap;

  calcstaples_wilson_nosum(GC, geo, param, r, i, stap);

  acc=0;

  for(hits=0; hits<numhits; hits++)
     {
     // compute old action
     for(i_stap=1; i_stap<2*(STDIM-1)+1; i_stap++)
      {
       times(&tmp_matrix, &(GC->lattice[r][i]), &stap[i_stap]);
       tmp = 2.0*retr(&tmp_matrix);
       action_old=-param->d_beta*(tmp*0.5);    // wilson term
       action_old-=param->d_adj_beta*(tmp*tmp/3.0);   // adjoint term
       tmp = 2.0*imtr(&tmp_matrix);
       action_old-=param->d_adj_beta*(tmp*tmp)/3.0;  // adjoint term
     }

     // compute the new link
     rand_matrix_p0(param->d_epsilon_metro, &rnd_matrix);   // rnd_matrix = Proj_on_the_group[ 1 + epsilon_metro*random_matrix ]
     if(casuale()<0.5)
       {
       times(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }
     else
       {
       times_dag1(&new_link, &rnd_matrix, &(GC->lattice[r][i]));
       }

     // new action
     for(i_stap=1; i_stap<2*(STDIM-1)+1; i_stap++)
       {
       times(&tmp_matrix, &new_link, &stap[i_stap]);
       tmp = 2.0*retr(&tmp_matrix);
       action_new=-param->d_beta*(tmp*0.5);    // wilson term
       action_new-=param->d_adj_beta*(tmp*tmp/3.0);   // adjoint term
       tmp = 2.0*imtr(&tmp_matrix);
       action_new-=param->d_adj_beta*(tmp*tmp)/3.0;  // adjoint term
       }

     if(casuale()< exp(action_old-action_new))
       {
       equal(&(GC->lattice[r][i]), &new_link);
       acc+=1;
       }
     else
       {
       acc+=0;
       }
     }

  return acc;
  }


// perform a complete update using Heatbath and Over relaxation
void update(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param)
   {
   for(int i=0; i<STDIM; i++)
      {
      if(param->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long r;
   int j, dir;
   const int over=5;

   // heatbath
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=0; r<(param->d_volume)/2; r++)
         {
         heatbath(GC, geo, param, r, dir);
         }

      #ifdef OPENMP_MODE
      #pragma omp parallel for num_threads(NTHREADS) private(r)
      #endif
      for(r=(param->d_volume)/2; r<(param->d_volume); r++)
         {
         heatbath(GC, geo, param, r, dir);
         }
      }

   // overrelax
   for(dir=0; dir<STDIM; dir++)
      {
      #ifdef THETA_MODE
      compute_clovers(GC, geo, param, dir);
      #endif

      for(j=0; j<over; j++)
         {
         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=0; r<(param->d_volume)/2; r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }

         #ifdef OPENMP_MODE
         #pragma omp parallel for num_threads(NTHREADS) private(r)
         #endif
         for(r=(param->d_volume)/2; r<(param->d_volume); r++)
            {
            overrelaxation(GC, geo, param, r, dir);
            }
         }
      }

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }

   GC->update_index++;
   }

// perform a complete update using Metropolis
double update_metropolis(Gauge_Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            int numhits)
   {
   for(int i=0; i<STDIM; i++)
      {
      if(param->d_size[i]==1)
        {
        fprintf(stderr, "Error: this functon can not be used in the completely reduced case (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
      }

   long r;
   double acc=0.0;
   int dir;

   if(fabs(param->d_adj_beta)<MIN_VALUE)    // Wilson action
     {
     // metropolis
     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef THETA_MODE
        compute_clovers(GC, geo, param, dir);
        #endif

        #ifdef OPENMP_MODE
        #pragma omp parallel for reduction(+:acc) num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<(param->d_volume)/2; r++)
           {
           acc += metropolis(GC, geo, param, r, dir, numhits);
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for reduction(+:acc) num_threads(NTHREADS) private(r)
        #endif
        for(r=(param->d_volume)/2; r<(param->d_volume); r++)
           {
           acc += metropolis(GC, geo, param, r, dir, numhits);
           }
         }
     }
   else    // Wilson plus adjoint action
     {
     // metropolis
     for(dir=0; dir<STDIM; dir++)
        {
        #ifdef THETA_MODE
        compute_clovers(GC, geo, param, dir);
        #endif

        #ifdef OPENMP_MODE
        #pragma omp parallel for reduction(+:acc) num_threads(NTHREADS) private(r)
        #endif
        for(r=0; r<(param->d_volume)/2; r++)
           {
           acc += metropolis(GC, geo, param, r, dir, numhits);
           }

        #ifdef OPENMP_MODE
        #pragma omp parallel for reduction(+:acc) num_threads(NTHREADS) private(r)
        #endif
        for(r=(param->d_volume)/2; r<(param->d_volume); r++)
           {
           acc += metropolis(GC, geo, param, r, dir, numhits);
           }
         }
       }
    acc /= (double)(numhits*param->d_volume*STDIM);

   // final unitarization
   #ifdef OPENMP_MODE
   #pragma omp parallel for num_threads(NTHREADS) private(r, dir)
   #endif
   for(r=0; r<(param->d_volume); r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         unitarize(&(GC->lattice[r][dir]));
         }
      }

   GC->update_index++;
   return acc;
   }


#endif
