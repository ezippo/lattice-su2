#ifndef YM_LOCAL_C
#define YM_LOCAL_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#ifdef OPENMP_MODE
  #include<omp.h>
#endif

#include"../include/gauge_conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void real_main(char *in_file)
    {
    Gauge_Conf GC;
    Geometry geo;
    GParam param;

    int count;
    double acc_metro=0.0;
    FILE *datafilep;
    time_t time1, time2;

    // to disable nested parallelism
    #ifdef OPENMP_MODE
      omp_set_nested(0);
    #endif

    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize geometry
    init_indexing_lexeo();
    init_geometry(&geo, &param);

    // initialize gauge configuration
    init_gauge_conf(&GC, &param);

    // montecarlo
    time(&time1);
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       acc_metro += update_metropolis(&GC, &geo, &param);

       if(count % param.d_measevery ==0 && count >= param.d_thermal)
         {
         perform_measures_localobs(&GC, &geo, &param, datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           write_conf_on_file(&GC, &param);

           // backup copy
           write_conf_on_file_back(&GC, &param);
           }
         }
       }
    time(&time2);
    // montecarlo end

    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &param);
      }

    // print simulation details
    print_parameters_local(&param, time1, time2, acc_metro/param.d_sample);

    // free gauge configuration
    free_gauge_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);
    }

// creates the input_file template
void print_template_input(void)
  {
  int i;
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size ");
    for (i=0; i<STDIM; i++)   fprintf(fp,"4 ");
    fprintf(fp,"\n");
    fprintf(fp, "beta     5.705\n");
    fprintf(fp, "adj_beta 0.0\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "measevery 1\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "\n");
    fprintf(fp, "epsilon_metro  0.1\n");
    fprintf(fp, "hits_metro     1\n");
    fprintf(fp,"\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)   // no input file, prints configuration details and creates the input_file template
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef OPENMP_MODE
        printf("\n\tusing OpenMP with %d threads\n", NTHREADS);
      #endif

      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else    // with input_file, makes the real main with parameters of the input_file
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
