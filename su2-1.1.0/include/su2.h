#ifndef SU2_H
#define SU2_H

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"macro.h"

//
// An Su2 matrix is represented as comp[0]+i\sum_{j=1}^3 comp[j]\sigma_j where
// sigma_j are Pauli matrices, comp[j] are real and \sum_{j=0}^3 comp[j]^2=1
//
typedef struct Su2 {
     double comp[4] __attribute__((aligned(DOUBLE_ALIGN)));
} Su2;



inline void init(Su2 * restrict A, double vec[4])
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=vec[0];
  A->comp[1]=vec[1];
  A->comp[2]=vec[2];
  A->comp[3]=vec[3];
  }


// A=1
inline void one(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=1.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=0
inline void zero(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=0.0;
  A->comp[1]=0.0;
  A->comp[2]=0.0;
  A->comp[3]=0.0;
  }


// A=B
inline void equal(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=B->comp[0];
  A->comp[1]=B->comp[1];
  A->comp[2]=B->comp[2];
  A->comp[3]=B->comp[3];
  }


// A=B^{dag}
inline void equal_dag(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0];
  A->comp[1]=-B->comp[1];
  A->comp[2]=-B->comp[2];
  A->comp[3]=-B->comp[3];
  }


// A+=B
inline void plus_equal(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A+=B^{dag}
inline void plus_equal_dag(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]+=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=B
inline void minus_equal(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]-=B->comp[1];
  A->comp[2]-=B->comp[2];
  A->comp[3]-=B->comp[3];
  }


// A-=(r*B)
inline void minus_equal_times_real(Su2 * restrict A, Su2 const * const restrict B, double r)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=(r*B->comp[0]);
  A->comp[1]-=(r*B->comp[1]);
  A->comp[2]-=(r*B->comp[2]);
  A->comp[3]-=(r*B->comp[3]);
  }


// A-=B^{dag}
inline void minus_equal_dag(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]-=B->comp[0];
  A->comp[1]+=B->comp[1];
  A->comp[2]+=B->comp[2];
  A->comp[3]+=B->comp[3];
  }


// A=b*B+c*C
inline void lin_comb(Su2 * restrict A,
                     double b, Su2 const * const restrict B,
                     double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] + c*C->comp[1];
  A->comp[2]= b*B->comp[2] + c*C->comp[2];
  A->comp[3]= b*B->comp[3] + c*C->comp[3];
  }


// A=b*B^{dag}+c*C
inline void lin_comb_dag1(Su2 * restrict A,
                          double b, Su2 const * const restrict B,
                          double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] + c*C->comp[1];
  A->comp[2]= -b*B->comp[2] + c*C->comp[2];
  A->comp[3]= -b*B->comp[3] + c*C->comp[3];
  }


// A=b*B+c*C^{dag}
inline void lin_comb_dag2(Su2 * restrict A,
                          double b, Su2 const * const restrict B,
                          double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= b*B->comp[0] + c*C->comp[0];
  A->comp[1]= b*B->comp[1] - c*C->comp[1];
  A->comp[2]= b*B->comp[2] - c*C->comp[2];
  A->comp[3]= b*B->comp[3] - c*C->comp[3];
  }


// A=b*B^{dag}+c*C^{dag}
inline void lin_comb_dag12(Su2 * restrict A,
                           double b, Su2 const * const restrict B,
                           double c, Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  b*B->comp[0] + c*C->comp[0];
  A->comp[1]= -b*B->comp[1] - c*C->comp[1];
  A->comp[2]= -b*B->comp[2] - c*C->comp[2];
  A->comp[3]= -b*B->comp[3] - c*C->comp[3];
  }


// A*=r
inline void times_equal_real(Su2 * restrict A, double r)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]*=r;
  A->comp[1]*=r;
  A->comp[2]*=r;
  A->comp[3]*=r;
  }


// A*=B
inline void times_equal(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]= a0*B->comp[0] - a1*B->comp[1] - a2*B->comp[2] - a3*B->comp[3];
  A->comp[1]= a0*B->comp[1] + a1*B->comp[0] - a2*B->comp[3] + a3*B->comp[2];
  A->comp[2]= a0*B->comp[2] + a2*B->comp[0] + a1*B->comp[3] - a3*B->comp[1];
  A->comp[3]= a0*B->comp[3] + a3*B->comp[0] - a1*B->comp[2] + a2*B->comp[1];
  }


// A*=B^{dag}
inline void times_equal_dag(Su2 * restrict A, Su2 const * const restrict B)
  {
  #ifdef DEBUG
  if(A==B)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
  #endif

  register double a0, a1, a2, a3;

  a0=A->comp[0];
  a1=A->comp[1];
  a2=A->comp[2];
  a3=A->comp[3];

  A->comp[0]=  a0*B->comp[0] + a1*B->comp[1] + a2*B->comp[2] + a3*B->comp[3];
  A->comp[1]= -a0*B->comp[1] + a1*B->comp[0] + a2*B->comp[3] - a3*B->comp[2];
  A->comp[2]= -a0*B->comp[2] + a2*B->comp[0] - a1*B->comp[3] + a3*B->comp[1];
  A->comp[3]= -a0*B->comp[3] + a3*B->comp[0] + a1*B->comp[2] - a2*B->comp[1];
  }


// A=B*C
inline void times(Su2 * restrict A,
                  Su2 const * const restrict B,
                  Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C
inline void times_dag1(Su2 * restrict A,
                       Su2 const * const restrict B,
                       Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]= B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B*C^{dag}
inline void times_dag2(Su2 * restrict A,
                       Su2 const * const restrict B,
                       Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] + B->comp[1]*C->comp[1] + B->comp[2]*C->comp[2] + B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] + B->comp[1]*C->comp[0] + B->comp[2]*C->comp[3] - B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] + B->comp[2]*C->comp[0] - B->comp[1]*C->comp[3] + B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] + B->comp[3]*C->comp[0] + B->comp[1]*C->comp[2] - B->comp[2]*C->comp[1];
  }


// A=B^{dag}*C^{dag}
inline void times_dag12(Su2 * restrict A,
                        Su2 const * const restrict B,
                        Su2 const * const restrict C)
  {
  #ifdef DEBUG
  if(A==B || A==C || B==C)
    {
    fprintf(stderr, "The same pointer is used twice in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
    __assume_aligned(&(B->comp), DOUBLE_ALIGN);
    __assume_aligned(&(C->comp), DOUBLE_ALIGN);
  #endif

  A->comp[0]=  B->comp[0]*C->comp[0] - B->comp[1]*C->comp[1] - B->comp[2]*C->comp[2] - B->comp[3]*C->comp[3];
  A->comp[1]= -B->comp[0]*C->comp[1] - B->comp[1]*C->comp[0] - B->comp[2]*C->comp[3] + B->comp[3]*C->comp[2];
  A->comp[2]= -B->comp[0]*C->comp[2] - B->comp[2]*C->comp[0] + B->comp[1]*C->comp[3] - B->comp[3]*C->comp[1];
  A->comp[3]= -B->comp[0]*C->comp[3] - B->comp[3]*C->comp[0] - B->comp[1]*C->comp[2] + B->comp[2]*C->comp[1];
  }


// random SU(2) matrix
void rand_matrix(Su2 *A);


// random SU(2) matrix with p0 given (used in the update)
void rand_matrix_p0(double p0, Su2 *A);


// sqrt of the determinant
inline double sqrtdet(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrt(A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1]
             + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3]);
  }


// l2 norm of the matrix
inline double norm(Su2 const * const restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  return sqrtdet(A);
  }


// real part of the trace /2
inline double retr(Su2 const * const restrict A)
  {
  return A->comp[0];
  }


// imaginary part of the trace /2
inline double imtr(Su2 const * const restrict A)
  {
  (void) A; // to suppress compilation warning of unused variable
  return 0.0;
  }


// unitarize the matrix
inline void unitarize(Su2 * restrict A)
  {
  #ifdef __INTEL_COMPILER
    __assume_aligned(&(A->comp), DOUBLE_ALIGN);
  #endif

  double p;

  p=A->comp[0]*A->comp[0] + A->comp[1]*A->comp[1] + A->comp[2]*A->comp[2] + A->comp[3]*A->comp[3];
  p=1.0/sqrt(p);

  A->comp[0]*=p;
  A->comp[1]*=p;
  A->comp[2]*=p;
  A->comp[3]*=p;
  }


// print on screen
void print_on_screen(Su2 const * const A);


// print on file
int print_on_file(FILE *fp, Su2 const * const A);


// print on binary file without changing endiannes
int print_on_binary_file_noswap(FILE *fp, Su2 const * const A);


// print on binary file changing endiannes
int print_on_binary_file_swap(FILE *fp, Su2 const * const A);


// print on binary file in big endian format
int print_on_binary_file_bigen(FILE *fp, Su2 const * const A);


// read from file
int read_from_file(FILE *fp, Su2 *A);


// read from binary file without changing endiannes
int read_from_binary_file_noswap(FILE *fp, Su2 *A);


// read from binary file changing endiannes
int read_from_binary_file_swap(FILE *fp, Su2 *A);


// read from binary file written in big endian
int read_from_binary_file_bigen(FILE *fp, Su2 *A);


#endif
