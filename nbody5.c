/*Retirer 2 iéme boucle dans move particle
 * retirer un accès au tableau 
  * retirer la puissance pour un x*x*x
  * 
  */ 
 #include <immintrin.h>
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

  f32 *x, *y, *z;
  f32 *vx, *vy, *vz;
  
} particle_t;

static inline f32 _hadd_ps_mm256(const __m256 a) {
  __m256 t1 = _mm256_hadd_ps(a, a);
  __m256 t2 = _mm256_hadd_ps(t1, t1);
  __m128 t3 = _mm256_extractf128_ps(t2, 1);
  __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2), t3);
  return _mm_cvtss_f32(t4);
}
//
void init(particle_t *p, u64 n)
{

    p->vx = malloc(sizeof(f32) * n);
    p->vy = malloc(sizeof(f32) * n);
    p->vz = malloc(sizeof(f32) * n);

    p->x = malloc(sizeof(f32) * n);
    p->y = malloc(sizeof(f32) * n);
    p->z = malloc(sizeof(f32) * n);

  for (u64 i = 0; i < n; i++)
    {
      //
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p->x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->y[i] = (f32)rand() / (f32)RAND_MAX;
      p->z[i] = sign * (f32)rand() / (f32)RAND_MAX;

      //
      p->vx[i] = (f32)rand() / (f32)RAND_MAX;
      p->vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p->vz[i] = (f32)rand() / (f32)RAND_MAX;
    }
}

//
void move_particles(particle_t *restrict p, const f32 dt,const u64 n)
{
  //
   const f32 softening = 1e-20;
    __m256 softeninga = _mm256_set1_ps(softening);
    __m256 dta = _mm256_set1_ps(dt);

  for (u64 i = 0; i < n; ++i)
    {
      //
       __m256 fx = _mm256_setzero_ps();
       __m256 fy = _mm256_setzero_ps();
       __m256 fz = _mm256_setzero_ps();
     __m256 dxia = _mm256_set1_ps(p->x[i]);
    __m256 dyia = _mm256_set1_ps(p->y[i]);
    __m256 dzia = _mm256_set1_ps(p->z[i]);

      for (u64 j = 0; j < n; j= j+8)
	{   


	  //Newton's law
    __m256 dx = _mm256_sub_ps(_mm256_loadu_ps(p->x+j),dxia);
    __m256 dy = _mm256_sub_ps(_mm256_loadu_ps(p->y+j),dyia);
    __m256 dz = _mm256_sub_ps(_mm256_loadu_ps(p->z+j),dzia);

    __m256 test = _mm256_fmadd_ps(dz, dz,_mm256_fmadd_ps(dy, dy,_mm256_fmadd_ps(dx, dx, softeninga)));

    __m256 d_2 = _mm256_rsqrt_ps (test);
   __m256 d_3_over_2=  _mm256_mul_ps(_mm256_mul_ps(d_2,d_2),d_2);
	  //Net force

          fx = _mm256_fmadd_ps(dx,d_3_over_2, fx);
      fy = _mm256_fmadd_ps(dy, d_3_over_2, fy);
      fz = _mm256_fmadd_ps(dz, d_3_over_2, fz);


	}

    p->vx[i] += dt * _hadd_ps_mm256(fx);
    p->vy[i] += dt * _hadd_ps_mm256(fy);
    p->vz[i] += dt * _hadd_ps_mm256(fz);  


    }

  //3 floating-point operations
  for (u64 i = 0; i < n; i=i+8)
    {

      __m256 tmpx = _mm256_fmadd_ps(dta,_mm256_loadu_ps(p->vx+i),_mm256_loadu_ps(p->x+i));
     __m256 tmpy = _mm256_fmadd_ps(dta,_mm256_loadu_ps(p->vy+i),_mm256_loadu_ps(p->y+i));
      __m256 tmpz  = _mm256_fmadd_ps(dta,_mm256_loadu_ps(p->vz+i),_mm256_loadu_ps(p->z+i));

      _mm256_storeu_ps(p->x + i, tmpx);
    _mm256_storeu_ps(p->y + i, tmpy);
    _mm256_storeu_ps(p->z + i, tmpz); 

     /* printf("%f %f %f\n",p->x[i],p->y[i],p->z[i]);      
      printf("%f %f %f\n",p->x[i+1],p->y[i+1],p->z[i+1]);      
      printf("%f %f %f\n",p->x[i+2],p->y[i+2],p->z[i+2]);      
      printf("%f %f %f\n",p->x[i+3],p->y[i+3],p->z[i+3]);      
      printf("%f %f %f\n",p->x[i+4],p->y[i+4],p->z[i+4]);      
      printf("%f %f %f\n",p->x[i+5],p->y[i+5],p->z[i+5]);      
      printf("%f %f %f\n",p->x[i+6],p->y[i+6],p->z[i+6]);      
      printf("%f %f %f\n",p->x[i+7],p->y[i+7],p->z[i+7]);      */


    }


}

//
int main(int argc, char **argv)
{
    //srand(0);

  
  //
  const u64 n = 16384 ;
  const u64 steps= 10;
  const f32 dt = 0.01;

  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t p ;// malloc(sizeof(particle_t) * n);

  //
  init(&p, n);

  const u64 s = sizeof(f32) * n * 6;
  
 printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  //
  for (u64 i = 0; i < steps; i++)
    {

      //Measure
      const f64 start = omp_get_wtime();

      move_particles(&p, dt, n);

      const f64 end = omp_get_wtime();

      //Number of interactions/iterations
      const f32 h1 = (f32)(n) * (f32)(n-1);

      //GFLOPS
      const f32 h2 = (20 * h1 + 6 * (f32)n) * 1e-9;
      
      if (i >= warmup)
	{
	  rate += h2 / (end - start);
	  drate += (h2 * h2) / ((end - start) * (end - start));
	}

      //
      printf("%5llu %10.3e %10.3e %8.1f %s\n",
	     i,
	     (end - start),
	     h1 / (end - start),
	     h2 / (end - start),
	     (i < warmup) ? "*" : "");
      
      fflush(stdout);
    }

  //
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, drate);
  printf("-----------------------------------------------------\n");
  
  //
  free(p.vx);
  free(p.vy);
  free(p.vz);
  free(p.x);
  free(p.y);
  free(p.z);


  //
  return 0;
}
