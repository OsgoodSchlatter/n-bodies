/*
** nbody_brute_force.c - nbody simulation using the brute-force algorithm (O(n*n))
**
**/

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include <helper_cuda.h>

#ifdef DISPLAY
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#include "ui.h"
#include "nbody.h"
#include "nbody_tools.h"

FILE *f_out = NULL;

int nparticles = 1500; /* number of particles */
float T_FINAL = 1.0;   /* simulation end time */
particle_t *particles;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

void init()
{
#ifdef DISPLAY
    Display *theDisplay; /* These three variables are required to open the */
    GC theGC;            /* particle plotting window.  They are externally */
    Window theMain;      /* declared in ui.h but are also required here.   */
#endif
}


/********************** kernel **************************/
__global__ void kernel_compute_force(particle_t *particles, int nparticles)
{
    /* Calcul de l'indice i*/
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;

    if (i<nparticles && j<nparticles) {
        //LECTURE GPU GLOBAL

        //particle_t* pi = &particles[i];
        // pi.x_force = 0;
        // pi.y_force = 0;
        //
        // //LECTURE GPU GLOBAL
        // particle_t* pj = &particles[j];
        // double x_sep, y_sep, dist_sq, grav_base;
        //
        // x_sep = pj->x_pos - pi->x_pos;
        // y_sep = pj->y_pos - pi->y_pos;
        // dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);
        //
        // /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
        // grav_base = GRAV_CONSTANT * (pi->mass) * (pj->mass) / dist_sq;
        //
        // //ECRITURE LOCAL -> REGISTRE
        // pi->x_force += grav_base * x_sep;
        // pi->y_force += grav_base * y_sep;
        //
        // //ecriture Global i
        // particles[i] = pi;

        particles[i].x_force = 0;
        particles[i].y_force = 0;

        double x_sep, y_sep, dist_sq, grav_base;

        x_sep = particles[j].x_pos - particles[i].x_pos;
        y_sep = particles[j].y_pos - particles[i].y_pos;
        dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);

        /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
        grav_base = GRAV_CONSTANT * (particles[i].mass) * (particles[j].mass) / dist_sq;

        //RENDRE ATOMIC
        particles[i].x_force += grav_base * x_sep;
        particles[i].y_force += grav_base * y_sep;

    }
}

// __global__ void kernel_move_particle(particle_t *pTab, int nparticles, double step)
// {
//     //HostToDevice => pTab, speed_sq, cur_acc, cur_speed
//
//     /* Calcul de l'indice i*/
//     int i = blockIdx.x*blockDim.x + threadIdx.x;
//     if (i<nparticles){
//         particle_t* p = &pTab[i];
//
//         p->x_pos += (p->x_vel) * step;
//         p->y_pos += (p->y_vel) * step;
//         double x_acc = p->x_force / p->mass;
//         double y_acc = p->y_force / p->mass;
//         p->x_vel += x_acc * step;
//         p->y_vel += y_acc * step;
//
//         /* compute statistics */
//         double cur_acc = (x_acc * x_acc + y_acc * y_acc);
//         cur_acc = sqrt(cur_acc);
//         double speed_sq = (p->x_vel) * (p->x_vel) + (p->y_vel) * (p->y_vel);
//         double cur_speed = sqrt(speed_sq);
//
//         //BARRIER
//         //sum_speed_sq += speed_sq;
//         //max_acc = MAX(max_acc, cur_acc);
//         //max_speed = MAX(max_speed, cur_speed);
//         //BARRIER
//
//         //DeviceToHost => pTab, speed_sq, cur_acc, cur_speed
//         pTab[i]=p;
//     }
// }
// __global__ void kernel_all_move_particles(particle_t *particles, double step,int nparticles)
// {
//     /* Calcul de l'indice i*/
//     int i = blockIdx.x*blockDim.x + threadIdx.x;
//     if (i<nparticles){
//         //LECTURE GPU GLOBAL
//         particle_t* pi = &particles[i];
//         int j;
//         pi.x_force = 0;
//         pi.y_force = 0;
//         for (j = 0; j < nparticles; j++)
//         {
//             //LECTURE GPU GLOBAL
//             particle_t *pj = &particles[j];
//             double x_sep, y_sep, dist_sq, grav_base;
//
//             x_sep = pj->x_pos - pi->x_pos;
//             y_sep = pj->y_pos - pi->y_pos;
//             dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);
//
//             /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
//             grav_base = GRAV_CONSTANT * (pi->mass) * (pj->mass) / dist_sq;
//
//             //ECRITURE LOCAL -> REGISTRE
//             pi->x_force += grav_base * x_sep;
//             pi->y_force += grav_base * y_sep;
//         }
//
//         pi>x_pos += (pi->x_vel) * step;
//         pi->y_pos += (pi->y_vel) * step;
//         double x_acc = pi->x_force / p->mass;
//         double y_acc = pi->y_force / p->mass;
//         p->x_vel += x_acc * step;
//         p->y_vel += y_acc * step;
//
//         /* compute statistics */
//         cur_acc = (x_acc * x_acc + y_acc * y_acc);
//         cur_acc = sqrt(cur_acc);
//         speed_sq = (p->x_vel) * (p->x_vel) + (p->y_vel) * (p->y_vel);
//         cur_speed = sqrt(speed_sq);
//
//         //Barriere
//         //ECRITURE GPU GLOBALE mais chacun son slot
//         particles[i]=pi;
//
//
//     }
// }



/* compute the force that a particle with position (x_pos, y_pos) and mass 'mass'
 * applies to particle p
 */
void compute_force(particle_t *p, double x_pos, double y_pos, double mass)
{
  double x_sep, y_sep, dist_sq, grav_base;

  x_sep = x_pos - p->x_pos;
  y_sep = y_pos - p->y_pos;
  dist_sq = MAX((x_sep * x_sep) + (y_sep * y_sep), 0.01);

  /* Use the 2-dimensional gravity rule: F = d * (GMm/d^2) */
  grav_base = GRAV_CONSTANT * (p->mass) * (mass) / dist_sq;

  p->x_force += grav_base * x_sep;
  p->y_force += grav_base * y_sep;
}

/* compute the new position/velocity */
void move_particle(particle_t *p, double step)
{

  p->x_pos += (p->x_vel) * step;
  p->y_pos += (p->y_vel) * step;
  double x_acc = p->x_force / p->mass;
  double y_acc = p->y_force / p->mass;
  p->x_vel += x_acc * step;
  p->y_vel += y_acc * step;

  /* compute statistics */
  double cur_acc = (x_acc * x_acc + y_acc * y_acc);
  cur_acc = sqrt(cur_acc);
  double speed_sq = (p->x_vel) * (p->x_vel) + (p->y_vel) * (p->y_vel);
  double cur_speed = sqrt(speed_sq);

  sum_speed_sq += speed_sq;
  max_acc = MAX(max_acc, cur_acc);
  max_speed = MAX(max_speed, cur_speed);
}

/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step)
{
  /* First calculate force for particles. */
  //nparticles*nparticles*1 calculs
  int i;
  // for (i = 0; i < nparticles; i++)
  // {
  //   int j;
  //   particles[i].x_force = 0;
  //   particles[i].y_force = 0;
  //   for (j = 0; j < nparticles; j++)
  //   {
  //     particle_t *p = &particles[j];
  //     /* compute the force of particle j on particle i */
  //     compute_force(&particles[i], p->x_pos, p->y_pos, p->mass);
  //   }
  //run_simulation
  //
  // }
  int blockSize = 32;
  dim3 dimBlock (blockSize,blockSize,1);
  dim3 dimGrid (ceil(nparticles/blockSize),ceil(nparticles/blockSize),1);
  kernel_compute_force<<<dimGrid, dimBlock>>>(particles,nparticles);
  cudaDeviceSynchronize();
  //Barrier
    printf("pos_x = %f \n",particles[0].x_pos);

  //nparticules*1*1 calculs
  /* then move all particles and return statistics */
  for (i = 0; i < nparticles; i++)
  {
    move_particle(&particles[i], step);
  }
  //Calcul du max et sum
}

// /* display all the particles */
// void draw_all_particles()
// {
//   int i;
//   for (i = 0; i < nparticles; i++)
//   {
//     int x = POS_TO_SCREEN(particles[i].x_pos);
//     int y = POS_TO_SCREEN(particles[i].y_pos);
//     draw_point(x, y);
//   }
// }

void print_all_particles(FILE *f)
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    particle_t *p = &particles[i];
    fprintf(f, "particle={pos=(%f,%f), vel=(%f,%f)}\n", p->x_pos, p->y_pos, p->x_vel, p->y_vel);
  }
}

void run_simulation()
{
  double t = 0.0, dt = 0.01;

  while (t < T_FINAL && nparticles > 0)
  {
      printf("t = %f\n",t);
    /* Update time. */
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt);

//    /* Adjust dt based on maximum speed and acceleration--this
//       simple rule tries to insure that no velocity will change
//       by more than 10% */
//      for(int i = 0; i < nparticles; i++)
//      {
//          printf("t = %f / i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / x_force = %f / y_force = %f\n",
//                 t,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].x_force,particles[i].y_force);
//      }

    dt = 0.1 * max_speed / max_acc;

    /* Plot the movement of the particle */
#if DISPLAY
    clear_display();
    draw_all_particles();
    flush_display();
#endif
  }
}

/*
  Place particles in their initial positions.
*/
void all_init_particles(int num_particles, particle_t *particles)
{
    int    i;
    double total_particle = num_particles;

    for (i = 0; i < num_particles; i++) {
        particle_t *particle = &particles[i];
#if 0
        particle->x_pos = ((rand() % max_resolution)- (max_resolution/2))*2.0 / max_resolution;
    particle->y_pos = ((rand() % max_resolution)- (max_resolution/2))*2.0 / max_resolution;
    particle->x_vel = particle->y_pos;
    particle->y_vel = particle->x_pos;
#else
        particle->x_pos = i*2.0/nparticles - 1.0;
        particle->y_pos = 0.0;
        particle->x_vel = 0.0;
        particle->y_vel = particle->x_pos;
#endif
        particle->mass = 1.0 + (num_particles+i)/total_particle;
        particle->node = NULL;

        //insert_particle(particle, root);
    }
}
/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char **argv)
{
  if (argc >= 2)
  {
    nparticles = atoi(argv[1]);
  }
  if (argc == 3)
  {
    T_FINAL = atof(argv[2]);
  }

  init();

  //Allocate managed memory
  cudaMallocManaged(&particles,sizeof(particle_t) * nparticles);
  /* Allocate global shared arrays for the particles data set. */
  //particles = malloc(sizeof(particle_t) * nparticles);
  all_init_particles(nparticles, particles);

  /* Initialize thread data structures */
#ifdef DISPLAY
  /* Open an X window to display the particles */
  simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

  struct timeval t1, t2;
  gettimeofday(&t1, NULL);


  /* Main thread starts simulation ... */
  run_simulation();



  gettimeofday(&t2, NULL);

  double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

#ifdef DUMP_RESULT
  FILE *f_out = fopen("particles.log", "w");
  assert(f_out);
  print_all_particles(f_out);
  fclose(f_out);
#endif

  printf("-----------------------------\n");
  printf("nparticles: %d\n", nparticles);
  printf("T_FINAL: %f\n", T_FINAL);
  printf("-----------------------------\n");
  printf("Simulation took %lf s to complete\n", duration);

#ifdef DISPLAY
  clear_display();
  draw_all_particles();
  flush_display();

  printf("Hit return to close the window.");

  getchar();
  /* Close the X window used to display the particles */
  XCloseDisplay(theDisplay);
#endif
  //FREE
  cudaFree(particles);
  return 0;
}
