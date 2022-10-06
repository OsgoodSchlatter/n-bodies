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
#include <mpi.h>

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

double sum_speed_sq_global = 0;
double max_acc_global = 0;
double max_speed_global = 0;
double max_acc_local = 0;
double max_speed_local = 0;
double sum_speed_sq_local = 0;

int comm_rank, comm_size;

void init()
{
    if (comm_rank==0){
#ifdef DISPLAY
        Display *theDisplay; /* These three variables are required to open the */
        GC theGC;            /* particle plotting window.  They are externally */
        Window theMain;      /* declared in ui.h but are also required here.   */
#endif
    }

}



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

/* compute the new position/velocity Appel => move_particle(&particles[i], step); */
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


  sum_speed_sq_local += speed_sq;
  max_acc_local = MAX(max_acc_local, cur_acc);
  max_speed_local = MAX(max_speed_local, cur_speed);
 // printf("comm_rank %d : max_speed = %f\n",comm_rank,max_speed_local);
}

/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step)
{
  /* First calculate force for particles. */
  int i;
  /*  Séparer les nparticles entre les n_process
      Chaque process a besoin de :
       - particles
       -
  */
  //printf("comm_rank %d : %d <i< %d\n",comm_rank,comm_rank*nparticles/comm_size, (comm_rank+1)*nparticles/comm_size);
  for (i = (int) (comm_rank*nparticles/comm_size) ; i < (int) ((comm_rank+1)*nparticles/comm_size) ; i++)
  {
    //printf("nparticule %d par %d / %d\n",i,comm_rank,comm_size);
    int j;
    particles[i].x_force = 0;
    particles[i].y_force = 0;
    for (j = 0; j < nparticles; j++)
    {
      particle_t *p = &particles[j];
      /* compute the force of particle j on particle i // Ne modifie que particles[i]*/
      compute_force(&particles[i], p->x_pos, p->y_pos, p->mass);
    }
  }

  /* then move all particles and return statistics */
  for (i = (int) (comm_rank*nparticles/comm_size) ; i < (int) ((comm_rank+1)*nparticles/comm_size) ; i++)
  {
    move_particle(&particles[i], step);
  }
}

/* display all the particles */
void draw_all_particles()
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    int x = POS_TO_SCREEN(particles[i].x_pos);
    int y = POS_TO_SCREEN(particles[i].y_pos);
    draw_point(x, y);
  }
}

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
    int nParticulePerProcess = ((int) ((comm_rank+1)*nparticles/comm_size))-((int) ((comm_rank)*nparticles/comm_size));
    //printf("%d : %d",comm_rank,nParticulePerProcess);
    //ALLTOALL particles
    double* my_values = malloc(nParticulePerProcess*6*sizeof(double));
    // Gestion des destinataires et de l'emplacement des valeurs envoyées
    //int* counts_send = malloc(comm_size*sizeof(int));
    //int* displacements_send = malloc(comm_size*sizeof(int));
/*    for(int i = 0; i < comm_size; i++)
    {
        if (comm_rank==comm_size-1){
            counts_send[i]=(nparticles-((int) ((comm_size-1)*nparticles/comm_size))) *6;
        }
        else {
            counts_send[i] = ((int) (nparticles/comm_size))*6;
        }

        displacements_send[i] = (comm_rank*((int) (nparticles/comm_size))) *6;
    }*/
    // Gestion des sources et de la destination des valeurs reçus
    double* buffer_recv= malloc(nparticles*6*sizeof(double));
    int* counts_recv= malloc(comm_size*sizeof(int));
    int* displacements_recv = malloc(comm_size*sizeof(int));
    for(int i = 0; i < comm_size; i++)
    {
        counts_recv[i] = ((int) (nparticles/comm_size)) *6;
        displacements_recv[i] = (i*((int) (nparticles/comm_size))) *6;
    }
    counts_recv[comm_size-1]=(nparticles-((int) ((comm_size-1)*nparticles/comm_size))) *6;
        // AFFICHAGE DES TABLEAU DE GESTION D'ENVOI/RECEPTION
//        MPI_Barrier(MPI_COMM_WORLD);
//        printf("\nComm_rank %d\n",comm_rank);
////        printf("    counts_send \n");
////        for (int i = 0 ;i<comm_size ;i++){
////            printf("%d ",counts_send[i]);
////        }
////        printf("\n    displacements_send \n");
////        for (int i = 0 ;i<comm_size ;i++){
////            printf("%d ",displacements_send[i]);
////        }
//        printf("\n    counts_recv \n");
//        for (int i = 0 ;i<comm_size ;i++){
//            printf("%d ",counts_recv[i]);
//        }
//        printf("\n    displacements_recv \n");
//        for (int i = 0 ;i<comm_size ;i++){
//            printf("%d ",displacements_recv[i]);
//        }
//        printf("\n\n");
//        MPI_Barrier(MPI_COMM_WORLD);

    while (t < T_FINAL && nparticles > 0)
    {
        //printf("comm_rank %d : t = %f / dt = %f / max_acc = %f / max_speed = %f\n",comm_rank,t,dt,max_acc_global,max_speed_global);
//        MPI_Barrier(MPI_COMM_WORLD);
//        for(int i = 0; i < nparticles; i++)
//        {
//            printf("comm_rank %d : t = %f / i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / x_force = %f / y_force = %f\n",
//                   comm_rank,t,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].x_force,particles[i].y_force);
//        }
//        MPI_Barrier(MPI_COMM_WORLD);

        t += dt;

        max_acc_local = max_acc_global;
        max_speed_local = max_speed_global;
        sum_speed_sq_local = 0;

        all_move_particles(dt);

        //init_my_value avec particles x_pos y_pos x_vel y_vel x_force y_force //INIT QUE LA OU IL Y A BESOIN
        int j=(int) ((comm_rank)*nparticles/comm_size);
        for(int i = 0; i < nParticulePerProcess; i++)
        {
            //printf("comm_rank %d : t = %f / i = %d / j = %d\n",comm_rank,t,i,j);
            my_values[i*6] = particles[j].x_pos;
            my_values[i*6+1] = particles[j].y_pos;
            my_values[i*6+2] = particles[j].x_vel;
            my_values[i*6+3] = particles[j].y_vel;
            my_values[i*6+4] = particles[j].x_force;
            my_values[i*6+5] = particles[j].y_force;
            //printf("t %f j %d SEND %f\n",t,j,my_values[6*i+5]);
            j+=1;
        }
        MPI_Allgatherv(my_values,nParticulePerProcess*6, MPI_DOUBLE, buffer_recv, counts_recv, displacements_recv, MPI_DOUBLE, MPI_COMM_WORLD);
        //Mise à jour de particles
        for(int i = 0; i < nparticles; i++)
        {
            particles[i].x_pos =buffer_recv[6*i];
            particles[i].y_pos =buffer_recv[6*i+1];
            particles[i].x_vel =buffer_recv[6*i+2];
            particles[i].y_vel =buffer_recv[6*i+3];
            particles[i].x_force =buffer_recv[6*i+4];
            particles[i].y_force =buffer_recv[6*i+5];
            //printf("t %f i %d RECV %f\n",t,i,buffer_recv[6*i+5]);
        }
        // FIN ALLTOALL particles

        //ALLREDUCTION max_speed / max_acc
        MPI_Allreduce(&max_acc_local,&max_acc_global,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&max_speed_local,&max_speed_global,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        //REDUCTION sum_speed_sq
        MPI_Reduce(&sum_speed_sq_local,&sum_speed_sq_global,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//        MPI_Barrier(MPI_COMM_WORLD);
        if (comm_rank==0){
            dt = 0.1 * max_speed_global / max_acc_global;
            /* Plot the movement of the particle */
#if DISPLAY
            clear_display();
            draw_all_particles();
            flush_display();

#endif
        }
//        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
}

/*
  Simulate the movement of nparticles particles.
*/
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);



  if (argc >= 2)
  {
    nparticles = atoi(argv[1]);
  }
  if (argc == 3)
  {
    T_FINAL = atof(argv[2]);
  }

  init();

  /* Allocate global shared arrays for the particles data set. */
  particles = malloc(sizeof(particle_t) * nparticles);
  //
  double* bcast_buff = malloc(sizeof(double) * 5 * nparticles); //x_pos,y_pos,x_vel,y_vel,mass
  if (comm_rank==0){
      all_init_particles(nparticles, particles);

      for (int i = 0; i < nparticles; i++) {
          bcast_buff[i*5+0]=particles[i].x_pos;
          bcast_buff[i*5+1]=particles[i].y_pos;
          bcast_buff[i*5+2]=particles[i].x_vel;
          bcast_buff[i*5+3]=particles[i].y_vel;
          bcast_buff[i*5+4]=particles[i].mass;
      }
  }
//    MPI_Barrier(MPI_COMM_WORLD);
//    printf("\n Avant Bcast\n");
//    for(int i = 0; i < nparticles; i++)
//    {
//        printf("comm_rank %d :i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / mass = %f\n",
//               comm_rank,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].mass);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(bcast_buff,5 * nparticles,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (comm_rank!=0){
        for (int i = 0; i < nparticles; i++) {
            particles[i].x_pos=bcast_buff[i*5+0];
            particles[i].y_pos=bcast_buff[i*5+1];
            particles[i].x_vel=bcast_buff[i*5+2];
            particles[i].y_vel=bcast_buff[i*5+3];
            particles[i].mass=bcast_buff[i*5+4];
        }
    }

//    MPI_Barrier(MPI_COMM_WORLD);
//    printf("\n Après Bcast\n");
//    for(int i = 0; i < nparticles; i++)
//    {
//        printf("comm_rank %d :i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / mass = %f\n",
//               comm_rank,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].mass);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);


  if (comm_rank==0){
    printf("\n");
    /* Initialize thread data structures */
#ifdef DISPLAY
    /* Open an X window to display the particles */
    simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif
  }

  //Pk pas que process 0
  struct timeval t1, t2;
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank==0){
    gettimeofday(&t1, NULL);
  }

  //printf("Hi from %d / %d\n ",comm_rank,comm_size);
  /* Main thread starts simulation ... */
  run_simulation();

  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank==0){
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

  }

  MPI_Finalize();
  return 0;
}
