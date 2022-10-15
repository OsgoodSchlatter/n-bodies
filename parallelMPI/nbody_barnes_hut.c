/*
** nbody_barnes_hut.c - nbody simulation that implements the Barnes-Hut algorithm (O(nlog(n)))
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

node_t *root;

int comm_rank, comm_size;

double sum_speed_sq = 0;
double max_acc = 0;
double max_speed = 0;

double sum_speed_sq_global = 0;
double max_acc_global = 0;
double max_speed_global = 0;
double max_acc_local = 0;
double max_speed_local = 0;
double sum_speed_sq_local = 0;

int nParticulePerProcess;
//printf("%d : %d",comm_rank,nParticulePerProcess);
//ALLGATHERV particles
int n_caracteristic_shared;
double* my_values;
// Gestion des sources et de la destination des valeurs reçus
double* buffer_recv;
int* counts_recv;
int* displacements_recv ;
int actual_n_particule=0;


void init()
{
  init_alloc(8 * nparticles);
  root = malloc(sizeof(node_t));
  init_node(root, NULL, XMIN, XMAX, YMIN, YMAX);
#ifdef DISPLAY
  Display *theDisplay; /* These three variables are required to open the */
  GC theGC;            /* particle plotting window.  They are externally */
  Window theMain;      /* declared in ui.h but are also required here.   */
#endif
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

/* compute the force that node n acts on particle p */
void compute_force_on_particle(node_t *n, particle_t *p)
{
  if (!n || n->n_particles == 0)
  {
    return;
  }
  if (n->particle)
  {
    /* only one particle */
    assert(n->children == NULL);

    /*
      If the current node is an external node (and it is not body b),
      calculate the force exerted by the current node on b, and add
      this amount to b's net force.
    */
    compute_force(p, n->x_center, n->y_center, n->mass);
  }
  else
  {
    /* There are multiple particles */

#define THRESHOLD 2
    double size = n->x_max - n->x_min; // width of n
    double diff_x = n->x_center - p->x_pos;
    double diff_y = n->y_center - p->y_pos;
    double distance = sqrt(diff_x * diff_x + diff_y * diff_y);

#if BRUTE_FORCE
    /*
      Run the procedure recursively on each of the current
      node's children.
      --> This result in a brute-force computation (complexity: O(n*n))
    */
    int i;
    for (i = 0; i < 4; i++)
    {
      compute_force_on_particle(&n->children[i], p);
    }
#else
    /* Use the Barnes-Hut algorithm to get an approximation */
    if (size / distance < THRESHOLD)
    {
      /*
  The particle is far away. Use an approximation of the force
      */
      compute_force(p, n->x_center, n->y_center, n->mass);
    }
    else
    {
      /*
        Otherwise, run the procedure recursively on each of the current
  node's children.
      */
      int i;
      for (i = 0; i < 4; i++)
      {
        compute_force_on_particle(&n->children[i], p);
      }
    }
#endif
  }
}

void compute_force_in_node(node_t *n)
{
  if (!n)
    return;

  if (n->particle)
  {
    particle_t *p = n->particle;
    p->x_force = 0;
    p->y_force = 0;
    compute_force_on_particle(root, p);
  }
  if (n->children)
  {
    int i;
    for (i = 0; i < 4; i++)
    {
            compute_force_in_node(&n->children[i]);
    }
  }
}

/* compute the new position/velocity */
void move_particle(particle_t *p, double step, node_t *new_root)
{

  assert(p->node != NULL);
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

  p->node = NULL;
  if (p->x_pos < new_root->x_min ||
      p->x_pos > new_root->x_max ||
      p->y_pos < new_root->y_min ||
      p->y_pos > new_root->y_max)
  {
    nparticles--;
  }
  else
  {
    insert_particle(p, new_root);
  }
}

/* compute the new position of the particles in a node */
void move_particles_in_node(node_t *n, double step, node_t *new_root)
{
  if (!n)
    return;

  if (n->particle)
  {
    particle_t *p = n->particle;
    move_particle(p, step, new_root);
  }
  if (n->children)
  {
    int i;
    for (i = 0; i < 4; i++)
    {
        //if (comm_rank==i){
            move_particles_in_node(&n->children[i], step, new_root);
        //}
    }
  }
}

/* compute the new position of the particles in a node */
void remplirMyValues(node_t *n)
{
    if (!n)
        return;

    if (n->particle)
    {
        my_values[actual_n_particule*n_caracteristic_shared] = n->particle->x_pos;
        my_values[actual_n_particule*n_caracteristic_shared+1] = n->particle->y_pos;
        my_values[actual_n_particule*n_caracteristic_shared+2] = n->particle->x_vel;
        my_values[actual_n_particule*n_caracteristic_shared+3] = n->particle->y_vel;
        my_values[actual_n_particule*n_caracteristic_shared+4] = n->particle->x_force;
        my_values[actual_n_particule*n_caracteristic_shared+5] = n->particle->y_force;
        actual_n_particule+=1;
    }
    if (n->children)
    {
        int i;
        for (i = 0; i < 4; i++)
        {
            remplirMyValues(&n->children[i]);
        }
    }
}

void recvSendBuffer(node_t *n)
{
    if (!n)
        return;

    if (n->particle)
    {
        n->particle->x_pos =buffer_recv[n_caracteristic_shared*actual_n_particule];
        n->particle->y_pos =buffer_recv[n_caracteristic_shared*actual_n_particule+1];
        n->particle->x_vel =buffer_recv[n_caracteristic_shared*actual_n_particule+2];
        n->particle->y_vel =buffer_recv[n_caracteristic_shared*actual_n_particule+3];
        n->particle->x_force =buffer_recv[n_caracteristic_shared*actual_n_particule+4];
        n->particle->y_force =buffer_recv[n_caracteristic_shared*actual_n_particule+5];
        actual_n_particule+=1;
    }
    if (n->children)
    {
        int i;
        for (i = 0; i < 4; i++)
        {
            recvSendBuffer(&n->children[i]);

        }
    }
}

/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step) {
    //Un process par enfants
//    for (int i=0;i<4;i++) {
//        if (comm_rank == i) {
//            compute_force_in_node(&root->children[i]);
//        }
//    }
    compute_force_in_node(&root->children[comm_rank]);

    //changer les tableaux counts_recv  displacements_recv
    nParticulePerProcess = root->children[comm_rank].n_particles;
    free(my_values);
    my_values = malloc(nParticulePerProcess*n_caracteristic_shared*sizeof(double));
    for(int i = 0; i < comm_size; i++)
    {
        //RECV n particles informations from
        counts_recv[i] = root->children[i].n_particles * n_caracteristic_shared;
        if (i==0){
            displacements_recv[i] = 0;
        }
        else {
            displacements_recv[i] = displacements_recv[i-1]+ counts_recv[i-1];
        }
    }
    actual_n_particule=0;
    remplirMyValues(&root->children[comm_rank]);

//    // AFFICHAGE DES TABLEAU DE GESTION DE LA RECEPTION
//        MPI_Barrier(MPI_COMM_WORLD);
//        printf("\nComm_rank %d nParticules %d\n",comm_rank,nParticulePerProcess);
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

////  AFFICHAGE DE MY_VALUES
//    MPI_Barrier(MPI_COMM_WORLD);
//    for(int i = 0; i < nParticulePerProcess; i++)
//    {
//        printf("comm_rank %d : i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / x_force = %f / y_force = %f\n",
//               comm_rank,i,my_values[i*6+0],my_values[i*6+1],my_values[i*6+2],my_values[i*6+3],my_values[i*6+4],my_values[i*6+5]);
//    }
//    printf("\n");
//    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allgatherv(my_values,
                   nParticulePerProcess*n_caracteristic_shared,
                    MPI_DOUBLE,
                    buffer_recv,
                    counts_recv,
                    displacements_recv,
                    MPI_DOUBLE,
                    MPI_COMM_WORLD);
    //Deplacer donné de buffer send à particles
    actual_n_particule=0;
    for (int i=0;i<4;i++){
        recvSendBuffer(&root->children[i]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < nparticles; i++)
    {
        printf("comm_rank %d : i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / x_force = %f / y_force = %f\n",
               comm_rank,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].x_force,particles[i].y_force);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //CHACUN A paticles à jour
    MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank==0) {
      node_t *new_root = alloc_node();
      init_node(new_root, NULL, XMIN, XMAX, YMIN, YMAX);

      /* then move all particles and return statistics */
      move_particles_in_node(root, step, new_root);

      free_node(root);
      root = new_root;
  }
  //ENVOYER ROOT POUR CHACUN

}

void run_simulation()
{
  double t = 0.0, dt = 0.01;
    n_caracteristic_shared=6;
    counts_recv= malloc(comm_size*sizeof(int));
    displacements_recv = malloc(comm_size*sizeof(int));
    // Gestion des sources et de la destination des valeurs reçus
    buffer_recv= malloc(nparticles*n_caracteristic_shared*sizeof(double));


  while (t < T_FINAL && nparticles > 0)
  {
    /* Update time. */
    printf("comm_rank %d : t = %f / dt = %f / max_acc = %f / max_speed = %f\n",comm_rank,t,dt,max_acc_global,max_speed_global);
//      MPI_Barrier(MPI_COMM_WORLD);
//      for(int i = 0; i < nparticles; i++)
//      {
//          printf("comm_rank %d : t = %f / i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / x_force = %f / y_force = %f\n",
//                 comm_rank,t,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].x_force,particles[i].y_force);
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt);

    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */

    dt = 0.1 * max_speed / max_acc;

    /* Plot the movement of the particle */
#if DISPLAY
    node_t *n = root;
    clear_display();
    draw_node(n);
    flush_display();
#endif
  }
}

//PLUS TARD : PASSER UN QUART DE MAP si process puis allgather + root.first have children chacun de ceux la
/* create a quad-tree from an array of particles */
void insert_all_particles(int nparticles, particle_t *particles, node_t *root)
{
  int i;
  for (i = 0; i < nparticles; i++)
  {
    insert_particle(&particles[i], root);
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
// //OK TOUT LE MONDE EST A JOUR
//        MPI_Barrier(MPI_COMM_WORLD);
//    printf("\n Après Bcast\n");
//    for(int i = 0; i < nparticles; i++)
//    {
//        printf("comm_rank %d :i = %d / x_pos = %f / y_pos = %f / x_vel= %f / y_vel = %f / mass = %f\n",
//               comm_rank,i,particles[i].x_pos,particles[i].y_pos,particles[i].x_vel,particles[i].y_vel,particles[i].mass);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);


  insert_all_particles(nparticles, particles, root);

    if (comm_rank==0) {
        /* Initialize thread data structures */
#ifdef DISPLAY
        /* Open an X window to display the particles */
        simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

        struct timeval t1;
        gettimeofday(&t1, NULL);
    }
    MPI_Barrier(MPI_COMM_WORLD);

  /* Main thread starts simulation ... */
  run_simulation();


    MPI_Finalize();
    if (comm_rank==0) {
        struct timeval t2;
        gettimeofday(&t2, NULL);

        //double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

#ifdef DUMP_RESULT
        FILE *f_out = fopen("particles.log", "w");
        assert(f_out);
        print_particles(f_out, root);
        fclose(f_out);
#endif

        printf("-----------------------------\n");
        printf("nparticles: %d\n", nparticles);
        printf("T_FINAL: %f\n", T_FINAL);
        printf("-----------------------------\n");
        //printf("Simulation took %lf s to complete\n", duration);

#ifdef DISPLAY
        node_t *n = root;
        clear_display();
        draw_node(n);
        flush_display();

        printf("Hit return to close the window.");

        getchar();
        /* Close the X window used to display the particles */
        XCloseDisplay(theDisplay);
#endif
    }
  return 0;
}
