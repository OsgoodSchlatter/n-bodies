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
#include <stdbool.h>

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



int n_caracteristic_shared;
double** pToShare;
double** pToRecv;
int indexPToShare[16];
int indexPToRecv[16];
int carac_to_share=7; //x_pos, y_pos;		/* position of the particle */ x_vel, y_vel;		/* velocity of the particle */ x_force, y_force;	/* gravitational forces that apply against this particle */ mass;
int indexCommSend;
int indexCommRecv;

void init()
{
    init_alloc(8 * nparticles);
    root = malloc(sizeof(node_t));
    init_node(root, NULL, XMIN, XMAX, YMIN, YMAX);

    pToShare = malloc(sizeof(double)*16);
    for(int i = 0; i < 16; i++)
    {
        pToShare[i] = malloc(nparticles*sizeof(double)*carac_to_share);
    }

    pToRecv = malloc(sizeof(double)*16);
    for(int i = 0; i < 16; i++)
    {
        pToRecv[i] = malloc(nparticles*sizeof(double)*carac_to_share);
    }


#ifdef DISPLAY
  Display *theDisplay; /* These three variables are required to open the */
  GC theGC;            /* particle plotting window.  They are externally */
  Window theMain;      /* declared in ui.h but are also required here.   */
#endif
}

int isInWhichNode(particle_t *p){
    if (p->x_pos < XMIN ||
        p->x_pos > XMAX ||
        p->y_pos < YMIN ||
        p->y_pos > YMAX)
    {
        //pas dans la grille
        return -1;
    }
    else
    {
        if ( p->x_pos > 0){
            if (p->y_pos > 0){
                return 3;
            } else{
                return 1;
            }
        } else{
            if (p->y_pos > 0){
                return 2;
            } else{
                return 0;
            }
        }

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
  if (!n) {
      return;
  }

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
  // VERIFIER SI AU BONNE ENDROIT
  // ENVOYER OU PAS AU BON PROCESS MPI
  // CE PROCESS L'insrt dans son root
  p->node = NULL;
  int indexNode = isInWhichNode(p);
  if (indexNode<0)
  {
    nparticles--;
  }
  else if (indexNode == comm_rank){
      insert_particle(p, new_root);
  }
  else
  {
      //printf("[%d/%d] p->x_pos %f p->y_pos %f\n",comm_rank,comm_size,p->x_pos,p->y_pos );
      indexCommSend = comm_rank*4+indexNode;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+0]=p->x_pos;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+1]=p->y_pos;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+2]=p->x_vel;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+3]=p->y_vel;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+4]=p->x_force;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+5]=p->y_force;
      pToShare[indexCommSend][indexPToShare[indexCommSend]*carac_to_share+6]=p->mass;

      indexPToShare[indexCommSend]+=1;
      //printf("++ [%d/%d] indexPToShare[%d]=%d\n",comm_rank,comm_size,indexNode,indexPToShare[indexNode]);
  }
}

/* compute the new position of the particles in a node */
void move_particles_in_node(node_t *n, double step, node_t *new_root)
{
  if (!n)
  {
    return;
  }

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
      // if (comm_rank==i){
      move_particles_in_node(&n->children[i], step, new_root);
      //}
    }
  }
}

void recvRecvBuffer2(){
  // for(int rank=0;rank<4;rank++){
  //   if (comm_rank!=rank){
  //     int indexCommSend=comm_rank*4+rank;
  //     MPI_Sendrecv_replace(&pToShare[indexCommSend], indexPToShare[indexCommSend]*carac_to_share, MPI_DOUBLE, rank, 0, comm_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//comm_rank->rank
  //     int indexCommRecv=rank*4+comm_rank;
  //     MPI_Sendrecv_replace(&pToShare[indexCommRecv], indexPToShare[indexCommRecv]*carac_to_share, MPI_DOUBLE, comm_rank, 0, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//rank->comm_rank
  //   }
  // }

  if (comm_rank==0){
    //send
    MPI_Sendrecv_replace(&pToShare[1], indexPToShare[1]*carac_to_share, MPI_DOUBLE, 1, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->1
    MPI_Sendrecv_replace(&pToShare[2], indexPToShare[2]*carac_to_share, MPI_DOUBLE, 2, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->2
    MPI_Sendrecv_replace(&pToShare[3], indexPToShare[3]*carac_to_share, MPI_DOUBLE, 3, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->3

    //recv
    MPI_Sendrecv_replace(&pToShare[4], indexPToShare[4]*carac_to_share, MPI_DOUBLE, 0, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->0
    MPI_Sendrecv_replace(&pToShare[8], indexPToShare[8]*carac_to_share, MPI_DOUBLE, 0, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->0
    MPI_Sendrecv_replace(&pToShare[12], indexPToShare[12]*carac_to_share, MPI_DOUBLE, 0, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->0
  }
  if (comm_rank==1){
    MPI_Sendrecv_replace(&pToShare[1], indexPToShare[1]*carac_to_share, MPI_DOUBLE, 1, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->1
    //send
    MPI_Sendrecv_replace(&pToShare[4], indexPToShare[4]*carac_to_share, MPI_DOUBLE, 0, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->0
    MPI_Sendrecv_replace(&pToShare[6], indexPToShare[6]*carac_to_share, MPI_DOUBLE, 2, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->2
    MPI_Sendrecv_replace(&pToShare[7], indexPToShare[7]*carac_to_share, MPI_DOUBLE, 3, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->3

    MPI_Sendrecv_replace(&pToShare[9], indexPToShare[9]*carac_to_share, MPI_DOUBLE, 1, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->1
    MPI_Sendrecv_replace(&pToShare[13], indexPToShare[13]*carac_to_share, MPI_DOUBLE, 1, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->1
  }
  if (comm_rank==2){
    MPI_Sendrecv_replace(&pToShare[2], indexPToShare[2]*carac_to_share, MPI_DOUBLE, 2, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->2
    MPI_Sendrecv_replace(&pToShare[6], indexPToShare[6]*carac_to_share, MPI_DOUBLE, 2, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->2

    //send
    MPI_Sendrecv_replace(&pToShare[8], indexPToShare[8]*carac_to_share, MPI_DOUBLE, 0, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->0
    MPI_Sendrecv_replace(&pToShare[9], indexPToShare[9]*carac_to_share, MPI_DOUBLE, 1, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->1
    MPI_Sendrecv_replace(&pToShare[11], indexPToShare[11]*carac_to_share, MPI_DOUBLE, 3, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->3

    MPI_Sendrecv_replace(&pToShare[14], indexPToShare[14]*carac_to_share, MPI_DOUBLE, 2, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->2
  }
  if (comm_rank==3){
    MPI_Sendrecv_replace(&pToShare[3], indexPToShare[3]*carac_to_share, MPI_DOUBLE, 3, 0, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//0->3
    MPI_Sendrecv_replace(&pToShare[7], indexPToShare[7]*carac_to_share, MPI_DOUBLE, 3, 0, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//1->3
    MPI_Sendrecv_replace(&pToShare[11], indexPToShare[11]*carac_to_share, MPI_DOUBLE, 3, 0, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//2->3

    //send
    MPI_Sendrecv_replace(&pToShare[12], indexPToShare[12]*carac_to_share, MPI_DOUBLE, 0, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->0
    MPI_Sendrecv_replace(&pToShare[13], indexPToShare[13]*carac_to_share, MPI_DOUBLE, 1, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->1
    MPI_Sendrecv_replace(&pToShare[14], indexPToShare[14]*carac_to_share, MPI_DOUBLE, 2, 0, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);//3->2
  }
}
void recvRecvBuffer(){
  MPI_Allgather(&indexPToShare[comm_rank*4], 4, MPI_INT, &indexPToRecv, 16, MPI_INT, MPI_COMM_WORLD);

  for(int rank=0;rank<4;rank++){
    if (comm_rank!=rank){
      indexCommSend=comm_rank*4+rank;//comm_rank->rank
      printf("%d SEND %d DATA TO %d : %d\n",comm_rank, indexPToRecv[indexCommSend],rank,indexCommSend);
      MPI_Send(&pToShare[indexCommSend], indexPToRecv[indexCommSend]*carac_to_share, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
    }
  }
  MPI_Request requests[3];
  int i=0;
  for(int rank=0;rank<4;rank++){
    if (comm_rank!=rank){
      indexCommRecv=rank*4+comm_rank;//rank->comm_rank
      printf("%d RECV %d DATA FROM %d : %d\n",comm_rank, indexPToRecv[indexCommRecv],rank,indexCommRecv );
      MPI_Irecv(&pToRecv[indexCommRecv], indexPToRecv[indexCommRecv]*carac_to_share, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, &requests[i]);
      i++;
    }
  }
  MPI_Waitall(3, requests, MPI_STATUSES_IGNORE);
}

void insertNewParticule(node_t *new_root){
    //printf("[%d/%d] insertNewParticule IN\n",comm_rank,comm_size);
    particle_t *p;
    p = malloc(sizeof(particle_t));
    p->node=NULL;
    for (int rank = 0; rank<4 ;rank++){
        //printf("[%d/%d] indexBuffer[%d] = %d\n",comm_rank,comm_size,i,indexBuffer[i]);
        if (comm_rank!=rank){
          indexCommRecv=rank*4+comm_rank;//recu de rank
            for (int j=0;j<indexPToShare[indexCommRecv];j++){

                p->x_pos=pToRecv[indexCommRecv][carac_to_share*j+0];
                p->y_pos=pToRecv[indexCommRecv][carac_to_share*j+1];
                p->x_vel=pToRecv[indexCommRecv][carac_to_share*j+2];
                p->y_vel=pToRecv[indexCommRecv][carac_to_share*j+3];
                p->x_force=pToRecv[indexCommRecv][carac_to_share*j+4];
                p->y_force=pToRecv[indexCommRecv][carac_to_share*j+5];
                p->mass=pToRecv[indexCommRecv][carac_to_share*j+6];
                //printf("[%d/%d] new_root->x_max = %f new_root->x_min = %f , p->x_pos %f p->y_pos %f\n",comm_rank,comm_size,new_root->x_max,new_root->x_min,p->x_pos,p->y_pos);
                insert_particle(p, new_root);
              //  printf("[%d/%d] sender = %d, current_p = %d, p->x_pos = %f\n",comm_rank,comm_size,i,j,p->x_pos);

            }
        }
    }
}


/*
  Move particles one time step.

  Update positions, velocity, and acceleration.
  Return local computations.
*/
void all_move_particles(double step)
{

  MPI_Barrier(MPI_COMM_WORLD);

    // On pourrait faire passer que x_force et y force dans cette partie
  compute_force_in_node(root);

  //printf("[%d/%d] compute_force_in_node pass\n",comm_rank,comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

    node_t *new_root = alloc_node();
    init_node(new_root, NULL,root->x_min, root->x_max, root->y_min, root->y_max);

    //printf("[%d/%d] init_node pass\n",comm_rank,comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    /* then move all particles and return statistics */

    move_particles_in_node(root, step, new_root);

    //printf("[%d/%d] move_particule_in_node pass\n",comm_rank,comm_size);

    MPI_Barrier(MPI_COMM_WORLD);


    //S'envoyer les particules non intégrées + inserer dans les noeuds
    //S'envoyer les tailles des buffers indexPToShare

    recvRecvBuffer();

    MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d/%d] recvRecvBuffer pass\n",comm_rank,comm_size);
    insertNewParticule(new_root);
    //dechargerRecvBuffer() //use tableau[3] de recvBuffer
    //insertSharedParticules()

    //printf("[%d/%d] insertNewParticule pass\n",comm_rank,comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

    free_node(root);
    root = new_root;
    //printf("[%d/%d] root=new_root pass\n",comm_rank,comm_size);
    MPI_Barrier(MPI_COMM_WORLD);

}

void run_simulation()
{
  double t = 0.0, dt = 0.01;

  while (t < T_FINAL && nparticles > 0)
  {
    /* Update time. */
    //      MPI_Barrier(MPI_COMM_WORLD);
    printf("[%d/%d] t = %f\n",comm_rank,comm_size,t);
    for (int i=0;i<4;i++){
        indexPToShare[i]=0;
    }
    t += dt;
    /* Move particles with the current and compute rms velocity. */
    all_move_particles(dt);
    printf("[%d/%d] all_move_particule pass\n",comm_rank,comm_size);
    //Gérer max_acc max_speed
    //MPI_Bcast(&max_acc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //MPI_Bcast(&max_speed, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Adjust dt based on maximum speed and acceleration--this
       simple rule tries to insure that no velocity will change
       by more than 10% */
    dt = 0.1 ;//* max_speed / max_acc;

    /* Plot the movement of the particle */
    if (comm_rank == 0)
    {
#if DISPLAY
      node_t *n = root;
      clear_display();
      draw_node(n);
      flush_display();
#endif
    }
  }
}

// PLUS TARD : PASSER UN QUART DE MAP si process puis allgather + root.first have children chacun de ceux la
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
  double *bcast_buff = malloc(sizeof(double) * 5 * nparticles); // x_pos,y_pos,x_vel,y_vel,mass
  if (comm_rank == 0)
  {
    all_init_particles(nparticles, particles);
    for (int i = 0; i < nparticles; i++)
    {
      bcast_buff[i * 5 + 0] = particles[i].x_pos;
      bcast_buff[i * 5 + 1] = particles[i].y_pos;
      bcast_buff[i * 5 + 2] = particles[i].x_vel;
      bcast_buff[i * 5 + 3] = particles[i].y_vel;
      bcast_buff[i * 5 + 4] = particles[i].mass;
    }
  }
  MPI_Bcast(bcast_buff, 5 * nparticles, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (comm_rank != 0)
  {
    for (int i = 0; i < nparticles; i++)
    {
      particles[i].x_pos = bcast_buff[i * 5 + 0];
      particles[i].y_pos = bcast_buff[i * 5 + 1];
      particles[i].x_vel = bcast_buff[i * 5 + 2];
      particles[i].y_vel = bcast_buff[i * 5 + 3];
      particles[i].mass = bcast_buff[i * 5 + 4];
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
  *root=root->children[comm_rank];

  printf("[%d/%d] x_max %f x_min %f y_max %f y_min %f\n",comm_rank,comm_size,root->x_max,root->x_min,root->y_max,root->y_min);

  struct timeval t1, t2;
  if (comm_rank == 0)
  {
    /* Initialize thread data structures */
#ifdef DISPLAY
    /* Open an X window to display the particles */
    simple_init(100, 100, DISPLAY_SIZE, DISPLAY_SIZE);
#endif

    gettimeofday(&t1, NULL);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* Main thread starts simulation ... */
  run_simulation();

  MPI_Finalize();

  if (comm_rank == 0)
  {
    gettimeofday(&t2, NULL);

    double duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

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
    printf("Simulation took %lf s to complete\n", duration);

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
