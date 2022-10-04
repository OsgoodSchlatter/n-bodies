#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "nbody_alloc.h"

static void* mem_new_buffer(struct memory_t *mem, int nb_blocks) {
  Bloc *ptr;
  ptr = calloc(nb_blocks, mem->block_size);
  assert(ptr != 0);
  /* Decoupage du bloc retourne par malloc en nbMaxBloc constituant une liste*/
  /* chainee                                                                 */
  char* block_ptr = (char*)ptr;
  Bloc*p = (Bloc*)block_ptr;
  for (int i=0 ; i<(nb_blocks-1) ; i++) {
    p = (Bloc*) block_ptr;
    block_ptr += mem->block_size;
    p->suivant = (Bloc*) block_ptr;
  }
  p->suivant = NULL;
  return ptr;
}

/**************************************************************************/
/* Routine d'initialisation de l'allocateur : constitue une liste chainee */
/* de nbMaxAlloc blocs de taille TAILLEBLOC, le premier bloc etant pointe */
/* par debutListe et le dernier pointant sur NULL                         */
/**************************************************************************/
void mem_init(struct memory_t *mem, size_t block_size, int nb_blocks)
{
  /* Memorisation du debut de la liste chainee */
  mem->block_size = block_size;
  mem->debutListe = mem_new_buffer(mem, nb_blocks);
  mem->nb_allocated = nb_blocks;
  mem->nb_free = nb_blocks;
}

/**************************************************************************/
/* Fonction renvoyant un pointeur sur une zone memoire                    */
/* NB : on ne peut pas preciser la taille, puisque la taille est          */
/*      predefinie                                                        */
/**************************************************************************/
void *mem_alloc(struct memory_t*mem)
{
  Bloc *ptr;
  if(mem->debutListe == NULL) {
    /* we need to allocate more memory blocks */
    mem->debutListe = mem_new_buffer(mem, mem->nb_allocated);
    mem->nb_free += mem->nb_allocated;
    mem->nb_allocated *=2;
  }

  ptr = mem->debutListe;
  mem->debutListe = mem->debutListe->suivant;
  memset(ptr, 0, mem->block_size);

  mem->nb_free--;
  return (void *)ptr;
}

/**************************************************************************/
/* Fonction liberant la zone memoire pointee par ptr                      */
/**************************************************************************/
void mem_free(struct memory_t* mem, void *ptr)
{
  Bloc *pBloc = (Bloc*)ptr;
  pBloc->suivant = mem->debutListe;
  mem->debutListe = pBloc;
  mem->nb_free++;
}
