#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "contig_generation.h"

shared_hash_table_t* create_shared_hash_table(int64_t nEntries, shared_memory_heap_t *memory_heap, int nKmers)
{
   shared_hash_table_t *result;
   int64_t n_buckets = nEntries * LOAD_FACTOR;

   result = (shared_hash_table_t*) malloc(sizeof(shared_hash_table_t));
   result->size = n_buckets;
   result->table = (shared shared_bucket_t*) upc_all_alloc(n_buckets , sizeof(shared_bucket_t));
   
   if (result->table == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %llu bytes\n", n_buckets, sizeof(shared_bucket_t));
      exit(1);
   }
   
   memory_heap->heap = (shared kmer_t *) upc_all_alloc(THREADS, sizeof(shared kmer_t)*ceil(nKmers/THREADS));
   if (memory_heap->heap == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      exit(1);
   }
   memory_heap->posInHeap = MYTHREAD * ceil(nKmers/THREADS);
   
   return result;
}

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
int64_t lookup_kmer_shared(shared uint64_t *hashtable, const unsigned char *kmer, int64_t size, shared kmer_t* memory_heap, shared uint64_t* collisions)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(size, (char*) packedKmer);
   uint64_t result;
   
   result = (int) hashtable[hashval];
   char comp_kmer[KMER_PACKED_LENGTH];
    upc_memget(comp_kmer, memory_heap[result].kmer,KMER_PACKED_LENGTH);

   if ( memcmp(packedKmer, comp_kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
       return result;
   }
   else {
        while ( memcmp(packedKmer, comp_kmer, KMER_PACKED_LENGTH * sizeof(char)) != 0 ) {
            result = (int) collisions[result];
            upc_memget(comp_kmer, memory_heap[result].kmer,KMER_PACKED_LENGTH);
            
        }
       return result;
    }
}

#endif // KMER_HASH_H
