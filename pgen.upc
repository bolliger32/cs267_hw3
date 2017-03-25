#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_collective.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


shared int nKmers,total_chars_to_read,*totalStart;
shared int64_t cur_start_pos[THREADS],start_pos[THREADS];
int main(int argc, char *argv[]){
    
	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
    /* Read the input file name */
    char *input_UFX_name = argv[1];
    if (MYTHREAD == 0) {
        nKmers = getNumKmersInUFX(input_UFX_name);
    }
     init_LookupTable();
	///////////////////////////////////////////
	upc_barrier;
	inputTime += gettime();
    
    
    
	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
    /* Create a hash table */
    
//    shared_hash_table_t* hashtable = create_shared_hash_table(nKmers,&memory_heap, nKmers);
    int64_t n_buckets = (nKmers+LOAD_FACTOR-1)*LOAD_FACTOR;
    shared uint64_t* hashtable = (shared uint64_t*) upc_all_alloc(n_buckets, sizeof(uint64_t));
    upc_forall(int i = 0; i<n_buckets; i++; hashtable[i]) {
        hashtable[i]=-1;
    }
    
    /* Create collisions table for indices when there is a collision */
    shared uint64_t* collisions = (shared uint64_t*) upc_all_alloc(nKmers,sizeof(uint64_t));
    upc_forall(int i = 0; i<nKmers; i++; i) {
        collisions[i]= (uint64_t) -1;
    }
    
    /* Create memory heap */
    int mh_block_size = (nKmers+(THREADS-1))/THREADS;
    shared kmer_t* memory_heap = (shared kmer_t*) upc_all_alloc(THREADS, sizeof(kmer_t)*mh_block_size);
    
    // Create start kmers list
    int64_t* startKmersList_local = (int64_t*) malloc(nKmers/THREADS*sizeof(int64_t));

    int64_t cur_start_pos_local=0;

    
   /* Read the kmers from the input file and store them in the  memory_heap */
    total_chars_to_read = nKmers * LINE_SIZE;

    int start_ix = mh_block_size*MYTHREAD;
    int stop_ix = mh_block_size*(MYTHREAD+1);
    if (MYTHREAD==(THREADS-1)) stop_ix = nKmers;
    int lines_to_read = stop_ix - start_ix;
    FILE *inputFile = fopen(input_UFX_name,"r");
    fseek(inputFile, sizeof(unsigned char)*start_ix*LINE_SIZE,SEEK_SET);
    
    unsigned char wb[LINE_SIZE];
    char packedKmer[KMER_PACKED_LENGTH];
    for (int mh_ix = start_ix; mh_ix < stop_ix; mh_ix++) {
        fread(wb, sizeof(unsigned char), LINE_SIZE, inputFile);
        char left_ext = (char) wb[KMER_LENGTH+1];
        char right_ext = (char) wb[KMER_LENGTH+2];
        /* Pack a k-mer sequence appropriately */
        packSequence(wb, (unsigned char*) packedKmer, KMER_LENGTH);
        int64_t hashval = hashkmer(n_buckets, (char*) packedKmer);

        /* place into heap */
        upc_memput(memory_heap[mh_ix].kmer, packedKmer,KMER_PACKED_LENGTH*sizeof(char));
        
        memory_heap[mh_ix].r_ext = right_ext;
        memory_heap[mh_ix].l_ext = left_ext;
        
        /* place into hashtable or collisions table */
        int res = bupc_atomicI_cswap_strict(hashtable+hashval,-1,mh_ix);
        while (res != -1) {
            res = bupc_atomicI_cswap_strict(collisions + res,-1,mh_ix);
        }
        
        // Add to start list
        if (left_ext == 'F') {
            startKmersList_local[cur_start_pos_local]=mh_ix;
            cur_start_pos_local++;
            
        }
    }
    fclose(inputFile);
    cur_start_pos[MYTHREAD] = cur_start_pos_local;
    bupc_atomicI_fetchadd_relaxed(totalStart,cur_start_pos_local);
    upc_all_prefix_reduceI(start_pos,cur_start_pos, UPC_ADD, THREADS, 1, NULL, UPC_IN_NOSYNC | UPC_OUT_ALLSYNC);
//    upc_all_reduceI(totalStart,cur_start_pos,UPC_ADD,THREADS,1, NULL, UPC_IN_NOSYNC | UPC_OUT_ALLSYNC);
    upc_barrier;
    shared [] int64_t* startKmersList = (shared [] int64_t*) upc_all_alloc(1, *totalStart*sizeof(int64_t));

    upc_memput(&startKmersList[start_pos[MYTHREAD]-cur_start_pos[MYTHREAD]],startKmersList_local,cur_start_pos_local*sizeof(int64_t));
    free(startKmersList_local);
	constrTime += gettime();
    
	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
    char buff[10];
    int r = sprintf(buff,"%d",MYTHREAD);
    char* buf = (char*) malloc((r + 9) * sizeof(char));
   sprintf(buf, "pgen_%d.out", MYTHREAD); // puts string into buffer
   FILE* outputFile = fopen(buf, "w");

   /* Pick start nodes from the startKmersList */
    
    char unpackedKmer[KMER_LENGTH+1];
    unpackedKmer[KMER_LENGTH] = '\0';
    char pKmer[KMER_PACKED_LENGTH];
    int64_t contigID = 0, totBases = 0;

    upc_forall(int i = 0; i < *totalStart; i++; i) {
        uint64_t cur_start_ix = startKmersList[i];
        unpackSequenceShared((shared unsigned char*) memory_heap[cur_start_ix].kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);

      /* Initialize current contig with the seed content */
      char cur_contig[MAXIMUM_CONTIG_SIZE];
      memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
      int posInContig = KMER_LENGTH;
      char right_ext = memory_heap[cur_start_ix].r_ext;

      /* Keep adding bases while not finding a terminal node */
      while (right_ext != 'F') {
         cur_contig[posInContig] = right_ext;
         posInContig++;

         /* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
          int64_t cur_kmer_ptr = lookup_kmer_shared(hashtable, (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH], n_buckets,memory_heap, collisions);
         right_ext = memory_heap[cur_kmer_ptr].r_ext;
      }
    
    
    
    
      /* Print the contig since we have found the corresponding terminal node */
      cur_contig[posInContig] = '\0';
      fprintf(outputFile,"%s\n", cur_contig);
      contigID++;
      totBases += strlen(cur_contig);
    }
   fclose(outputFile);
    free(buf);

	////////////////////////////////////////////////////////////
	traversalTime += gettime();
    upc_barrier;
    upc_free(hashtable);
//    upc_free(startKmersList); // crashes for some reason
    upc_free(collisions);
    upc_free(memory_heap);
	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
