#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>    /* for getopt */
#include <stdint.h>
#include <pthread.h>

#include <zlib.h>
#include <seqtk/kseq.h>

//init kseq.h
KSEQ_INIT(gzFile, gzread);

// Take a string of nucleotides, and return an index into the kmer_table
inline uint16_t nucleotide_sequence_to_index(char* kmer, int length){
  uint16_t to_return = 0;

  int i;
  int positions_to_move;
  for (i=0; i<length; i++){
    positions_to_move = 2*(length-i-1);
    switch(kmer[i]){
      case 'A':
      case 'a':
        break; //i.e. add 0
      case 'C':
      case 'c':
        to_return += 1 << positions_to_move;
        break;
      case 'G':
      case 'g':
        to_return += 2 << positions_to_move;
        break;
      case 'T':
      case 't':
        to_return += 3 << positions_to_move;
        break;
      default:
        fprintf(stderr, "Detected unexcepted nucleotide: %c\n", kmer[i]);
        exit(1);
    }
  }
  if (to_return >= 1<<16){
    printf("Found kmer that was too big for the table! Fail.\n");
    exit(1);
  }

  return to_return;
}

typedef struct seqtab {
  char* s;
  int l;
  uint16_t* kmer_table;
} seq_and_table;

void* process_sequence(void* arg){
  seq_and_table* seq_table = (seq_and_table*) arg;
  if (seq_table->s < (char*)100){
    printf("Found abberant char*, failure imminent\n");
    exit(1);
  }
  int position = 0;
  for (position = 0; seq_table->l-position >= 8; position++){
    //find index converting nucleotide string to bit array
    uint16_t index = nucleotide_sequence_to_index(seq_table->s+position, 8);

    //stm_transact array lookup and increment
    __transaction_atomic {
      seq_table->kmer_table[index] += 1;
    }
  }
  return NULL;
}

int main(int argc, char *argv[]){
  char c;
  int num_threads = 1;

  while ((c = getopt (argc, argv, "ht:")) != -1){
    switch (c){
      exit(0);
      case 'h':
        printf("\n  Usage: %s [-t <num_threads>] <seq_file>\n\n",argv[0]);
        printf("  Count kmers, output some of them. The <seq_file> can be a FASTA or FASTQ file, gzipped or uncompressed.\n\n");
        printf("  Options:\n");
        printf("   -t THREADS  number of threads to use [default: %i]\n", num_threads);
        printf("   -h          show this help\n");
        printf("\n");
        exit(0);
      case 't':
        num_threads = atoi(optarg);
        break;
    }
  }
  printf("Using %i threads\n", num_threads);

  //calloc sufficient space for the kmer counts
  uint16_t* kmer_table = (uint16_t*) calloc(1<<16, sizeof(uint16_t));
  printf("biggest=%i\n", nucleotide_sequence_to_index("TTTTTTTT",8));
  printf("smallest=%i\n", nucleotide_sequence_to_index("AAAAAAAA",8));
  printf("middle=%i\n", nucleotide_sequence_to_index("AAAAAATA",8));

  //open input file
  //setup kseq reading
  gzFile fp;
  char* path = argv[optind];
  if (path == NULL){
    printf("Reading from stdin\n");
    fp = gzopen("/dev/stdin","r");
  } else {
    printf("Reading from %s\n", path);
    fp = gzopen(path, "r");
  }

  kseq_t *seq;
  int l;
  seq = kseq_init(fp);

  pthread_t* threads = malloc(num_threads*sizeof(pthread_t));
  seq_and_table* passing_structs = calloc(num_threads, sizeof(seq_and_table));
  //int pthread_create(pthread_t * pth, pthread_attr_t *att, void * (*function), void * arg);

  //foreach sequence
  int thread_number = 0;
  int sequence_number = 0;
  int pthread_return;
  seq_and_table* st;
  while ((l = kseq_read(seq)) >= 0) {
    pthread_return = pthread_join(threads[thread_number], NULL);
    if (sequence_number > thread_number && pthread_return != 0){
      printf("pthread_create returned with strange status %d, exiting", pthread_return);
      exit(1);
    }

    st = passing_structs+thread_number;
    if (sequence_number <= thread_number){
      st->kmer_table = kmer_table;
      st->s = (char *) malloc(seq->seq.l*sizeof(char));
      strncpy(st->s, seq->seq.s, seq->seq.l);
      st->l = seq->seq.l;
    } else {
      st->s = realloc(st->s, seq->seq.l);
      strncpy(st->s, seq->seq.s, seq->seq.l);
      st->l = seq->seq.l;
    }

    pthread_return = pthread_create(threads+thread_number, NULL, process_sequence, st);
    if (pthread_return != 0){
      printf("Unexpected pthread_create: %i\n",pthread_return);
      exit(1);
    }

    sequence_number += 1;
    thread_number += 1;
    if (thread_number >= num_threads)
      thread_number = 0;
  }

  //wait for threads to finish
  for (l=0; l<num_threads; l++){
    pthread_join(threads[l], NULL);
  }

  //cleanup
  kseq_destroy(seq);
  gzclose(fp);

  //print example kmer index (AAAA AAAA) which appears multiple times
  printf("Found AAAA AAAA %i times\n", kmer_table[0]);
  printf("Found AAAA AAAC %i times\n", kmer_table[1]);
  printf("Found AAAA AAAG %i times\n", kmer_table[2]);

  return 0;
}
