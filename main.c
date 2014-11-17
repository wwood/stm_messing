#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>    /* for getopt */
#include <stdint.h>

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
    //printf("Moving %c %i\n",kmer[i], positions_to_move);
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

  return to_return;
}

int main(int argc, char *argv[]){
  //calloc sufficient space for the kmer counts
  uint16_t* kmer_table = (uint16_t*) calloc(1<<16, sizeof(uint16_t));
  printf("biggest=%i\n", nucleotide_sequence_to_index("TTTTTTTT",8));
  printf("smallest=%i\n", nucleotide_sequence_to_index("AAAAAAAA",8));
  printf("middle=%i\n", nucleotide_sequence_to_index("AAAAAATA",8));

  //open input file
  //setup kseq reading
  gzFile fp;
  char* path = argv[1];
  if (path == NULL){
    fp = gzopen("/dev/stdin","r");
  } else {
    fp = gzopen(path, "r");
  }

  kseq_t *seq;
  int l;
  seq = kseq_init(fp);
  int position;

  //foreach sequence
  while ((l = kseq_read(seq)) >= 0) {
    //while there is another kmer in the read
    position = 0;
    for (position = 0; seq->seq.l-position >= 8; position++){
      //find index converting nucleotide string to bit array
      //printf("looking at %s\n", seq->seq.s+position);
      uint16_t index = nucleotide_sequence_to_index(seq->seq.s+position, 8);
      //printf("Found index %i\n",index);

      //stm_transact array lookup and increment
      kmer_table[index] += 1;
    }
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
