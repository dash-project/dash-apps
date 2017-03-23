#ifndef CONTIG_GENERATION_H
#define CONTIG_GENERATION_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

#ifndef MAXIMUM_CONTIG_SIZE
#define MAXIMUM_CONTIG_SIZE 100000
#endif

#ifndef KMER_LENGTH
#define KMER_LENGTH 19
#endif

#ifndef LOAD_FACTOR
#define LOAD_FACTOR 1
#endif

#ifndef LINE_SIZE
#define LINE_SIZE (KMER_LENGTH + 4)
#endif

#ifndef KMER_PACKED_LENGTH
#define KMER_PACKED_LENGTH 5
#endif

static double gettime(void)
{
  struct timeval tv;
  if (gettimeofday(&tv, NULL)) {
    perror("gettimeofday");
    abort();
  }
  return ((double)tv.tv_sec) + tv.tv_usec / 1000000.0;
}

/* K-mer data structure */
typedef struct kmer_t kmer_t;
struct kmer_t {
  char kmer[KMER_PACKED_LENGTH];
  char l_ext;
  char r_ext;
  kmer_t *next;
};

/* Start k-mer data structure */
typedef struct start_kmer_t start_kmer_t;
struct start_kmer_t {
  kmer_t *kmerPtr;
  start_kmer_t *next;
};

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t {
  kmer_t *head;  // Pointer to the first entry of that bucket
};

/* Hash table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
  int64_t size;     // Size of the hash table
  bucket_t *table;  // Entries of the hash table are pointers to buckets
};

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
  kmer_t *heap;
  int64_t posInHeap;
};

/* Returns the number of UFX kmers in a file */
int64_t getNumKmersInUFX(const char *filename)
{
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "Could not open %s for reading!\n", filename);
    return -1;
  }
  char firstLine[LINE_SIZE + 1];
  firstLine[LINE_SIZE] = '\0';
  if (fread(firstLine, sizeof(char), LINE_SIZE, f) != LINE_SIZE) {
    fprintf(stderr, "Could not read %d bytes!\n", LINE_SIZE);
    return -2;
  }
  // check structure and size of kmer is correct!
  if (firstLine[LINE_SIZE] != '\0') {
    fprintf(stderr,
            "UFX text file is an unexpected line length for kmer length %d\n",
            KMER_LENGTH);
    return -3;
  }
  if (firstLine[KMER_LENGTH] != ' ' && firstLine[KMER_LENGTH] != '\t') {
    fprintf(stderr, "Unexpected format for firstLine '%s'\n", firstLine);
    return -4;
  }

  struct stat buf;
  int fd = fileno(f);
  if (fstat(fd, &buf) != 0) {
    fprintf(stderr, "Could not stat %s\n", filename);
    return -5;
  }
  int64_t totalSize = buf.st_size;
  if (totalSize % LINE_SIZE != 0) {
    fprintf(stderr,
            "UFX file is not a multiple of %d bytes for kmer length %d\n",
            LINE_SIZE, KMER_LENGTH);
    return -6;
  }
  fclose(f);
  int64_t numKmers = totalSize / LINE_SIZE;
  printf("Detected %j kmers in text UFX file: %s\n", numKmers, filename);
  return numKmers;
}

#endif  // CONTIG_GENERATION_H
