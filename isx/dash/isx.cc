/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions 
are met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above 
      copyright notice, this list of conditions and the following 
      disclaimer in the documentation and/or other materials provided 
      with the distribution.
    * Neither the name of Intel Corporation nor the names of its 
      contributors may be used to endorse or promote products 
      derived from this software without specific prior written 
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
*/

#include <mpi.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <unistd.h> // sleep()
#include <sys/stat.h>
#include <stdint.h>

#include <libdash.h>

#include "params.h"
#include "isx.h"
#include "timer.h"
#include "pcg_basic.h"

#define ROOT_PE 0

uint64_t NUM_PES; // Number of parallel workers
uint64_t TOTAL_KEYS; // Total number of keys across all PEs
uint64_t NUM_KEYS_PER_PE; // Number of keys generated on each PE
uint64_t NUM_BUCKETS; // The number of buckets in the bucket sort
uint64_t BUCKET_WIDTH; // The size of each bucket
uint64_t MAX_KEY_VAL; // The maximum possible generated key value

int my_rank;
int comm_size;


#ifdef PERMUTE
int * permute_array;
#endif

int main(int argc,  char ** argv)
{
  //start_pes(0);
  dash::init(&argc, &argv);
  my_rank = dash::myid();
  comm_size = dash::size();

  char * log_file = parse_params(argc, argv);

  int err =  bucket_sort();

  log_times(log_file);

  dash::finalize();
  return err;
}


// Parses all of the command line input and definitions in params.h
// to set all necessary runtime values and options
static char * parse_params(const int argc, char ** argv)
{
  if(argc != 3)
  {
    if( my_rank == 0){
      printf("Usage:  \n");
      printf("  ./%s <total num keys(strong) | keys per pe(weak)> <log_file>\n",argv[0]);
    }
    exit(1);
  }

  NUM_PES = (uint64_t) comm_size;
  MAX_KEY_VAL = DEFAULT_MAX_KEY;
  NUM_BUCKETS = NUM_PES;
  BUCKET_WIDTH = (uint64_t) ceil((double)MAX_KEY_VAL/NUM_BUCKETS);
  char * log_file = argv[2];
  char scaling_msg[64];

  switch(SCALING_OPTION){
    case STRONG:
      {
        TOTAL_KEYS = (uint64_t) atoi(argv[1]);
        NUM_KEYS_PER_PE = (uint64_t) ceil((double)TOTAL_KEYS/NUM_PES);
        sprintf(scaling_msg,"STRONG");
        break;
      }

    case WEAK:
      {
        NUM_KEYS_PER_PE = (uint64_t) (atoi(argv[1]));
        sprintf(scaling_msg,"WEAK");
        break;
      }

    case WEAK_ISOBUCKET:
      {
        NUM_KEYS_PER_PE = (uint64_t) (atoi(argv[1]));
        BUCKET_WIDTH = ISO_BUCKET_WIDTH; 
        MAX_KEY_VAL = (uint64_t) (NUM_PES * BUCKET_WIDTH);
        sprintf(scaling_msg,"WEAK_ISOBUCKET");
        break;
      }

    default:
      {
        if(my_rank == 0){
          printf("Invalid scaling option! See params.h to define the scaling option.\n");
        }
        exit(1);
        break;
      }
  }

  assert(MAX_KEY_VAL > 0);
  assert(NUM_KEYS_PER_PE > 0);
  assert(NUM_PES > 0);
  assert(MAX_KEY_VAL > NUM_PES);
  assert(NUM_BUCKETS > 0);
  assert(BUCKET_WIDTH > 0);

  if(my_rank == 0){
    printf("ISx MPI 2 sided v%1d.%1d\n",MAJOR_VERSION_NUMBER,MINOR_VERSION_NUMBER);
#ifdef PERMUTE
    printf("Random Permute Used in ATA.\n");
#endif
    printf("  Number of Keys per PE: %" PRIu64 "\n", NUM_KEYS_PER_PE);
    printf("  Max Key Value: %" PRIu64 "\n", MAX_KEY_VAL);
    printf("  Bucket Width: %" PRIu64 "\n", BUCKET_WIDTH);
    printf("  Number of Iterations: %u\n", NUM_ITERATIONS);
    printf("  Number of PEs: %" PRIu64 "\n", NUM_PES);
    printf("  %s Scaling!\n",scaling_msg);
    }

  return log_file;
}


/*
 * The primary compute function for the bucket sort
 * Executes the sum of NUM_ITERATIONS + BURN_IN iterations, as defined in params.h
 * Only iterations after the BURN_IN iterations are timed
 * Only the final iteration calls the verification function
 */
static int bucket_sort(void)
{
  int err = 0; 

  init_timers(NUM_ITERATIONS);

#ifdef PERMUTE
  create_permutation_array();
#endif

  for(unsigned int i = 0; i < (NUM_ITERATIONS + BURN_IN); ++i)
  {

    // Reset timers after burn in 
    if(i == BURN_IN){ init_timers(NUM_ITERATIONS); }

    dash::barrier();

    timer_start(&timers[TIMER_TOTAL]);

    if (my_rank == 0) {
      std::cout << "Allocating bucket_sizes(" << NUM_PES
                << ", " << NUM_BUCKETS << "); total = "
                << NUM_PES * NUM_BUCKETS << std::endl;
    }
    // two-dimensional array holding the local bucket sizes for each process
    timer_start(&timers[TIMER_BSIZE_ALLOC]);
    dash::NArray<int, 2> bucket_sizes(NUM_PES, NUM_BUCKETS, dash::BLOCKED, dash::NONE);
    timer_stop(&timers[TIMER_BSIZE_ALLOC]);

    uninitialized_vector<KEY_TYPE> my_keys = make_input();

    count_local_bucket_sizes(my_keys, bucket_sizes);

    // TODO: We need dash::accumulate here!
    int max_bucket_size = *dash::max_element(
                             bucket_sizes.begin(),
                             bucket_sizes.end());

    if (my_rank == 0) {
      std::cout << "Allocating bucket_sizes(" << NUM_PES
                << ", " << NUM_BUCKETS
                << ", " << max_bucket_size << "); total = "
                << NUM_PES * NUM_BUCKETS * max_bucket_size << std::endl;
    }
    // 3-dimensional array holding the bucketed elements per unit
    timer_start(&timers[TIMER_BUCK_ALLOC]);
    dash::NArray<int, 3> bucketed_keys(NUM_PES, NUM_BUCKETS, max_bucket_size,
                                       dash::BLOCKED, dash::NONE, dash::NONE);
    timer_stop(&timers[TIMER_BUCK_ALLOC]);
    if (my_rank == 0) {
      std::cout << "Done allocating buckets" << std::endl;
    }

    bucketize_local_keys(my_keys, bucketed_keys, bucket_sizes);
    // release the allocated memory
    my_keys.free();

//    if (my_rank == 0) {
//      for (int i = 0; i < comm_size; ++i) {
//        for (int j = 0; j < comm_size; ++j) {
//          for (int k = 0; k < max_bucket_size; ++k) {
//            std::cout << "bucketed_keys[" << i << "][" << j << "][" << k << "] = "
//                      << (int)bucketed_keys[i][j][k] << std::endl;
//          }
//        }
//      }
//    }

    long long int  my_bucket_size;
    uninitialized_vector<KEY_TYPE> my_bucket_keys = exchange_keys(bucketed_keys,
                                                         bucket_sizes,
                                                         my_bucket_size);


    std::vector<CNT_TYPE> my_local_key_counts = count_local_keys(my_bucket_keys, my_bucket_size);

    dash::barrier();

    timer_stop(&timers[TIMER_TOTAL]);

    // verify the burn-in iteration
    if(i == 0) {
      err = verify_results(my_local_key_counts, my_bucket_keys, my_bucket_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
  return err;
}


/*
 * Generates uniformly random keys [0, MAX_KEY_VAL] on each rank using the time and rank
 * number as a seed
 */
static uninitialized_vector<KEY_TYPE> make_input(void)
{
  timer_start(&timers[TIMER_INPUT]);

  // TODO: this is a performance sink hole as (potentially GBs)
  //       are initialized without need.
  uninitialized_vector<KEY_TYPE> my_keys(NUM_KEYS_PER_PE);
  KEY_TYPE *__restrict mk = my_keys.data();

  pcg32_random_t rng = seed_my_rank();

  for(unsigned int i = 0; i < NUM_KEYS_PER_PE; ++i) {
    mk[i] = pcg32_boundedrand_r(&rng, MAX_KEY_VAL);
  }

  timer_stop(&timers[TIMER_INPUT]);

#ifdef DEBUG
  
  char msg[1024];
  sprintf(msg,"Rank %d: Initial Keys: ", my_rank);
  for(int i = 0; i < NUM_KEYS_PER_PE; ++i){
    if(i < PRINT_MAX)
    sprintf(msg + strlen(msg),"%d ", my_keys[i]);
  }
  sprintf(msg + strlen(msg),"\n");
  printf("%s",msg);
  fflush(stdout);
  
#endif
  return my_keys;
}


/*
 * Computes the size of each bucket by iterating all keys and incrementing
 * their corresponding bucket's size
 */
static inline int * count_local_bucket_sizes(
  const uninitialized_vector<KEY_TYPE> &my_keys,
  dash::NArray<int, 2>        &bucket_sizes)
{
  timer_start(&timers[TIMER_BCOUNT]);

  const KEY_TYPE *__restrict mk      = my_keys.data();
  int *__restrict local_bucket_sizes = bucket_sizes.lbegin();
  std::fill(local_bucket_sizes, local_bucket_sizes + NUM_BUCKETS, 0);

  for(unsigned int i = 0; i < NUM_KEYS_PER_PE; ++i){
    const uint32_t bucket_index = mk[i]/BUCKET_WIDTH;
    local_bucket_sizes[bucket_index]++;
  }

  timer_stop(&timers[TIMER_BCOUNT]);

#ifdef DEBUG
  
  char msg[1024];
  sprintf(msg,"Rank %d: local bucket sizes: ", my_rank);
  for(int i = 0; i < NUM_BUCKETS; ++i){
    if(i < PRINT_MAX)
    sprintf(msg + strlen(msg),"%d ", local_bucket_sizes[i]);
  }
  sprintf(msg + strlen(msg),"\n");
  printf("%s",msg);
  fflush(stdout);
  
#endif

  return local_bucket_sizes;
}

/*
 * Places local keys into their corresponding local bucket.
 * The contents of each bucket are not sorted.
 */
static inline void bucketize_local_keys(
  const uninitialized_vector<KEY_TYPE> &my_keys,
  dash::NArray<int, 3>        &buckets,
  dash::NArray<int, 2>        &bucket_sizes)
{
  timer_start(&timers[TIMER_BUCKETIZE]);

  const KEY_TYPE * __restrict mk = my_keys.data();

  // reset bucket sizes
  dash::fill(bucket_sizes.begin(), bucket_sizes.end(), 0);

  int *__restrict lbs = bucket_sizes.lbegin();
  auto bucket_size    = buckets.extent(2);
  int *__restrict lbp = buckets.lbegin();

  for(unsigned int i = 0; i < NUM_KEYS_PER_PE; ++i){
    const KEY_TYPE key = mk[i];
    const uint32_t bucket_id = key / BUCKET_WIDTH;
    uint32_t index = lbs[bucket_id]++;
    assert(index < NUM_KEYS_PER_PE);
//    lb(bucket_id, index) = key;
//    std::cout << my_rank << ": Putting key " << key << " to "
//              << bucket_id * bucket_size + index
//              << " (" << bucket_id << " * " << bucket_size << " + " << index
//              << std::endl;
    lbp[bucket_id * bucket_size + index] = key;
  }

  dash::barrier();

  timer_stop(&timers[TIMER_BUCKETIZE]);

#ifdef DEBUG
  
  char msg[1024];
  sprintf(msg,"Rank %d: local bucketed keys: ", my_rank);
  for(int i = 0; i < NUM_KEYS_PER_PE; ++i){
    if(i < PRINT_MAX)
    sprintf(msg + strlen(msg),"%d ", my_local_bucketed_keys[i]);
  }
  sprintf(msg + strlen(msg),"\n");
  printf("%s",msg);
  fflush(stdout);
  
#endif
}

/*
 * Each PE sends the contents of its local buckets to the PE that owns that bucket.
 */
static inline uninitialized_vector<KEY_TYPE> exchange_keys(
  dash::NArray<int, 3> &buckets,
  dash::NArray<int, 2> &bucket_sizes,
  long long int        &my_bucket_size)
{
  timer_start(&timers[TIMER_ATA_KEYS]);

  auto myid = dash::myid();

  uninitialized_vector<int> my_bucket_sizes(dash::size());

  // TODO: This is really naive (O(n^2) communications).
  //       Can we do this more elgantly? sth like dash::slice?
  // TODO: At least this should use async operations
  //       (requires .async on dash::NArray)
  for (int i = 0; i < dash::size(); i++) {
    my_bucket_sizes[i] = bucket_sizes(i, myid);
  }

  my_bucket_size = std::accumulate(
                     my_bucket_sizes.begin(), my_bucket_sizes.end(), 0);

  uninitialized_vector<KEY_TYPE> my_bucket_keys(my_bucket_size);

  int offset = 0;

  for (int i = 0; i < dash::size(); i++) {
//    std::cout << "Fetching "
//              << my_bucket_sizes[i]
//              << " elements from unit " << i
//              << " into offset " << offset << std::endl;
    // TODO: broken! wait for fixed version, use DART instead
//    dash::copy_async(
//      buckets[i][myid].begin(),
//      buckets[i][myid].begin() + my_bucket_sizes[i],
//      my_bucket_keys.data() + offset);
    dart_get(
      my_bucket_keys.data() + offset,
      buckets(i, myid, 0).dart_gptr(),
      my_bucket_sizes[i], DART_TYPE_INT);
    offset += my_bucket_sizes[i];
  }

  buckets.flush_all();

//  int i = 0;
//  for (auto val : my_bucket_keys) {
//    std::cout << my_rank << ": my_bucket_keys[" << i++ << "] = " << val << std::endl;
//  }

  dash::barrier();
  timer_stop(&timers[TIMER_ATA_KEYS]);

#ifdef DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  char msg[1024];
  sprintf(msg,"Rank %d: Bucket Size %lld | | Keys after exchange:", my_rank, *my_bucket_size);
  for(int i = 0; i < *my_bucket_size; ++i){
    if(i < PRINT_MAX)
    sprintf(msg + strlen(msg),"%d ", my_bucket_keys[i]);
  }
  sprintf(msg + strlen(msg),"\n");
  printf("%s",msg);
  fflush(stdout);
#endif

  return my_bucket_keys;
}


/*
 * Counts the occurence of each key in my bucket. 
 * Key indices into the count array are the key's value minus my bucket's 
 * minimum key value to allow indexing from 0.
 * my_bucket_keys: All keys in my bucket unsorted [my_rank * BUCKET_WIDTH, (my_rank+1)*BUCKET_WIDTH)
 */
static inline std::vector<CNT_TYPE>
count_local_keys(const uninitialized_vector<KEY_TYPE>& my_bucket_keys,
                 const long long int my_bucket_size)
{
  std::vector<CNT_TYPE> my_local_key_counts(BUCKET_WIDTH);

  timer_start(&timers[TIMER_SORT]);

  const int my_min_key = my_rank * BUCKET_WIDTH;
  const KEY_TYPE *__restrict mbk = my_bucket_keys.data();
        CNT_TYPE *__restrict mlk = my_local_key_counts.data();

  // Count the occurences of each key in my bucket
  for(int i = 0; i < my_bucket_size; ++i){
    const unsigned int key_index = mbk[i] - my_min_key;

    assert(mbk[i] >= my_min_key);
    assert(key_index < BUCKET_WIDTH);

    ++mlk[key_index];
  }
  timer_stop(&timers[TIMER_SORT]);

#ifdef DEBUG
  
  char msg[4096];
  sprintf(msg,"Rank %d: Bucket Size %lld | Local Key Counts:", my_rank, my_bucket_size);
  for(int i = 0; i < BUCKET_WIDTH; ++i){
    if(i < PRINT_MAX)
    sprintf(msg + strlen(msg),"%d ", my_local_key_counts[i]);
  }
  sprintf(msg + strlen(msg),"\n");
  printf("%s",msg);
  fflush(stdout);
  
#endif

  return my_local_key_counts;
}

/*
 * Verifies the correctness of the sort. 
 * Ensures all keys are within a PE's bucket boundaries.
 * Ensures the final number of keys is equal to the initial.
 */
static int verify_results(std::vector<CNT_TYPE>          &my_local_key_counts,
                          uninitialized_vector<KEY_TYPE> &my_local_keys,
                          const long long int my_bucket_size)
{

  MPI_Barrier(MPI_COMM_WORLD);

  int error = 0;

  const int my_min_key = my_rank * BUCKET_WIDTH;
  const int my_max_key = (my_rank+1) * BUCKET_WIDTH - 1;

  // Verify all keys are within bucket boundaries
  for(int i = 0; i < my_bucket_size; ++i){
    const int key = my_local_keys[i];
    if((key < my_min_key) || (key > my_max_key)){
      printf("Rank %d Failed Verification!\n",my_rank);
      printf("Key: %d is outside of bounds [%d, %d]\n", key, my_min_key, my_max_key);
      error = 1;
    }
  }

  // Verify the sum of the key population equals the expected bucket size
  int bucket_size_test = 0;
  for(unsigned int i = 0; i < BUCKET_WIDTH; ++i){
    bucket_size_test += my_local_key_counts[i];
  }
  if(bucket_size_test != my_bucket_size){
      printf("Rank %d Failed Verification!\n",my_rank);
      printf("Actual Bucket Size: %d Should be %lld\n", bucket_size_test, my_bucket_size);
      error = 1;
  }

  // Verify the final number of keys equals the initial number of keys
  long long int total_num_keys = 0;
  MPI_Allreduce(&my_bucket_size, &total_num_keys, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

  if(total_num_keys != (long long int)(NUM_KEYS_PER_PE * NUM_PES)){
    if(my_rank == ROOT_PE){
      printf("Verification Failed!\n");
      printf("Actual total number of keys: %lld", total_num_keys );
      printf(" Expected %" PRId64 "\n", NUM_KEYS_PER_PE * NUM_PES );
      error = 1;
    }
  }
  return error;
}

/*
 * Gathers all the timing information from each PE and prints
 * it to a file. All information from a PE is printed as a row in a tab seperated file
 */
static void log_times(char * log_file)
{
  FILE * fp = NULL;

  for(int i = 0; i < TIMER_NTIMERS; ++i){
    timers[i].all_times = gather_rank_times(&timers[i]);
    timers[i].all_counts = gather_rank_counts(&timers[i]);
  }

  if(my_rank == ROOT_PE)
  {
    int print_names = 0;
    if(file_exists(log_file) != 1){
      print_names = 1;
    }

    if((fp = fopen(log_file, "a+b"))==NULL){
      perror("Error opening log file:");
      exit(1);
    }

    if(print_names == 1){
      print_run_info(fp);
      print_timer_names(fp);
    }
    print_timer_values(fp);

    report_summary_stats();

    fclose(fp);
  }

}

/*
 * Computes the average total time and average all2all time and prints it to the command line
 */
static void report_summary_stats(void)
{
  
  if(timers[TIMER_TOTAL].seconds_iter > 0) {
    const uint32_t num_records = NUM_PES * timers[TIMER_TOTAL].seconds_iter;
    double temp = 0.0;
    for(unsigned int i = 0; i < num_records; ++i){
      temp += timers[TIMER_TOTAL].all_times[i];
    }
      printf("Average total time (per PE): %f seconds\n", temp/num_records);
  }

  if(timers[TIMER_ATA_KEYS].seconds_iter >0) {
    const uint32_t num_records = NUM_PES * timers[TIMER_ATA_KEYS].seconds_iter;
    double temp = 0.0;
    for(unsigned int i = 0; i < num_records; ++i){
      temp += timers[TIMER_ATA_KEYS].all_times[i];
    }
    printf("Average all2all time (per PE): %f seconds\n", temp/num_records);
  }
}

/*
 * Prints all the labels for each timer as a row to the file specified by 'fp'
 */
static void print_timer_names(FILE * fp)
{
  for(int i = 0; i < TIMER_NTIMERS; ++i){
    if(timers[i].seconds_iter > 0){
      fprintf(fp, "%s (sec)\t", timer_names[i]);
    }
    if(timers[i].count_iter > 0){
      fprintf(fp, "%s_COUNTS\t", timer_names[i]);
    }
  }
  fprintf(fp,"\n");
}

/*
 * Prints all the relevant runtime parameters as a row to the file specified by 'fp'
 */
static void print_run_info(FILE * fp)
{
  fprintf(fp,"DASH\t");
  fprintf(fp,"NUM_PES %" PRIu64 "\t", NUM_PES);
  fprintf(fp,"Max_Key %" PRIu64 "\t", MAX_KEY_VAL); 
  fprintf(fp,"Num_Iters %u\t", NUM_ITERATIONS);

  switch(SCALING_OPTION){
    case STRONG: {
        fprintf(fp,"Strong Scaling: %" PRIu64 " total keys\t", NUM_KEYS_PER_PE * NUM_PES);
        break;
      }
    case WEAK: {
        fprintf(fp,"Weak Scaling: %" PRIu64 " keys per PE\t", NUM_KEYS_PER_PE);
        break;
      }
    case WEAK_ISOBUCKET: {
        fprintf(fp,"Weak Scaling Constant Bucket Width: %" PRIu64 " keys per PE \t", NUM_KEYS_PER_PE);
        fprintf(fp,"Constant Bucket Width: %" PRIu64 "\t", BUCKET_WIDTH);
        break;
      }
    default:
      {
        fprintf(fp,"Invalid Scaling Option!\t");
        break;
      }

  }

#ifdef PERMUTE
    fprintf(fp,"Randomized All2All\t");
#elif INCAST
    fprintf(fp,"Incast All2All\t");
#else
    fprintf(fp,"Round Robin GlobBucketSeq\t");
#endif

    fprintf(fp,"\n");
}

/*
 * Prints all of the timining information for an individual PE as a row
 * to the file specificed by 'fp'. 
 */
static void print_timer_values(FILE * fp)
{
  unsigned int num_records = NUM_PES * NUM_ITERATIONS; 

  for(unsigned int i = 0; i < num_records; ++i) {
    for(int t = 0; t < TIMER_NTIMERS; ++t){
      if(!timers[t].all_times.empty()){
        fprintf(fp,"%f\t", timers[t].all_times[i]);
      }
      if(!timers[t].all_counts.empty()){
        fprintf(fp,"%u\t", timers[t].all_counts[i]);
      }
    }
    fprintf(fp,"\n");
  }
}

/* 
 * Aggregates the per PE timing information
 */ 
static std::vector<double> gather_rank_times(_timer_t * const timer)
{
  if(timer->seconds_iter > 0) {

    const unsigned int num_records = NUM_PES * timer->seconds_iter;

    std::vector<double> all_times(num_records);

    MPI_Allgather(timer->seconds, timer->seconds_iter, MPI_DOUBLE,
                  all_times.data(), timer->seconds_iter, MPI_DOUBLE,
                  MPI_COMM_WORLD);

    return all_times;
  }
  else{
    return std::vector<double>();
  }
}

/*
 * Aggregates the per PE timing 'count' information 
 */
static std::vector<unsigned int> gather_rank_counts(_timer_t * const timer)
{
  if(timer->count_iter > 0){
    const unsigned int num_records = NUM_PES * timer->count_iter;

    std::vector<unsigned int> all_counts(num_records);

    MPI_Allgather(timer->count, timer->count_iter, MPI_UNSIGNED,
                  all_counts.data(), timer->count_iter, MPI_UNSIGNED,
                  MPI_COMM_WORLD);

    return all_counts;
  }
  else{
    return std::vector<unsigned int>();
  }
}
/*
 * Seeds each rank based on the rank number and time
 */
static inline pcg32_random_t seed_my_rank(void)
{
  pcg32_random_t rng;
  pcg32_srandom_r(&rng, (uint64_t) my_rank, (uint64_t) my_rank );
  return rng;
}


/*
 * Tests whether or not a file exists. 
 * Returns 1 if file exists
 * Returns 0 if file does not exist
 */
static int file_exists(char * filename)
{
  struct stat buffer;

  if(stat(filename,&buffer) == 0){
    return 1;
  }
  else {
    return 0;
  }

}

#ifdef PERMUTE
/*
 * Creates a randomly ordered array of PEs used in the exchange_keys function
 */
static void create_permutation_array()
{

  permute_array = malloc( NUM_PES * sizeof(int) );

  for(int i = 0; i < NUM_PES; ++i){
    permute_array[i] = i;
  }

  shuffle(permute_array, NUM_PES, sizeof(int));
}

/*
 * Randomly shuffles a generic array
 */
static void shuffle(void * array, size_t n, size_t size)
{
  char tmp[size];
  char * arr = array;
  size_t stride = size * sizeof(char);
  if(n > 1){
    for(int i = 0; i < (n - 1); ++i){
      size_t rnd = (size_t) rand();
      size_t j = i + rnd/(RAND_MAX/(n - i) + 1);
      memcpy(tmp, arr + j*stride, size);
      memcpy(arr + j*stride, arr + i*stride, size);
      memcpy(arr + i*stride, tmp, size);
    }
  }
}
#endif

