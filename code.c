/**

 * University of Pittsburgh

 * Department of Computer Science

 * CS1645: Introduction to HPC Systems

 * Student: 

 * Instructor Esteban Meneses, PhD (emeneses@pitt.edu)

 * MPI particle-interaction code. 

 */

#
include "mpi.h"

#
include < stdio.h >

  #include < stdlib.h >

  #include < math.h >

  #include < iostream >

  #include < fstream >

  //#include <stdio.h>

  //#include <stdlib.h>

  //#include <mpi.h>

  #include < stddef.h >

  using namespace std;

#
define TAG 7

# define CONSTANT 777

// Particle-interaction constants

# define A 10250000.0

# define B 726515000.5

# define MASS 0.1

# define DELTA 1

// Random initialization constants

# define POSITION 0

# define VELOCITY 1

// Structure for shared properties of a particle (to be included in messages)

struct Particle {

  float x;

  float y;

  float mass;

  float fx;

  float fy;

};

// Headers for auxiliar functions

float random_value(int type);

void print_particles(struct Particle * particles, int n);

void interact(struct Particle * source, struct Particle * destination);

void compute_interaction(struct Particle * source, struct Particle * destination, int limit);

void compute_self_interaction(struct Particle * set, int size);

void merge(struct Particle * first, struct Particle * second, int limit);

int read_file(struct Particle * set, int size, char * file_name);

// Main function

main(int argc, char * * argv) {

  int myRank; // Rank of process

  int p; // Number of processes

  int n; // Number of total particles

  int previous; // Previous rank in the ring

  int next; // Next rank in the ring

  int tag = TAG; // Tag for message

  int number; // Number of local particles

  struct Particle * globals; // Array of all particles in the system

  struct Particle * locals; // Array of local particles

  struct Particle * remotes; // Array of foreign particles

  char * file_name; // File name

  MPI_Status status; // Return status for receive

  int j, rounds, initiator, sender;

  double start_time, end_time;

  // checking the number of parameters

  if (argc < 2) {

    printf("ERROR: Not enough parameters\n");

    printf("Usage: %s <number of particles> [<file>]\n", argv[0]);

    exit(1);

  }

  // getting number of particles

  n = atoi(argv[1]);

  // initializing MPI structures and checking p is odd

  MPI_Init( & argc, & argv);

  MPI_Comm_rank(MPI_COMM_WORLD, & myRank);

  MPI_Comm_size(MPI_COMM_WORLD, & p);

  if (p % 2 == 0) {

    p = p - 1;

    if (myRank == p) {

      MPI_Finalize();

      return 0;

    }

  }

  srand(myRank + myRank * CONSTANT);

  // acquiring memory for particle arrays

  number = n / p;

  locals = (struct Particle * ) malloc(number * sizeof(struct Particle));

  remotes = (struct Particle * ) malloc(number * sizeof(struct Particle));

  if (argc == 3) {

    if (myRank == 0) {

      globals = (struct Particle * ) malloc(n * sizeof(struct Particle));

      // YOUR CODE GOES HERE (reading particles from file)

      read_file(globals, n, argv[2]);

    }

    MPI_Scatter(globals, number * 5, MPI_FLOAT, locals, number * 5, MPI_FLOAT, 0, MPI_COMM_WORLD);

  } else {

    // random initialization of local particle array

    for (j = 0; j < number; j++) {

      locals[j].x = random_value(POSITION);

      locals[j].y = random_value(POSITION);

      locals[j].fx = 0.0;

      locals[j].fy = 0.0;

      locals[j].mass = MASS;

    }

  }

  if (myRank == 0) {

    start_time = MPI_Wtime();

  }

  next = (myRank + p + 1) % p;

  previous = (myRank + p - 1) % p;

  for (int i = 0; i < (p - 1) / 2; i++)

  {

    if (i == 0)

    {

      MPI_Send(locals, number * 5, MPI_FLOAT, next, tag, MPI_COMM_WORLD);

    } else

    {

      MPI_Send(remotes, number * 5, MPI_FLOAT, next, tag, MPI_COMM_WORLD);

    }

    MPI_Recv(remotes, number * 5, MPI_FLOAT, previous, tag, MPI_COMM_WORLD, & status);

    compute_interaction(locals, remotes, number);

  }

  int val = ((p - 1) / 2);

  MPI_Send(remotes, number * 5, MPI_FLOAT, ((myRank + p - val) % p), tag, MPI_COMM_WORLD);

  MPI_Recv(remotes, number * 5, MPI_FLOAT, ((myRank + p + val) % p), tag, MPI_COMM_WORLD, & status);

  merge(locals, remotes, number);

  compute_self_interaction(locals, number);

  if (myRank == 0) {

    end_time = MPI_Wtime();

    printf("Duration: %f seconds\n", (end_time - start_time));

  }

  if (argc == 3) {

    MPI_Gather(locals, number * 5, MPI_FLOAT, globals, number * 5, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (myRank == 0) {

      print_particles(globals, n);

    }

  }

  MPI_Finalize();

}

// Function for random value generation

float random_value(int type) {

  float value;

  switch (type) {

  case POSITION:

    value = (float) rand() / (float) RAND_MAX * 100.0;

    break;

  case VELOCITY:

    value = (float) rand() / (float) RAND_MAX * 10.0;

    break;

  default:

    value = 1.1;

  }

  return value;

}

// Function for printing out the particle array

void print_particles(struct Particle * particles, int n) {

  int j;

  printf("Index\tx\ty\tmass\tfx\tfy\n");

  for (j = 0; j < n; j++) {

    printf("%d\t%f\t%f\t%f\t%f\t%f\n", j, particles[j].x, particles[j].y, particles[j].mass, particles[j].fx, particles[j].fy);

  }

}

// Function for printing out the particle array

void print_particles_input(struct Particle * particles, int n) {

  int j;

  printf("x\ty\tmass\n");

  for (j = 0; j < n; j++) {

    printf("%f\t%f\t%f\n", particles[j].x, particles[j].y, particles[j].mass);

  }

}

// Function for computing interaction among two particles

// There is an extra test for interaction of identical particles, in which case there is no effect over the destination

void interact(struct Particle * first, struct Particle * second) {

  float rx, ry, r, fx, fy, f;

  // computing base values

  rx = first - > x - second - > x;

  ry = first - > y - second - > y;

  r = sqrt(rx * rx + ry * ry);

  //printf("%f\n",r);

  if (r == 0.0)

    return;

  f = A / pow(r, 6) - B / pow(r, 12);

  fx = f * rx / r;

  fy = f * ry / r;

  // updating sources's structure

  first - > fx = first - > fx + fx;

  first - > fy = first - > fy + fy;

  // updating destination's structure

  second - > fx = second - > fx - fx;

  second - > fy = second - > fy - fy;

}

// Function for computing interaction between two sets of particles

void compute_interaction(struct Particle * first, struct Particle * second, int size) {

  int j, k;

  //printf("%f\n",first[0].x);

  for (j = 0; j < size; j++) {

    for (k = 0; k < size; k++) {

      interact( & first[j], & second[k]);

    }

  }

}

// Function for computing interaction between two sets of particles

void compute_self_interaction(struct Particle * set, int size) {

  int j, k;

  for (j = 0; j < size; j++) {

    for (k = j + 1; k < size; k++) {

      interact( & set[j], & set[k]);

    }

  }

}

// Function to merge two particle arrays

// Permanent changes reside only in first array

void merge(struct Particle * first, struct Particle * second, int limit) {

  int j;

  for (j = 0; j < limit; j++) {

    first[j].fx += second[j].fx;

    first[j].fy += second[j].fy;

  }

}

// Reads particle information from a text file

int read_file(struct Particle * set, int size, char * file_name) {

  ifstream file(file_name);

  if (file.is_open()) {

    // reading particle values

    for (int i = 0; i < size; i++) {

      file >> set[i].x;

      file >> set[i].y;

      file >> set[i].mass;

      set[i].fx = 0.0;

      set[i].fy = 0.0;

    }

    // closing file

    file.close();

  } else {

    cout << "Error opening file: " << file_name << endl;

    return 1;

  }

  return 0;

}