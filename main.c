#include <stdio.h>
#include  <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#define PI 3.141592654

int main(int argc, char* argv[] )
{
    int alpha;
    float sigma;
    int M;
    int N;
    long int seed;
        
    alpha = atoi(argv[1]);
    sigma = atof(argv[2]);
    M = atoi(argv[3]);
    N = atoi(argv[4]);
    
    double prob[101]; // Declare array to hold of probabilities of being within a a bound of convergence

    int i;
    double P[101]; // Declare array for recording occurances of a probability condition
    for(i = 0; i < 101; ++i)
    {
        P[i] = 0.0;
    }
    
    int TFin = 10; // Final Time
    double dt = (double)TFin/(double)M; // Create size of time step
    
    if (argc > 5) // Intial seed for random numbers is given
    {
        seed = atol(argv[5]);
    }
    else // Intial seed not given, so generate one from /dev/urandom
    {
        FILE* urand = fopen("/dev/urandom", "r");
        fread( &seed, sizeof(long int), 1, urand);
        fclose(urand);
    }

    for (i = 0; i < N; ++i) // Iterate over number of trials
    {
        srand48(seed + i); // Seed uniform random number generator, add i to the seed for each trial so each trial is independent
        double U1, U2; // Declare two uniform random number variables
        double* random = (double*) malloc(sizeof(double)*M); // Dynanmically allocate memory for array of gaussian r.n.
        int i;
        for (i=0; i<(M-1); i+=2) // Loop through number of time steps
        {
            U1 = drand48();// Create uniform random number U1
            U2 = drand48();// Create uniform random number U2
            random[i] = sqrt(-2 * log(U1))*cos(2*PI*U2); // Use U1 and U2 to fill array with gaussian r.n.
            random[i+1] = sqrt(-2 * log(U1))*cos(2*PI*U2); // Using Box-Muler algorithm
        }
        
        double X[M]; // Declare array of X-coord. of length M
        double Y[M]; // Declare array of Y-coord. of length M
        X[0] = 0.1*random[0]; // Initial condition X_t0 = 0.1*N(0,1)
        Y[0] = 0.0; // Initial condition Y_t0 = 0

        int j = 1;
        for (j=1; j<M; ++j) // Loop through number of timesteps M
        {
          double dW = sqrt(dt)*random[j]; // Declare random noise coeefficient dW = squareroot(dt)*N(0,1)

          // Update new X-coord. then update new Y-coord
          X[j] = X[j-1] + Y[j-1]*dt;
          Y[j] = Y[j-1] + ((alpha*alpha - X[j-1]*X[j-1])*X[j-1] - Y[j-1])*dt + sigma*X[j-1]*dW;
        }

        int k = 1;
        int loc = 0; // Declare temporary holding variable location
        for(k=1; k<101; ++k)
        {
            loc = k*10 - 1; //Define location for iterations over k, location is for corresponding location in (X, Y) coord.

            if( ( ((X[loc] + alpha)*(X[loc] + alpha) +  Y[loc]*Y[loc])  <= 0.25 ) || ( ((X[loc] - alpha)*(X[loc] - alpha) + Y[loc]*Y[loc]) <= 0.25 ) )
            {
                P[k]++; // Condition met, add one to p(k)
            }
        }
        
        free((void*)random);  // Free dynamically allocated memory of the random numbers
    }

    for(i = 0; i<101; ++i) // Calculate the probability at time t = i/10 for i = 0, 1, ..., 100
    {
      prob[i] = ((double)P[i])/((double)N);
    }

    // Read out file of probabilites to Prob.out
    FILE *fp;
    fp = fopen( "Prob.out", "w");
    fwrite(&prob, sizeof(double), 101, fp);
    fclose(fp);

    return 0;
}
