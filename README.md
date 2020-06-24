# Serial-VanDerPol-Oscillator


    /*
    int main(int argc, char* argv[])

    Inputs:
            int argc: Contains the number of arguements given to the program,
                            function takes 4 arguements with and optional 5th
                            arguement for the seed of the random number generator.
            
            char* argv[]: An array of strings containing the function's arguements,
                            arguements should be given in the following order.

                    alpha: 1 - Half distance between two stable attractors of the spring mass system
                    sigma: 0.5 - Strength of the random noise coefficient
                    M: 1000 - Number of time steps
                    N: 100000 - Number of trials
                    seed: (optional) - seed for the random number generator

    Outputs:

            Prob: Returns an array for the probability of being within...

            Also prints to command prompt the arguements given to the function, the random seed used, and total run time

           Solves the VDP Oscillator using the Euler-Maruyama method,
           and record the updates the number of occurances where the soultion is within specified bounds
           Takes as arguements alpha, sigma, # of time steps M, the final time run time TFin, and record P
    */
    


// Solves the VDP Oscillator using the Euler-Maruyama method,
// and record the updates the number of occurances where the soultion is within specified bounds
// Takes as arguements alpha, sigma, # of time steps M, the final time run time TFin, and record P
