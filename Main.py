import numpy as np
import math


def vdpOsc(alpha, sigma, M, N, TFin):
    PI = math.pi
    P = [0] * 101  # Declare array for recording occurances of a probability condition
    dt = TFin / M

    for i in range(0, N): # Iterate over number of trials
        R = [0] * M
        X = [0] * M
        Y = [0] * M

        for j in range(0, M - 1, 2): # Iterate over number of timesteps
            U1, U2 = np.random.uniform(0, 1), np.random.uniform(0, 1) # Create uniform random number U1 and U2

            # Use U1 and U2 to fill array with gaussian r.n. using Box-Muller method
            R[j] = math.sqrt(-2 * math.log(U1)) * math.cos(2 * PI * U2)
            R[j + 1] = math.sqrt(-2 * math.log(U1)) * math.cos(2 * PI * U2)

        X[0] = 0.1 * R[0] #Initial condition X_t0 = 0.1*N(0,1)
        Y[0] = 0.0 # Initial condition Y_t0 = 0

        for j in range(1, M): # Iterate over number of timesteps
            dW = math.sqrt(dt) * R[j] #Declare random noise coeefficient dW = squareroot(dt)*N(0,1)

            # Update new X and Y co-ords.
            X[j] = X[j - 1] + Y[j - 1] * dt
            Y[j] = Y[j - 1] + ((alpha * alpha - X[j - 1] * X[j - 1]) * X[j - 1] - Y[j - 1]) * dt + sigma * X[j - 1] * dW

        for k in range(1, 101): # Iterate over solutions to see if within the bounds of convergence
            loc = k * 10 - 1
            P[k] += (((X[loc] + alpha) * (X[loc] + alpha) + Y[loc] * Y[loc]) <= 0.25) or (((X[loc] - alpha) * (X[loc] - alpha) + Y[loc] * Y[loc]) <= 0.25)

    # Calculate the probability of being within the bounds of convergence at 100 evenly spaced points in time
    Prob = [x / N for x in P] 

    return Prob


ans = vdpOsc(1, 0.5, 1000, 1000, 10)

print(ans)
