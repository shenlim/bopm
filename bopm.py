# Binomial Options Pricing Model (BOPM) - Nov 2018. Author: Shen Lim.
# Copyright 2018, Shen Lim, All Rights Reserved.
import numpy as np
import time # Optional to measure elapsed time.

start_time = time.time()
def bopm():
    S = float(input("Spot stock price? "))
    K = float(input("Strike price? "))
    T = float(input("Time to maturity (years)? "))
    r = float(input("Risk-free rate (decimals)? "))
    sigma = float(input("Volatility (decimals)? "))
    N = int(input("Number of steps? "))
    call = input("Call or Put? ")
    q = float(input("Dividends? Enter 0 if none. "))
    A = input("American or European? ")
    
    deltaT = T / N
    u = np.exp(sigma * np.sqrt(deltaT))
    d = 1 / u
    p = (np.exp((r - q)* deltaT) - d) / (u - d)

    np.set_printoptions(formatter={'float': lambda x: "{0:0.4f}".format(x)})

    # Stock tree
    tree = np.zeros((N + 1, N + 1), dtype = float)
    for i in range(0, N + 1):
        for j in range(0, i + 1):
            tree[i, j] = (S * u ** (j) * d ** (i - j))
    print("\nStock Tree:\n", tree, "\n")

    # Exercise tree
    exC = np.zeros((N + 1, N + 1), dtype = float)
    for i in range(0, N + 1):
        for j in range(0, i + 1):
            exC[i, j] = np.maximum(((S * u ** (j) * d ** (i - j)) - K), 0)

    exP = np.zeros((N + 1, N + 1), dtype = float)
    for i in range(0, N + 1):
        for j in range(0, i + 1):
            exP[i, j] = np.maximum((K - (S * u ** (j) * d ** (i - j))), 0)
    
    if (call.lower() == "c") or (call.lower() == "call"):
        print("Exercise Tree:\n", exC, "\n")

    elif (call.lower() == "p") or (call.lower() == "put"):
        print("Exercise Tree:\n", exP, "\n")

    else:
        return ("Invalid input.")

    # Option tree
    opTree = np.zeros((N + 1, N + 1), dtype = float)
    
    if (call.lower() == "c") or (call.lower() == "call"):
        opTree[N,] = exC[N,] # Fix option value at final node (already  computed in the exercise tree).

    elif (call.lower() == "p") or (call.lower() == "put"):
        opTree[N,] = exP[N,]

    else:
        return ("Invalid input.")

    if (A.lower() == "e") or (A.lower() == "european"):
        for i in range(N - 1, -1, -1): # Recursively calculate option value at each node (except final).
            for j in range(0, i + 1):
                opTree[i, j] = (opTree[i + 1, j + 1] * p + opTree[i + 1, j] * (1 - p)) / np.exp(r * deltaT)

        delta = (opTree[1, 1] - opTree[1, 0]) / (tree [1, 1] - tree[1, 0])
        
        print("Option Tree:\n", opTree, "\n",
              "\nP           : ", p, "\n",
              "\nDelta       : ", delta, "\n",
              "\nOption Price: ", opTree[0, 0])

    elif (A.lower() == "a") or (A.lower() == "american"):
        for i in range(N - 1, -1, -1):
            for j in range(0, i + 1):
                if (call.lower() == "c") or (call.lower() == "call"): exPayoff = exC[i, j]
                elif (call.lower() == "p") or (call.lower() == "put"): exPayoff = exP[i, j]
                else: return("Invalid input.")
                holdPayoff = (opTree[i + 1, j + 1] * p + opTree[i + 1, j] * (1 - p)) / np.exp(r * deltaT)
                opTree[i, j] = np.maximum(exPayoff, holdPayoff)

        delta = (opTree[1, 1] - opTree[1, 0]) / (tree [1, 1] - tree[1, 0])
        
        print("Option Tree:\n", opTree, "\n",
              "\nP           : ", p, "\n",
              "\nDelta       : ", delta, "\n",
              "\nOption Price: ", opTree[0, 0])

    else:
        return("Invalid input.")

    print("\nElapsed time: %s seconds" % ((time.time() - start_time) / 1000))
