import itertools
import numpy as np

H_X = ...        # r x n binary array -> grab from QR codes directly (generate using MAGMA online if necessary)
z_log = ...      # length-n binary array for logical Z -> (pg 2 or 3 of the paper)

n = H_X.shape[1]
w_max = 50 # might have to lower to 30

w_min = None
A_wmin = 0

for w in range(1, w_max+1):
    found_this_weight = False
    for positions in itertools.combinations(range(n), w):
        e = np.zeros(n, dtype=int)
        e[list(positions)] = 1

        if np.all((H_X @ e) % 2 == 0):         # undetected
            if (z_log @ e) % 2 == 1:          # harmful
                found_this_weight = True
                if w_min is None:
                    w_min = w
                    A_wmin = 1
                elif w == w_min:
                    A_wmin += 1
    if found_this_weight:
        break
