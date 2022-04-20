# Constants for CNT Generation

import numpy as np

EPS = 1.0e-8

# Graphene 

A_CC = 1.418
A = np.sqrt(3.0) * A_CC

D = np.array([-A_CC, 0.0])

A1 = 0.5 * A * np.array([np.sqrt(3), 1.0])
A2 = 0.5 * A * np.array([np.sqrt(3), -1.0])
