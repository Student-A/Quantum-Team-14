import numpy as np
import math

TRULY_RANDOM_CHOICES = np.fromfile("./res/random_1M.dat", dtype=bool)

def get_random_choice_generator():
    for i in range(len(TRULY_RANDOM_CHOICES)):
        yield TRULY_RANDOM_CHOICES[i]

RANDOM_GENERATOR = get_random_choice_generator()
        
def get_random_choice():
        next(RANDOM_GENERATOR)
        
def get_random_number(): # between 0 and 1
    return float(int("0b"+"".join([str(int(get_random_choice())) for i in range(32)]), 2))/(0x11111111111111111111111111111111)
        
def classical_random_walk1D(n, coinflips, rightBias):
    current_pos = 0
    for i in range(n):
        if get_random_number() > rightBias:
            current_pos += 1
        else:
            current_pos -= 1
    return current_pos

def quantum_random_walk1D(n, coinflips):
    F = 1/math.sqrt(2) * np.array([[1, 1], [1, -1]])
    
print(get_random_choice())
print(get_random_choice())
print(get_random_choice())
print(get_random_number())
print(get_random_number())
