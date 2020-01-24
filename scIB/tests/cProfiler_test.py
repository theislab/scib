import math
import random
from joblib import Parallel, delayed
import cProfile
from pstats import Stats

prof = cProfile.Profile()

def is_prime(num):
    if (num==2) or (num==3):
        return True
    if num==1:
        return False
    if (num % 2) == 0:
        return False

    root = int(math.sqrt(num))

    for elt in range(5,root):
        if (num % elt) == 0:
            return False
    return True


res = [random.randrange(100000000) for i in range(10000000)]

#def checkprime(res):
#    for i in res:
#        is_prime(i)

def runPrime(nJ):
    n2 = Parallel(n_jobs=nJ)(delayed(is_prime)(i) for i in res)

#prof.runcall(runPrime,1)
#print(f'runtime_1: {round(Stats(prof).total_tt,0)} s')

prof.runcall(runPrime,20)
print(f'runtime_2: {round(Stats(prof).total_tt,0)} s')


prof = cProfile.Profile()
prof.runcall(runPrime,40)
print(f'runtime_4: {round(Stats(prof).total_tt,0)} s')
