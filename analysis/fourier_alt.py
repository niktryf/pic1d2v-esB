import numpy as np
import matplotlib.pyplot as plt

### Find grid size and timesteps from input file
inputFile = open("../input.txt")
for line in inputFile:
    s = line.split()
    if s:
        if s[0] == 'S':
            size = int(s[1])
        if s[0] == 'T':
            totalTime = float(s[1])
            timestep = float(s[2])
            interval = float(s[3])
            steps = int( (totalTime/timestep)/interval + 2) 

print "steps: ", steps
print "size: ", size
timedata_x = np.zeros(steps)
f = open('../output/E1D.txt')
t = 0
for line in f:
    s = line.split()
    if s:
        p = int(s[0])
        timedata_x[t] = float(s[1])
        print p, timedata_x[t]
    else:
        t+=1
        print " "
