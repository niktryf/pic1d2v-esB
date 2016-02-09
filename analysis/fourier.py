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
            steps = int( (totalTime/timestep)/interval + 1)
            if (totalTime/timestep)%interval != 0:
                steps+=1

### Allocate array for grid data:
### rows -> gridpoints
### each column is a different output
timedata_x = np.zeros(shape = (size, steps))
#timedata_y = np.zeros(shape = (size, steps))

f = open('../output/E1D.txt')
t = 0
for line in f:
    s = line.split()
    if s:
        p = int(s[0])
        timedata_x[p][t] = float(s[1])
    else:
        t+=1

#print timedata_x

t = np.arange(0, totalTime+timestep, timestep*interval) # range: (start, end-not including-, step)

# Prepare figure
fig = plt.figure()
ax = plt.subplot(111)

# Do FFT
for cell in range(0, size):
    sp = np.fft.fft(timedata_x[cell], len(t)) # fft
    sp_abs = np.absolute(sp) # Calculate absolute values
    sp_abs_norm = sp_abs/len(sp) # Normalize

    freq = np.fft.fftfreq(len(t), d=timestep*interval) #frequencies (1.0/t)
    omega = 2.0*np.pi*freq #rad freq
    plt.plot(omega, sp_abs_norm)
    
# Axes Labels
ax.set_xlabel('$\omega$ (rad/sec)')
ax.set_ylabel('Amplitude')
# Make Grids Visible
plt.grid(b=True, which='major', color='k', linestyle='-')
plt.grid(b=True, which='minor', color='k', linestyle='--')
plt.minorticks_on()
plt.show()

############### Test: Print highest frequency component etc (UNDER CONSTRUCTION)
#print sp_abs_norm.argmax() # max argument
#print np.amax(sp_abs_norm) # max of array
#print omega[0:50] #50 first omegas
#print omega[sp_abs_norm.argmax()] #index of omega corresponding to max value of sp

