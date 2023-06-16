
import numpy as np
import matplotlib.pyplot as plt

# Generate n pseudo-random numbers with(mu,sigma,n)
x = np.random.normal(0, 100, 100) 

x_min, x_max = np.min(x), np.max(x)

N_MIN = 4  # Min number of bins (integer), must be > 1
N_MAX = 50 # Max number of bins (integer)

N = np.arange(N_MIN,N_MAX) # number of bins
D = (x_max-x_min)/N        # Bin size vector
C = np.zeros(np.size(D))

# Computation of the cost function
for i in range(np.size(N)):
    edges = np.linspace(x_min,x_max,N[i]+1) # Bin edges
ki = plt.hist(x,edges,alpha=0.5)[0]     # Count number of events in bins
k = np.mean(ki)                         # Mean of event count
v = np.sum((ki-k)**2)/N[i]              # Variance of event count
C[i] = (2*k-v)/((D[i])**2)              # Cost Function

# Optimal bin size Selection
cmin = np.min(C)
idx  = np.where(C == cmin)[0][0]
optD = D[idx] 

# Plotting
edges = np.linspace(x_min,x_max,N[idx]+1)
plt.figure(figsize=(12,4))
plt.subplot(121)
plt.hist(x,edges)
plt.title("Histogram")
plt.ylabel("Frequency")
plt.xlabel("Value")
#plt.savefig('Hist.png')

plt.subplot(122)
plt.plot(D,C,'.',alpha=0.5)
plt.plot(optD,cmin,'*r');
#plt.savefig('Fobj.png')
