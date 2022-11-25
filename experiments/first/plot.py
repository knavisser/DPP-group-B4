import matplotlib.pyplot as plt
import numpy as np

# read in the data
datapoints = [1000, 10000, 100000, 1000000, 10000000]

data1 = np.array([])
data1Seq = np.array([])


with open('total3.1.txt', 'r') as f:
    for line in f:
        data1 = np.append(data1, float(line))

with open('total3.1_seq.txt', 'r') as f:
    for line in f:
        data1Seq = np.append(data1Seq, float(line))


data1 = data1.reshape(5, 10)
data1Seq = data1Seq.reshape(5, 10)


print(data1)
print(data1Seq)

# calculate the mean and standard deviation
mean = np.mean(data1, axis=1)
std = np.std(data1, axis=1)

mean2 = np.mean(data1Seq, axis=1)
std2 = np.std(data1Seq, axis=1)




# calculate speedup by dividing sequential time by parallel time
speedup = mean2 / mean


# plot the speedup
plt.plot(datapoints, speedup, 'o-')
plt.xlabel('Number of Datapoints')
plt.ylabel('Speedup')
plt.title('Speedup of Parallel vs Sequential')
plt.legend()
plt.savefig('speedup.png')
plt.show()



