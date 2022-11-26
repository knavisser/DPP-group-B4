import matplotlib.pyplot as plt
import numpy as np

# read in the data
datapoints = [1000, 10000, 100000, 1000000, 10000000]

data1 = np.array([])
data1Seq = np.array([])

with open('totalP3.1_8processes.txt', 'r') as f:
    for line in f:
        data1 = np.append(data1, float(line))
with open('totalS3.1txt', 'r') as f:
    for line in f:
        data1Seq = np.append(data1Seq, float(line))

data1 = data1.reshape(5, 10)
data1Seq = data1Seq.reshape(5, 10)

# calculate the mean and standard deviation
mean = np.mean(data1, axis=1)
std8 = np.std(data1, axis=1)

mean2sequential = np.mean(data1Seq, axis=1)
std2sequential = np.std8(data1Seq, axis=1)

speedup = mean2sequential/mean
# normalise speedup so that it starts at 0 and ends at 1
speedup3_1 = (speedup - speedup[0]) / (speedup[-1] - speedup[0])


data2 = np.array([])
data2Seq = np.array([])
with open('totalP1.1.txt', 'r') as f:
    for line in f:
        data2 = np.append(data2, float(line))
with open('totalS1.1.txt', 'r') as f:
    for line in f:
        data2Seq = np.append(data2Seq, float(line))

data2 = data2.reshape(5, 10)
data2Seq = data2Seq.reshape(5, 10)

# calculate the mean and standard deviation
mean = np.mean(data2, axis=1)
std = np.std(data2, axis=1)

mean2 = np.mean(data2Seq, axis=1)
std2 = np.std(data2Seq, axis=1)

speedup = mean2/mean
# normalise speedup so that it starts at 0 and ends at 1
speedup1_1 = (speedup - speedup[0]) / (speedup[-1] - speedup[0])

data3 = np.array([])
data3Seq = np.array([])

with open('totalP1.2.txt', 'r') as f:
    for line in f:
        data3 = np.append(data3, float(line))
with open('totalS1.2.txt', 'r') as f:
    for line in f:
        data3Seq = np.append(data3Seq, float(line))

data3 = data3.reshape(5, 10)
data3Seq = data3Seq.reshape(5, 10)

# calculate the mean and standard deviation
mean = np.mean(data3, axis=1)
std = np.std(data3, axis=1)

mean2 = np.mean(data3Seq, axis=1)
std2 = np.std(data3Seq, axis=1)

speedup = mean2/mean
# normalise speedup so that it starts at 0 and ends at 1
speedup1_2 = (speedup - speedup[0]) / (speedup[-1] - speedup[0])


data4 = np.array([])
with open('totalP3.1_1proces.txt', 'r') as f:
    for line in f:
        data4 = np.append(data4, float(line))

data4 = data4.reshape(5, 10)

# calculate the mean and standard deviation
mean = np.mean(data4, axis=1)
std = np.std(data4, axis=1)

speedup = mean2sequential/mean

# normalise speedup so that it starts at 0 and ends at 1
speedup3_1_1 = (speedup - speedup[0]) / (speedup[-1] - speedup[0])

data5 = np.array([])
with open('totalP3.1_4processes.txt', 'r') as f:
    for line in f:
        data5 = np.append(data5, float(line))

data5 = data5.reshape(5, 10)

mean = np.mean(data5, axis=1)
std = np.std(data5, axis=1)

speedup = mean2sequential/mean
# normalise speedup so that it starts at 0 and ends at 1
speedup3_1_4proces = (speedup - speedup[0]) / (speedup[-1] - speedup[0])


# # plot the speedup
# plt.plot(datapoints, speedup3_1, label='MPI: 8 nodes | 8 processes', marker='o')
# plt.plot(datapoints, speedup1_1, label='Pthreads: 2 threads', marker='o')
# plt.plot(datapoints, speedup1_2, label='OpenMP: 2 threads', marker='o')
# plt.title('Speedup of a wave equation with several parallelisation methods')
# plt.xlabel('Number of datapoints')
# plt.ylabel('Speedup: Sequential time (s) / Parallel time (s)')
# plt.xscale('log')
# plt.legend()
# plt.savefig('speedupVersus.png')
# plt.show()

# errorplot
plt.errorbar(datapoints, speedup3_1, yerr=std, label='MPI: 8 nodes | 8 processes', marker='o')
plt.errorbar(datapoints, speedup3_1_1, yerr=std2, label='MPI: 8 nodes | 1 process', marker='o')





