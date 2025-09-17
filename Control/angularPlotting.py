import matplotlib.pyplot as plt
import numpy as np
import csv


allData = []

with open('testing.csv', 'r') as file:
    csv_reader = csv.reader(file)
    i = 0
    for row in csv_reader:
        if(i > 0):
            newArr = []
            for d in row:
                newArr.append(float(d))
            allData.append(newArr)
        i = i + 1
allData = np.array(allData)
ts = allData[:,0]
w1 = allData[:,1]
w2 = allData[:,2]
w3 = allData[:,3]
a1 = allData[:,4]
a2 = allData[:,5]
a3 = allData[:,6]

plt.plot(ts,w3)
plt.xlabel("Time")
plt.ylabel("OMEGA")
plt.xlim((0,16))
plt.ylim((-3,3))
plt.show()