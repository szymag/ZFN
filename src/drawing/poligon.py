import numpy as np

for i in np.array([0, 16, 32, 48, 54, 70, 90]):
    tmp = np.loadtxt('dys_' + str(i) + '.dat')[0][1:]
    a = np.argmin(np.abs(tmp - 5.2e9))
    print(a)