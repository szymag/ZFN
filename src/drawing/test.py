import matplotlib.pyplot as plt
plt.rcdefaults()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def rang(x):
    return [2 * i for i in range(x // 2)]

file0 = np.loadtxt('0tm5.txt') / 10e9
file1 = np.loadtxt('1tm5.txt') / 10e9
file2 = np.loadtxt('0tm6.txt') / 10e9
file3 = np.loadtxt('1tm6.txt') / 10e9
file4 = np.loadtxt('0tm7.txt') / 10e9
file5 = np.loadtxt('1tm7.txt') / 10e9
file6 = np.loadtxt('0tm8.txt') / 10e9
file7 = np.loadtxt('1tm8.txt') / 10e9

ax = plt.gca()

patches = []
# add a rectangle


for i in [0, 2, 4, 6]:
    tmp = eval('file' + str(i))
    tmp1 = eval('file' + str(i+1))
    for j in rang(800 + i*25):
        rect = mpatches.Rectangle((tmp[j], 5*i), (tmp1[j] - tmp[j])*1.01, 10)
        patches.append(rect)
        rect = mpatches.Rectangle((tmp1[j + 1], 5*i), (tmp[j + 1] - tmp1[j + 1])*1.01, 10)
        patches.append(rect)


collection = PatchCollection(patches,  alpha=1, color='green')
ax.add_collection(collection)
plt.ylim([0, 40])
plt.xlim([1.1, 3])
ax.axes.get_yaxis().set_visible(False)
plt.title('tytuł')
plt.xlabel('frequency [GHz]')
plt.text(1.2, 5, 'TM5', style='italic', bbox={'facecolor':'white', 'alpha':0.4, 'pad':10})
plt.text(1.2, 15, 'TM6', style='italic', bbox={'facecolor':'white', 'alpha':0.4, 'pad':10})
plt.text(1.2, 25, 'TM7', style='italic', bbox={'facecolor':'white', 'alpha':0.4, 'pad':10})
plt.text(1.2, 35, 'TM8', style='italic', bbox={'facecolor':'white', 'alpha':0.4, 'pad':10})
plt.savefig('tm.svg')



