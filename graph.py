import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import numpy as np


temp = "Temp.txt"
Omg = "Omg.txt"
Psi = "Psi.txt"
temp_data = []
Omg_data = []
Psi_data = []

with open(temp, 'r') as temp_file, open(Omg, 'r') as Omg_file, open(Psi, 'r') as Psi_file:
    for line in temp_file:
        temp_data.append(list(map(float, line.rstrip().split(' '))))
    for line in Omg_file:
        Omg_data.append(list(map(float, line.rstrip().split(' '))))
    for line in Psi_file:
        Psi_data.append(list(map(float, line.rstrip().split(' '))))

example_temp = "temp_all.txt"
example_data = []
with open(example_temp, 'r') as example:
    for line in example:
        example_data.append(list(map(float, line.rstrip().split(' ')))) 


# Make data
N = len(temp_data[0])
M = len(temp_data)
X = np.arange(0, 1, 1/N)
Y = np.arange(0, 1, 1/M)
X, Y = np.meshgrid(X, Y)
Z_temp = np.array(temp_data)
Z_omg = np.array(Omg_data)
Z_psi = np.array(Psi_data)

example_array = np.zeros((1000, M, N))
for t in range(len(example_data)):
    for i in range(N*M):
        example_array[t][i // M][i % N] = example_data[t][i]

TX = np.arange(0, 1, 1/N)
TY = np.arange(0, 1, 1/M)
TX, TY = np.meshgrid(TX, TY)
Z_example = example_array
level = 10
Z_temp_levels = np.linspace(np.min(Z_example), np.max(Z_example), level)

# Plot the surface
font = {'size'   : 10}
matplotlib.rc('font', **font)
fig, ax = plt.subplots(1, 3, subplot_kw={'projection': '3d'})

for item in ax:
    item.set_xlabel("$x$")
    item.set_ylabel("$y$")

ax[0].plot_surface(X, Y, Z_temp, cmap=cm.magma)
ax[0].set(title="$T$")

ax[1].plot_surface(X, Y, Z_psi, cmap=cm.magma)
ax[1].set(title="$psi$")

ax[2].plot_surface(X, Y, Z_omg, cmap=cm.magma)
ax[2].set(title="$omega$")


plt.show()

fig, ax = plt.subplots(1, 3, figsize=(10,10))

for item in ax:
    item.set_xlabel("$x$")
    item.set_ylabel("$y$")
    item.set_aspect('equal', adjustable='box')


level = 10
CS = ax[0].contour(X, Y, Z_temp, level)
ax[0].clabel(CS, inline=True)
ax[0].set(title="$T$")

CS = ax[1].contour(X, Y, Z_psi, level)
ax[1].clabel(CS, inline=True)
ax[1].set(title="$psi$")

CS = ax[2].contour(X, Y, Z_omg, level)
ax[2].clabel(CS, inline=True)
ax[2].set(title="$omega$")


plt.show()

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(projection = '3d')
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set(title="$T$")
ax.view_init(elev=90, azim=-90, roll=0)
for t in range(0, 1000):
    # optionally clear axes and reset limits
    plt.gca().cla()
    ax.plot_surface(TX, TY, Z_example[t], cmap=cm.hot)
    fig.canvas.draw()
    plt.pause(0.001)
    # plt.savefig(f"./anim/flow{t}.png")
plt.show()