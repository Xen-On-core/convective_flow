import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
from matplotlib.animation import FuncAnimation
import time

from PIL import Image

import sys
import io
import numpy as np


def help():
    print("%s hadles saved data in .txt files and save plane graph, \n" % sys.argv[0])
    print("     contours and gif animation (or displays real-time render).\n\n")
    print("Usage:\n")
    print("     %s OPTION       \n" % sys.argv[0])
    print("\nOptions:\n")
    print("     -R, --real-time             \n")
    print("     -F, --file-name FILENAME    \n")
    print("     -?, --help                  \n")


if len(sys.argv) < 2 or len(sys.argv) > 3 or sys.argv[1] in ['-?', '--help']:
    help()
    exit(0)

realtime = False
if sys.argv[1] in ['-R', '--real-time']:
    realtime = True
elif sys.argv[1] in ['-F', '--file-name']:
    save_files_name = sys.argv[2]
else:
    help()
    exit(0)

datadir = './output_data'

temp = f"{datadir}/temperature.txt"
Omg = f"{datadir}/omega.txt"
Psi = f"{datadir}/psi.txt"
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

example_temp = f"{datadir}/temperature_time_data.txt"
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

example_array = np.zeros((len(example_data), M, N))
for t in range(len(example_data)):
    for i in range(N*M):
        example_array[t][i // M][i % N] = example_data[t][i]

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

if realtime:
    plt.show()
else:
    plt.savefig(f'./images/{save_files_name}_plane.png')
plt.close()

fig, ax = plt.subplots(1, 3, figsize=(10,10))

for item in ax:
    item.set_xlabel("$x$")
    item.set_ylabel("$y$")
    item.set_aspect('equal', adjustable='box')

CS = ax[0].contour(X, Y, Z_temp, level)
ax[0].clabel(CS, inline=True)
ax[0].set(title="$T$")

CS = ax[1].contour(X, Y, Z_psi, level)
ax[1].clabel(CS, inline=True)
ax[1].set(title="$psi$")

CS = ax[2].contour(X, Y, Z_omg, level)
ax[2].clabel(CS, inline=True)
ax[2].set(title="$omega$")

if realtime:
    plt.show()
else:
    plt.savefig(f'./images/{save_files_name}_contour.png')
plt.close()

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection = '3d')
ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set(title="$T, time = 0.0$")
ax.view_init(elev=90, azim=-90, roll=0)

frames_num = len(example_data)

# def update(i):
#     start_time = time.time()
#     if i+1 % 100 == 0:
#         print(f"=== IMAGE ITERATION {i+1} ===")
#         print("--- %s seconds ---" % (time.time() - start_time))
#     plt.gca().cla()
#     ax.set(title=f"$T$, time = {i / frames_num}")
#     ax.plot_surface(TX, TY, Z_example[i], cmap=cm.hot)

# anim = FuncAnimation(fig, update, frames = frames_num, interval=10)
# anim.save(f'./anim/{save_files_name}.gif')

if realtime:
    # This code need for show real-time animation
    for t in range(0, frames_num):
        plt.gca().cla()
        ax.set(title=f"$T$, time = {t / frames_num}")
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.plot_surface(X, Y, Z_example[t], cmap=cm.magma)
        fig.canvas.draw()
        plt.pause(0.01)
    plt.show()
else:
    # This code need for save gif of real-time animation
    #   without displaying figure.
    print("=== SAVE IMAGES ===")
    image = []
    times = []
    start_time = time.time()
    for t in range(0, frames_num):
        if (t+1) % 100 == 0:
            print("=== IMAGE ITERATION %d ===" % (t+1))
            print("--- %f seconds ---" % (time.time() - start_time))
            times.append((time.time() - start_time))
            start_time = time.time()
        plt.gca().cla()
        ax.set(title=f"$T$, time = {t / frames_num}")
        ax.plot_surface(X, Y, Z_example[t], cmap=cm.magma)
        buf = io.BytesIO()
        fig.savefig(buf)
        buf.seek(0)
        image.append(Image.open(buf))
    print("--- GLOBAL SAVE TIME ---")
    print("--- %f seconds ---" % (sum(times)))

    print("=== SAVE GIF ===")
    image[0].save(
        f'./anim/{save_files_name}.gif', 
        save_all = True, 
        append_images = image[1:], 
        optimize = False, 
        duration = 10,
        loop = 0)
    buf.close()
    print("\tSuccessful!\n")