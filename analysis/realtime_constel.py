import matplotlib.pyplot as plt
import numpy as np
import threading
import sys
import queue

data_queue = queue.Queue() # thread safe

def read_stdin(): # input thread
    for line in sys.stdin:
        x,y = line.strip().split()
        x,y = float(x), float(y)
        data_queue.put((x,y))

input_thread = threading.Thread(target=read_stdin, daemon=True)
input_thread.start()

plt.ion()
fig, ax = plt.subplots()
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect("equal")
ax.grid(True)
plt.title("Real-time constellation")

theta = np.linspace(0, 2*np.pi, 256)
circle_x = np.cos(theta)
circle_y = np.sin(theta)
ax.plot(circle_x, circle_y, linestyle="--", color="lightgray")

points = []
sc = ax.scatter([], [])

while True:
    updated = False
    while not data_queue.empty():
        points.append(data_queue.get())
        updated = True
    if updated:
        xs,ys = zip(*points)
        sc.set_offsets(np.column_stack((xs,ys)))
        plt.draw()
    plt.pause(0.1)

    if not input_thread.is_alive() and data_queue.empty():
        break

plt.ioff()
plt.show()
