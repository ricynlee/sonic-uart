# !python3
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import numpy as np
import threading
import sys
import queue

data_queue = queue.Queue() # thread safe

if len(sys.argv)>1 and sys.argv[1]=="-test":
    def read_stdin(): # mock
        import random, time, math
        for _ in range(6000):
            for _ in range(10):
                theta, r = random.uniform(-1, 1), math.sqrt(random.uniform(0, 1))
                x, y = r * math.cos(2*math.pi*theta), r * math.sin(2*math.pi*theta)
                data_queue.put((x,y))
            time.sleep(0.01)
else:
    def read_stdin(): # input thread
        for line in sys.stdin:
            x,y = line.strip().split()
            x,y = float(x), float(y)
            data_queue.put((x,y))

input_thread = threading.Thread(target=read_stdin, daemon=True)
input_thread.start()

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
MAX_POINTS = 32

def update(_):
    updated = False
    while not data_queue.empty():
        points.append(data_queue.get())
        updated = True

    if len(points) > MAX_POINTS:
        points[:] = points[-MAX_POINTS:]

    if updated and points:
        xs, ys = zip(*points)
        sc.set_offsets(np.column_stack((xs, ys)))
    return sc,

animation = ani.FuncAnimation(fig, update, interval=100, blit=True, cache_frame_data=False)

plt.show()
