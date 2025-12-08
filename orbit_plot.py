import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# from main import MAX_FPS

index = 0
dt = 1000 # Updates the plot every second

X = [[], [], []]
V = [[], [], []]
A = [[], [], []]
T = []

def trim(v):
    if len(v) > 100:
        v.pop(0)

def animate(sat):
    pos = sat.pos
    vel = sat.vel
    acc = sat.acc

    trim(X[0]); trim(X[1]); trim(X[2])
    trim(V[0]); trim(V[1]); trim(V[2])
    trim(A[0]); trim(A[1]); trim(A[2])
    trim(T)

    global index
    index += 1

    X[0].append(pos[0]); X[1].append(pos[1]); X[2].append(pos[2])
    V[0].append(vel[0]); V[1].append(vel[1]); V[2].append(vel[2])
    A[0].append(acc[0]); A[1].append(acc[1]); A[2].append(acc[2])
    T.append(index)

    plt.cla()
    plt.plot(T, X[0], color = "red", label = "x")
    plt.plot(T, X[1], color = "green", label = "y")
    plt.plot(T, X[2], color = "blue", label = "z")


def plot(sat):
    ani = FuncAnimation(
        plt.gcf(),
        lambda frame: animate(sat),
        interval=dt
    )

    plt.style.use('dark_background')
    plt.tight_layout()
    plt.legend(loc='upper right')
    plt.show()