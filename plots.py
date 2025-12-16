import matplotlib.pyplot as plt
import numpy as np
from tuple_operations import *

# WGS-84 constants
EARTH_SEMI_MAJOR = 6378137.0
EARTH_SEMI_MINOR = 6356752.3142

RENDER_SCALE = 1e-6

def make_sphere(pos, radius, res=30):
    """Return coordinates for a sphere centered at pos with given radius."""
    u, v = np.mgrid[0:2*np.pi:res*1j, 0:np.pi:res*1j]
    x = radius * np.cos(u) * np.sin(v) + pos[0]
    y = radius * np.sin(u) * np.sin(v) + pos[1]
    z = radius * np.cos(v) + pos[2]
    return x, y, z

def plot(satellite, planet):
    plt.ion()  # interactive mode

    fig = plt.figure(figsize=(12, 8))

    # ----------------- 3D ORBIT -----------------
    ax3d = fig.add_subplot(2, 1, 1, projection='3d')
    L = 10 * planet.axesParams[0] * RENDER_SCALE
    ax3d.set_xlim(-L, L)
    ax3d.set_ylim(-L, L)
    ax3d.set_zlim(-L, L)
    ax3d.set_box_aspect([1, 1, 1])
    ax3d.set_title("3D Orbit")

    # Earth reference sphere
    xs, ys, zs = make_sphere((0, 0, 0), EARTH_SEMI_MAJOR * RENDER_SCALE)
    ax3d.plot_surface(xs, ys, zs, color='blue', alpha=0.4)

    # Satellite orbit and marker
    orbit_line, = ax3d.plot([], [], [], color='cyan', lw=1)
    sat_point, = ax3d.plot([], [], [], 'ro', markersize=4)

    # ----------------- 2D STATE PLOTS (SIDE BY SIDE) -----------------
    ax_speed = fig.add_subplot(2, 2, 3)
    ax_speed.set_title("Speed (m/s)")
    ax_speed.set_xlabel("Time (s)")
    ax_speed.set_ylabel("Speed")

    ax_alt = fig.add_subplot(2, 2, 4)
    ax_alt.set_title("Altitude (m)")
    ax_alt.set_xlabel("Time (s)")
    ax_alt.set_ylabel("Altitude")

    times = []
    speeds = []
    altitudes = []

    speed_line, = ax_speed.plot([], [], 'r', label='Speed')
    alt_line, = ax_alt.plot([], [], 'b', label='Altitude')

    ax_speed.legend()
    ax_alt.legend()

    import time
    start_time = time.time()
    dt_plot = 0.01

    while plt.fignum_exists(fig.number):
        # --- update 3D orbit ---
        if len(satellite.orbit) > 0:
            xs = [p[0] * RENDER_SCALE for p in satellite.orbit]
            ys = [p[1] * RENDER_SCALE for p in satellite.orbit]
            zs = [p[2] * RENDER_SCALE for p in satellite.orbit]
            orbit_line.set_data(xs, ys)
            orbit_line.set_3d_properties(zs)

        x, y, z = vectorScalar(satellite.pos, RENDER_SCALE)
        sat_point.set_data([x], [y])
        sat_point.set_3d_properties([z])

        # --- update 2D state ---
        t = time.time() - start_time
        times.append(t)
        speeds.append(vectorMag(satellite.vel))
        altitudes.append(satellite.getAltitude())

        speed_line.set_data(times, speeds)
        alt_line.set_data(times, altitudes)

        ax_speed.relim()
        ax_speed.autoscale_view()
        ax_alt.relim()
        ax_alt.autoscale_view()

        plt.pause(dt_plot)
