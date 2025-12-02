import pygame
import sys
import math
import random
from enum import Enum

#################### MISSION LOGIC ####################

######## ENUMERABLE VALUES ########
class MissionStatus(Enum):
    IN_ORBIT_AND_MONITORING = 0
    IN_ORBIT_AND_NOT_MONITORING = 1
    IN_ORBIT_AND_OUT_OF_RANGE = 2
    OUT_OF_ORBIT_AND_MONITORING = 3
    OUT_OF_ORBIT_AND_NOT_MONITORING = 4
    OUT_OF_ORBIT_AND_OUT_OF_RANGE = 5

class MissionAnomalies(Enum):
    J2_EFFECT = 0
    DRAG = 1
    SOLAR_RADIATION_PRESSURE = 2
    MAGNETIC_EFFECTS = 3

######## CONTROL LOGIC ########
class MissionManager:
    def __init__(self, status, satellite, J2=False, DRAG=False, SOLAR=False, MAG=False):
        self.status = status
        self.satellite = satellite

        self.keplerOrbit = 200.0  # pixels (display units)

        self.anomalies = [False] * len(MissionAnomalies)
        if J2:    self.anomalies[MissionAnomalies.J2_EFFECT.value] = True
        if DRAG:  self.anomalies[MissionAnomalies.DRAG.value] = True
        if SOLAR: self.anomalies[MissionAnomalies.SOLAR_RADIATION_PRESSURE.value] = True
        if MAG:   self.anomalies[MissionAnomalies.MAGNETIC_EFFECTS.value] = True

    def updateStatus(self):
        if self.status == MissionStatus.IN_ORBIT_AND_MONITORING:
            if not self.satellite.checkRange():
                self.satellite.toggleMonitoring()
                self.status = MissionStatus.IN_ORBIT_AND_NOT_MONITORING

        elif self.status == MissionStatus.IN_ORBIT_AND_NOT_MONITORING:
            if self.satellite.checkRange():
                self.satellite.toggleMonitoring()
                self.status = MissionStatus.IN_ORBIT_AND_MONITORING

        elif self.status == MissionStatus.OUT_OF_ORBIT_AND_MONITORING:
            if not self.satellite.checkRange():
                self.satellite.toggleMonitoring()
                self.status = MissionStatus.OUT_OF_ORBIT_AND_NOT_MONITORING
            self.satellite.fixOrbit(self.keplerOrbit)

        elif self.status == MissionStatus.OUT_OF_ORBIT_AND_NOT_MONITORING:
            if self.satellite.checkRange():
                self.satellite.toggleMonitoring()
                self.status = MissionStatus.OUT_OF_ORBIT_AND_MONITORING
            self.satellite.fixOrbit(self.keplerOrbit)


######## RIGID BODY IMPLEMENTATION ########
class Dynamics:
    def __init__(self, pos=(0.0, 0.0, 0.0), vel=(0.0, 0.0, 0.0), acc=(0.0, 0.0, 0.0)):
        self.pos = pos
        self.vel = vel
        self.acc = acc

    def updatePosition(self, radius, theta, phi):
        x = radius * math.cos(phi) * math.cos(theta)
        y = radius * math.cos(phi) * math.sin(theta)
        z = radius * math.sin(phi)
        self.pos = (x, y, z)

        # If satellite has orbit list, record point
        if hasattr(self, "orbit"):
            self.orbit.append((x, y, z, self.orbitColor))

    def updateVelocity(self, speed, theta, phi):
        vx = speed * math.cos(phi) * math.cos(theta)
        vy = speed * math.cos(phi) * math.sin(theta)
        vz = speed * math.sin(phi)
        self.vel = (vx, vy, vz)

    def updateAcceleration(self, rate, theta, phi):
        ax = rate * math.cos(phi) * math.cos(theta)
        ay = rate * math.cos(phi) * math.sin(theta)
        az = rate * math.sin(phi)
        self.acc = (ax, ay, az)


class Entity:
    def __init__(self, mass, type, sprite=None):
        self.mass = mass
        self.type = type
        self.sprite = sprite


class Satellite(Entity, Dynamics):
    def __init__(self, mass, type, sprite, pos=(0.0, 0.0, 0.0), vel=(0.0, 0.0, 0.0), acc=(0.0, 0.0, 0.0)):
        Entity.__init__(self, mass, type, sprite)
        Dynamics.__init__(self, pos, vel, acc)
        self.monitoring = False
        self.orbit = []
        self.orbitColor = (255, 255, 255)  # white

    # --------------------------------------------------
    # ORBIT LOGIC
    # --------------------------------------------------
    def checkOrbitClosure(self, threshold=3.0):
        x, y, z = self.pos
        for px, py, pz, _ in self.orbit[:-30]: # _ is the color
            dx = x - px
            dy = y - py
            dz = z - pz
            if dx * dx + dy * dy + dz * dz < threshold * threshold:
                return True
        return False

    def updateDynamics(self, radius, speed, rate, theta, phi):
        self.updatePosition(radius, theta, phi)
        self.updateVelocity(speed, theta, phi)
        self.updateAcceleration(rate, theta, phi)

        if self.checkOrbitClosure():
            self.orbitColor = (
                random.randint(50, 255),
                random.randint(50, 255),
                random.randint(50, 255)
            )

    # --------------------------------------------------
    # DRAW
    # --------------------------------------------------
    def draw(self, screen):
        # draw orbit
        for ox, oy, oz, color in self.orbit:
            sx = int(CENTER[0] + ox)
            sy = int(CENTER[1] - oy)
            pygame.draw.circle(screen, color, (sx, sy), 2)

        # draw satellite
        x, y, z = self.pos
        sx = int(CENTER[0] + x)
        sy = int(CENTER[1] - y)
        pygame.draw.circle(
            screen,
            (255,255,180) if self.monitoring else (200,80,80),
            (sx, sy),
            6
        )

        # status text
        r = math.sqrt(x*x + y*y + z*z)
        txt = font.render(
            f"Status: {manager.status.name}   r={r:.1f}",
            True, (220,220,220)
        )
        screen.blit(txt, (10, 10))

    # --------------------------------------------------
    # RANGE + CONTROL
    # --------------------------------------------------
    def checkRange(self):
        x, y, z = self.pos
        r = math.sqrt(x*x + y*y + z*z)
        return 180.0 <= r <= 220.0

    def toggleMonitoring(self):
        self.monitoring = not self.monitoring
        print("toggleMonitoring ->", self.monitoring)

    def fixOrbit(self, desired_radius):
        x, y, z = self.pos
        r = math.sqrt(x*x + y*y + z*z)
        if r == 0: return

        ux, uy, uz = x/r, y/r, z/r
        correction = (desired_radius - r) * 0.05
        self.pos = (x + ux*correction, y + uy*correction, z + uz*correction)
        self.vel = (0.0, 0.0, 0.0)


#################### SIMULATION LOGIC ####################

pygame.init()

SCREEN_WIDTH, SCREEN_HEIGHT, MAX_FPS = 800, 600, 60
CENTER = (SCREEN_WIDTH // 2, SCREEN_HEIGHT // 2)

clock = pygame.time.Clock()
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption("Satellite Demo (with Orbit Trail)")

# Create satellite and manager
sat = Satellite(mass=1.0, type="CubeSat", sprite=None)
manager = MissionManager(
    status=MissionStatus.IN_ORBIT_AND_MONITORING,
    satellite=sat
)

# initial spherical coordinates
r = 200.0
theta = 0.0
phi = 0.0
angular_speed = 0.6

running = True
t = 0.0

font = pygame.font.SysFont(None, 20)

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    dt = clock.tick(MAX_FPS) / 1000.0
    t += dt

    # circular motion
    theta += angular_speed * dt
    r += 5.0 * math.sin(0.2 * t) * dt  # small drift

    sat.updateDynamics(radius=r, speed=0, rate=0, theta=theta, phi=phi)
    manager.updateStatus()

    screen.fill((0, 0, 10))
    pygame.draw.circle(screen, (20, 60, 20), CENTER, 50)  # Earth

    sat.draw(screen)
    pygame.display.flip()

pygame.quit()
sys.exit()