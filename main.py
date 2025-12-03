import pygame
import sys
import math
import random
from tupleVectorMath import Vector
from enum import Enum

#################### MISSION LOGIC ####################

######## GLOBAL CONSTANTS ########

# PHYSICS
G = 6.6743 * 1e-11
M_EARTH = 5.97219 * 1e24
M_SUN = 1.989 * 1e30
M_MOON = 7.34767309 * 1e22

# SIMULATION
SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
MAX_FPS = 60
CENTER = (SCREEN_WIDTH // 2, SCREEN_HEIGHT // 2)

######## ENUMERABLE VALUES ########
class MissionStatus(Enum):
    LAUNCHING = -1
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

######## GLOBAL STORAGE ########
gEntities = [] # stores all entities
gMissions = [] # all satellite missions
    
######## CONTROL LOGIC ########
class MissionManager:
    def __init__(self, status, mass, target, J2=False, DRAG=False, SOLAR=False, MAG=False):
        self.status = status
        self.satellite = Satellite(mass, target)

        self.keplerOrbit = 200.0  # pixels (display units)

        self.anomalies = [False] * len(MissionAnomalies)
        if J2:    self.anomalies[MissionAnomalies.J2_EFFECT.value] = True
        if DRAG:  self.anomalies[MissionAnomalies.DRAG.value] = True
        if SOLAR: self.anomalies[MissionAnomalies.SOLAR_RADIATION_PRESSURE.value] = True
        if MAG:   self.anomalies[MissionAnomalies.MAGNETIC_EFFECTS.value] = True

    def getStatus(self):
        return self.status

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
class Entity:
    def __init__(self, mass, pos=(0.0, 0.0, 0.0), vel=(0.0, 0.0, 0.0), acc=(0.0, 0.0, 0.0), force=(0.0, 0.0, 0.0)):
        self.mass = mass
        gEntities.append(self)
        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.force = force

    # initialize destructor to remove object from entities global array
    def __del__(self):
        gEntities.remove(self)

    # DYNAMICS

    def simGravitation():
        """MUST BE RUN EVERY FRAME TO ACCOUNT FOR POSITIONAL CHANGES"""
        for i in range(len(gEntities)):
            target = gEntities[i]
            fTarget = (0, 0, 0)

            for j in range(len(gEntities)):
                body = gEntities[j]
                
                if i == j:
                    continue
                
                rVec = Vector.vectorSub(body.pos, target.pos)

                # a body won't apply gravitational force on itself
                # F = Gm1m2 / r^2
                F = G * target.mass * body.mass / Vector.vectorDot(rVec, rVec)

                # unit vector originating from target and going towards the body
                unitVec = Vector.vectorUnit(rVec)
                resolvedVec = Vector.vectorResolve(Vector.vectorScalar(unitVec, F))
                fTarget = Vector.vectorAdd(fTarget, resolvedVec)

            target.force = fTarget
            gEntities[i] = target

    def addForces(self, f):
        """ALWAYS TO BE RUN AFTER SIMGRAVITATION() SINCE THE FUNCTION RESETS FORCES TO ZERO"""
        self.forces = Vector.vectorAdd(self.forces, f)

    def updateDynamics(self, dt):
        f = self.force
        self.acc = Vector.vectorScalar(f, 1 / self.mass)

        self.pos = Vector.vectorAdd( # s = ut + 1/2*at^2
            Vector.vectorScalar(self.vel, dt),
            Vector.vectorScalar(self.acc, 0.5 * (dt ** 2))
        )

        self.vel = Vector.vectorAdd(self.vel, Vector.vectorScalar(self.acc, dt)) # v = u + at
        
        x, y, z = self.pos

        # If satellite has orbit list, record point
        if hasattr(self, "orbit"):
            self.orbit.append((x, y, z, self.orbitColor))

class Satellite(Entity):
    def __init__(self, mass, target, pos=(0.0, 0.0, 0.0), vel=(0.0, 0.0, 0.0), acc=(0.0, 0.0, 0.0), force=(0.0, 0.0, 0.0)):
        Entity.__init__(self, mass, pos, vel, acc, force)
        self.monitoring = False
        self.target = target # target is an entity
        self.orbit = []
        self.orbitColor = (
                random.randint(50, 255),
                random.randint(50, 255),
                random.randint(50, 255)
        )

    # Check if the target is within the visible range of the satellite
    def checkRange(self, gravitationBody):
        """
        parametric form of line => L(t) = P + tU
        P = current position of satellite
        U = direction vector
        t = parameter

        Equation of ellipsoid = ((x - cx) / a)^2 + ((y - cy) / b)^2 + ((z - cz) / c)^2 = 1
        substitute x = L(tx), y = L(ty), z = L(tz)
        if no solutions to t exist (D < 0) => return true
        if D = 0 => edge case (tangential) OR D > 0 => return false

        the obtained coefficients of the quadratic expression are:
        A = (Ux*Ux)/(a*a) + (Uy*Uy)/(b*b) + (Uz*Uz)/(c*c)
        B = 2 * ( (dx*Ux)/(a*a) + (dy*Uy)/(b*b) + (dz*Uz)/(c*c) )
        C = (dx*dx)/(a*a) + (dy*dy)/(b*b) + (dz*dz)/(c*c) - 1

        Where:
        dx = Px - cx
        dy = Py - cy
        dz = Pz - cz
        """
        Px, Py, Pz = self.pos
        Ux, Uy, Uz = Vector.vectorSub(self.target.pos, self.pos)
        cx, cy, cz = gravitationBody.centre
        dx, dy, dz = Px - cx, Py - cy, Pz - cz
        a, b, c = gravitationBody.axesParams

        A = Ux ** 2 / a ** 2 + Uy ** 2 / b ** 2 + Uz ** 2 / c ** 2
        B = 2 * ((dx * Ux)/a**2 + (dy * Uy)/b**2 + (dz * Uz)/c**2)
        C = dx ** 2 / a ** 2 + dy ** 2 / b ** 2 + dz ** 2 / c ** 2 - 1

        toCheck = B ** 2 - 4 * A * C
        return toCheck < 0

    # Toggle satellite monitoring of the target
    def toggleMonitoring(self):
        self.monitoring = not self.monitoring
        print("toggleMonitoring ->", self.monitoring)

    # TBD : Fix satellite's orbit by firing thrusters
    def fixOrbit(self, desired_radius):
        x, y, z = self.pos
        r = math.sqrt(x*x + y*y + z*z)
        if r == 0: return

        ux, uy, uz = x/r, y/r, z/r
        correction = (desired_radius - r) * 0.05
        self.pos = (x + ux*correction, y + uy*correction, z + uz*correction)
        self.vel = (0.0, 0.0, 0.0)

    # Draw satellite and orbit on the screen
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

        # draw monitoring line
        if (self.monitoring):
            x, y, z = self.pos
            h, k, l = self.target.pos
            pygame.draw.line("red", (x, y), (h, k), 3)

        # status text
        r = math.sqrt(x*x + y*y + z*z)
        txt = font.render(
            f"Status: {manager.status.name}   r={r:.1f}",
            True, (220,220,220)
        )
        screen.blit(txt, (10, 10))


#################### SIMULATION LOGIC ####################

pygame.init()

clock = pygame.time.Clock()
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption("Satellite Demo (with Orbit Trail)")

# Create satellites and assign their corresponding managers
sat = Satellite(mass=1.0)
manager = MissionManager(
    status=MissionStatus.LAUNCHING,
    satellite=sat
)

# game loop and text 
running = True
t = 0.0
font = pygame.font.SysFont(None, 20)

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    dt = clock.tick(MAX_FPS) / 1000.0
    t += dt

    sat.updateDynamics(radius=r, phi=phi, theta=theta, v=(0,0,0), f=(0,0,0))
    manager.updateStatus()

    screen.fill((0, 0, 10))
    pygame.draw.circle(screen, (20, 60, 20), CENTER, 50)  # Earth

    sat.draw(screen)
    pygame.display.flip()

pygame.quit()
sys.exit()