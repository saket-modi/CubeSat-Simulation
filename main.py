import pygame
import sys
import math
import random
import threading
import plots as plts
from tuple_operations import *
from enum import Enum

#################### MISSION LOGIC ####################

######## GLOBAL CONSTANTS ########

#### PHYSICS ####
G = 6.6743 * 1e-11
M_EARTH = 5.97219 * 1e24
M_SUN = 1.989 * 1e30
M_MOON = 7.34767309 * 1e22
INTEGRATOR = 'EULER'

# Data obtained from WGS-84
EARTH_SEMI_MAJOR = 6378137.0
EARTH_SEMI_MINOR = 6356752.3142
MOON_SEMI_MAJOR = 1737e3
MOON_SEMI_MINOR = 1737e3

#### SIMULATION ####
SCALE = 1e-5 # Scale physics to pixels
SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
MAX_FPS = 60
CENTER = (SCREEN_WIDTH // 2, SCREEN_HEIGHT // 2)

######## ENUMERABLE VALUES ########
class MissionStatus(Enum):
    IN_ORBIT = 0
    OUT_OF_ORBIT = 1

class MissionPerturbations(Enum):
    J2_EFFECT = 0
    DRAG = 1
    SOLAR_RADIATION_PRESSURE = 2
    MAGNETIC_EFFECTS = 3

######## GLOBAL STORAGE ########
gEntities = [] # stores all entities
gMissions = [] # all satellite missions

######## CONTROL LOGIC ########
class MissionManager:
    def __init__(self, satellite, J2=False, DRAG=False, SOLAR=False, MAG=False):
        self.satellite = satellite
        gMissions.append(self)

        self.anomalies = [False] * len(MissionPerturbations)
        if J2:    self.anomalies[MissionPerturbations.J2_EFFECT.value] = True
        if DRAG:  self.anomalies[MissionPerturbations.DRAG.value] = True
        if SOLAR: self.anomalies[MissionPerturbations.SOLAR_RADIATION_PRESSURE.value] = True
        if MAG:   self.anomalies[MissionPerturbations.MAGNETIC_EFFECTS.value] = True

    def __del__(self):
        gMissions.remove(self)

    ## RUNS EVERY FRAME
    def missionUpdate(self, dt):
        pass

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
    @staticmethod
    def simGravitation():
        """
        MUST BE RUN EVERY FRAME TO ACCOUNT FOR POSITIONAL CHANGES
        """
        for i in range(len(gEntities)):
            target = gEntities[i]
            fTarget = (0, 0, 0)

            for j in range(len(gEntities)):
                body = gEntities[j]
                
                if i == j:
                    continue
                
                rVec = vectorSub(body.pos, target.pos)
                rMag = vectorMag(rVec)

                # a body won't apply gravitational force on itself
                # F = Gm1m2 / r^2
                F = 0 if rMag == 0 else G * target.mass * body.mass / (rMag ** 2)

                # unit vector originating from target and going towards the body
                vec = vectorScalar(vectorUnit(rVec), F)
                fTarget = vectorAdd(fTarget, vec)

            target.force = fTarget
            gEntities[i] = target

    def addForces(self, f):
        """ALWAYS TO BE RUN AFTER SIMGRAVITATION() SINCE THE FUNCTION RESETS FORCES TO ZERO"""
        self.force = vectorAdd(self.force, f)

    # Euler implementation
    def updateDynamics(self, dt):
        if self.mass == 0:
            return

        f = self.force
        self.acc = vectorScalar(f, 1 / self.mass)

        if INTEGRATOR == 'EULER':
            self.pos = vectorAdd( # s = ut + 1/2*at^2
                self.pos,
                vectorAdd(
                    vectorScalar(self.vel, dt),
                    vectorScalar(self.acc, 0.5 * (dt ** 2))
                )
            )

            self.vel = vectorAdd(self.vel, vectorScalar(self.acc, dt)) # v = u + at
            
            x, y, z = self.pos

            # If satellite has orbit list, record point
            if hasattr(self, "orbit"):
                self.orbit.append((x, y, z, self.orbitColor))

    def draw(self, screen):
        # pygame.draw.circle(screen, "green", transform(), 3)
        pass

    def toScreen(self):
        x, y, z = self.pos
        sx = int(CENTER[0] + x * SCALE)
        sy = int(CENTER[1] - y * SCALE)
        return (sx, sy)

# Cosmic bodies such as Earth, the Moon, asteroids, etc.
class cosmicBody(Entity):
    def __init__(self, mass, pos, axesParams, vel = (0.0, 0.0, 0.0), acc = (0.0, 0.0, 0.0), force = (0.0, 0.0, 0.0)):
        Entity.__init__(self, mass, pos, vel, acc, force)
        self.pos = pos # x, y, z
        self.axesParams = axesParams # a, b, c

    def getLaunchLocation(self):
        """
        TBD
        """
        pass

    def draw(self, screen):
        x, y, z = self.pos
        a, b, c = self.axesParams

        width_px = int(2 * a * SCALE)
        height_px = int(2 * b * SCALE)

        screen_center = transform(vectorScalar((x, y, z), SCALE), CENTER)
        top_left = (int(screen_center[0] - width_px/2), int(screen_center[1] - height_px/2))

        referenceRect = pygame.Rect(top_left[0], top_left[1], max(1,width_px), max(1,height_px))
        pygame.draw.ellipse(screen, "yellow", referenceRect, width = 1)

class Satellite(Entity):
    def __init__(self, mass, target, primaryBody, targetAlt, pos=(0.0, 0.0, 0.0), vel=(0.0, 0.0, 0.0), acc=(0.0, 0.0, 0.0), force=(0.0, 0.0, 0.0)):
        Entity.__init__(self, mass, pos, vel, acc, force)
        self.monitoring = False
        self.target = target # target is an entity
        self.orbit = []
        self.orbitColor = (
                random.randint(50, 255),
                random.randint(50, 255),
                random.randint(50, 255)
        )

        # The body around which the satellite is revolving
        self.primary = primaryBody

        # The altitude at which the satellite will revolve with an error correction of +-5%
        self.targetAlt = targetAlt

        # demonstration
        if self.pos == (0.0, 0.0, 0.0):
            r = self.primary.axesParams[0]/2 + targetAlt
            self.pos = (r, 0.0, 0.0) # a/2 + altitude

            # demo: Initial tangential velocity for circular orbit
            v_mag = math.sqrt(1.5 * G * self.primary.mass / r)
            self.vel = (0, v_mag, 0)

        def rotateX(v, angle):
            x, y, z = v
            c = math.cos(angle)
            s = math.sin(angle)
            return (x, y*c - z*s, y*s + z*c)

        # tilt orbit by 30 degrees
        tilt = math.radians(30)
        self.pos = rotateX(self.pos, tilt)
        self.vel = rotateX(self.vel, tilt)

    # Check if the target is within the visible range of the satellite
    def checkRange(self):
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
        Ux, Uy, Uz = vectorSub(self.target.pos, self.pos)
        cx, cy, cz = self.primary.pos
        dx, dy, dz = Px - cx, Py - cy, Pz - cz
        a, b, c = self.primary.axesParams

        A = Ux ** 2 / a ** 2 + Uy ** 2 / b ** 2 + Uz ** 2 / c ** 2
        B = 2 * ((dx * Ux)/a**2 + (dy * Uy)/b**2 + (dz * Uz)/c**2)
        C = dx ** 2 / a ** 2 + dy ** 2 / b ** 2 + dz ** 2 / c ** 2 - 1

        toCheck = B ** 2 - 4 * A * C # discriminant D = b^2 - 4ac
        self.monitoring = toCheck < 0
        return toCheck < 0

    # Fix satellite's orbit by firing thrusters
    def fixOrbit(self):
        a, b, c = self.primary.axesParams
        current = self.getAltitude() # a, b >> altitude so this works fine
        acceptable_error = 0.05

        # TBD: Orbital correction logic
        if current < 0.95 * self.targetAlt:
            pass

        if current > 1.05 * self.targetAlt:
            pass 

    # get the current altitude of the satellite
    def getAltitude(self):
        # (a + b) / 2 roughly gives the radii since a, b >> altitude
        return vectorMag(vectorSub(
            self.pos,
            self.primary.pos
        )) - (self.primary.axesParams[0] + self.primary.axesParams[1])/2
    
    # return satellite parameters to be displayed on the screen
    def getString(self):
        # "Position: " + Misc.toStr(self.pos) + 
        # "Velocity: " + Misc.toStr(self.vel) +
        # "Acceleration: " + Misc.toStr(self.acc) +
        # "Monitoring: " + ("Yes" if self.monitoring else "No") +
        return (
            "Altitude: " + str(round(self.getAltitude()/1000, ndigits = 2)) + " km"
        )

    # Draw satellite and orbit on the screen
    def draw(self, screen):
        # draw orbit
        for ox, oy, oz, color in self.orbit:
            coords = transform(vectorScalar((ox, oy, oz), SCALE), CENTER)
            pygame.draw.circle(screen, color, coords, 1)

        # draw satellite
        x, y, z = self.pos
        coords = transform(vectorScalar((x, y, z), SCALE), CENTER)
        pygame.draw.circle(
            screen,
            (255,255,180) if self.monitoring else (200,80,80),
            coords,
            3
        )

        # draw monitoring line
        if self.checkRange():
            sat_screen = self.toScreen()  # satellite screen coordinates
            tx, ty, tz = self.target.pos
            coords = transform(vectorScalar((tx, ty, tz), SCALE), CENTER)
            pygame.draw.line(screen, "red", sat_screen, coords, 1)

        # status text
        txt = font.render(
            self.getString(),
            True, (220,220,220)
        )
        screen.blit(txt, (10, 10))

#################### SIMULATION LOGIC ####################

pygame.init()

clock = pygame.time.Clock()
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption("CubeSats Simulation")

"""
planetOne here simulates the Earth, which is an oblate spheroid
An oblate spheroid can be represented in 2D as an ellipse

As per WGS-84, the standard parameters for the Earth's ellipsoid shape is (in meters):
a = 6378137.0
b = 6378137.0
c = 6356752.3142

Which means c is 0.3364% smaller than a (or b)
"""
planetOne = cosmicBody(
    mass = M_EARTH, 
    pos = (0, 0, 0), 
    axesParams = (EARTH_SEMI_MAJOR, EARTH_SEMI_MINOR, (EARTH_SEMI_MAJOR + EARTH_SEMI_MINOR)/2)
)

# planetTwo = cosmicBody(
#     mass = M_EARTH * 100,
#     pos = (384400e3, 384400e3, 384400e3),
#     axesParams = (MOON_SEMI_MAJOR, MOON_SEMI_MINOR, (MOON_SEMI_MAJOR + MOON_SEMI_MINOR)/2),
#     vel = (-1024/0.707, -1024/0.707, 0)
# )

targetStation = Entity(
    mass = 0, 
    pos = (EARTH_SEMI_MAJOR, 0, 0)
)

# pos = planetOne.getLaunchLocation()
CubeSat_One = Satellite(
    mass = 1, 
    target = targetStation,
    primaryBody = planetOne,
    targetAlt = 4_000_000 # 4000 kilometres
)

CubeSat_One_MissionManager = MissionManager(
    satellite = CubeSat_One
)

# RUN PLOTS FROM plots.py

# plot orbit and satellite along with cosmic entities (Earth, Moon, etc.)
# plot the position, velocity, and acceleration of the satellite
threading.Thread(
    target=plts.plot,
    args=(CubeSat_One, planetOne),
    daemon=True
).start()

# game loop and text 
running = True
font = pygame.font.SysFont(None, 20)
t = 0

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # This should actually be divided by 1000 to be converted to seconds
    # however a large value of dt is presumed to speed up simulation as a simple "hacky" fix
    dt = clock.tick(MAX_FPS) 
    t += dt

    # Fill the screen with a dark background
    screen.fill((0, 0, 10))

    # Draw the POSITIVE x, y, and z axes. z-axis will be displayed at an angle pi/4 from both x and y
    # Note that pygame's coordinate system instantiates the origin at the TOP LEFT CORNER
    # Basically, add the coordinate to the x-axis, add -1 * coordinate for y-axis
    # So (x, y) in this reference is (0, 0) [CENTER] + (x, -y)

    # x-axis
    pygame.draw.line(screen, "red", CENTER, getCoords((SCREEN_WIDTH // 2, 0), CENTER))
    # y-axis
    pygame.draw.line(screen, "green", CENTER, getCoords((0, SCREEN_HEIGHT // 2), CENTER))
    # z-axis
    pygame.draw.line(screen, "blue", CENTER, getCoords((-SCREEN_WIDTH // 2, -SCREEN_HEIGHT // 2), CENTER))

    # simulates gravitational forces for each entity
    Entity.simGravitation()

    # updating attributes for satellites based on mission status
    CubeSat_One_MissionManager.missionUpdate(dt)

    # updating dynamics for all entities
    for i in range(len(gEntities)):
        gEntities[i].updateDynamics(dt)
        gEntities[i].draw(screen)
    
    pygame.display.flip()

pygame.quit()
sys.exit()