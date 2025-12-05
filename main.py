import pygame
import sys
import math
import random
from tupleVectorMath import Vector
from enum import Enum

#################### MISSION LOGIC ####################

######## GLOBAL CONSTANTS ########

#### PHYSICS ####
G = 6.6743 * 1e-11
M_EARTH = 5.97219 * 1e24
M_SUN = 1.989 * 1e30
M_MOON = 7.34767309 * 1e22

# Data obtained from WGS-84
EARTH_SEMI_MAJOR = 6378137.0
EARTH_SEMI_MINOR = 6356752.3142

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
                
                rVec = Vector.vectorSub(body.pos, target.pos)
                rMag = Vector.vectorMag(rVec)

                # a body won't apply gravitational force on itself
                # F = Gm1m2 / r^2
                F = 0 if rMag == 0 else G * target.mass * body.mass / (rMag ** 2)

                # unit vector originating from target and going towards the body
                vec = Vector.vectorScalar(Vector.vectorUnit(rVec), F)
                fTarget = Vector.vectorAdd(fTarget, vec)

            target.force = fTarget
            gEntities[i] = target

    def addForces(self, f):
        """ALWAYS TO BE RUN AFTER SIMGRAVITATION() SINCE THE FUNCTION RESETS FORCES TO ZERO"""
        self.force = Vector.vectorAdd(self.force, f)

    def updateDynamics(self, dt):
        if self.mass == 0:
            return

        f = self.force
        self.acc = Vector.vectorScalar(f, 1 / self.mass)

        self.pos = Vector.vectorAdd( # s = ut + 1/2*at^2
            self.pos,
            Vector.vectorAdd(
                Vector.vectorScalar(self.vel, dt),
                Vector.vectorScalar(self.acc, 0.5 * (dt ** 2))
            )
        )

        self.vel = Vector.vectorAdd(self.vel, Vector.vectorScalar(self.acc, dt)) # v = u + at
        
        x, y, z = self.pos

        # If satellite has orbit list, record point
        if hasattr(self, "orbit"):
            self.orbit.append((x, y, z, self.orbitColor))

    def draw(self, screen):
        pygame.draw.circle(screen, "green", self.toScreen(), 3)

    def toScreen(self):
        x, y, z = self.pos
        sx = int(CENTER[0] + x * SCALE)
        sy = int(CENTER[1] - y * SCALE)
        return (sx, sy)

# Cosmic bodies such as Earth, the Moon, asteroids, etc.
class cosmicBody(Entity):
    def __init__(self, mass, center, axesParams, vel = (0.0, 0.0, 0.0), acc = (0.0, 0.0, 0.0), force = (0.0, 0.0, 0.0)):
        Entity.__init__(self, mass, center, vel, acc, force)
        self.center = center # x, y, z
        self.axesParams = axesParams # a, b, c

    def getLaunchLocation(self):
        """
        TBD
        """
        pass

    # return satellite parameters to be displayed on the screen
    def getString(self):
        return (
            "Position: " + self.satellite.pos + 
            "\nVelocity: " + self.satellite.vel +
            "\nAcceleration: " + self.satellite.acc +
            "\nMonitoring: " + ("Yes" if self.satellite.monitoring else "No") +
            "\nStatus: " + self.status
        )

    def draw(self, screen):
        x, y, z = self.center
        a, b, c = self.axesParams
        sx = int(CENTER[0] + x * SCALE)
        sy = int(CENTER[1] - y * SCALE)

        print(a, b)
        referenceRect = pygame.Rect(sx, sy, int(a * SCALE / 2), int(b * SCALE / 2))
        pygame.draw.ellipse(screen, "blue", referenceRect)

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
            r = self.primary.axesParams[0] + targetAlt
            self.pos = (r, 0.0, 0.0) # a + altitude

            # demo: Initial tangential velocity for circular orbit
            v_mag = math.sqrt(G * self.primary.mass / r)
            self.vel = (0, v_mag, 0)

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
        Ux, Uy, Uz = Vector.vectorSub(self.target.pos, self.pos)
        cx, cy, cz = self.primary.center
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
        return Vector.vectorMag(Vector.vectorSub(
            self.pos,
            self.primary.center
        )) - (self.primary.axesParams[0] + self.primary.axesParams[1])/2 

    # Draw satellite and orbit on the screen
    def draw(self, screen):
        # draw orbit
        for ox, oy, oz, color in self.orbit:
            sx = int(CENTER[0] + ox * SCALE)
            sy = int(CENTER[1] - oy * SCALE)
            pygame.draw.circle(screen, color, (sx, sy), 2)

        # draw satellite
        x, y, z = self.pos
        sx = int(CENTER[0] + x * SCALE)
        sy = int(CENTER[1] - y * SCALE)
        pygame.draw.circle(
            screen,
            (255,255,180) if self.monitoring else (200,80,80),
            (sx, sy),
            6
        )

        # draw monitoring line
        if self.checkRange:
            sat_screen = self.toScreen()  # satellite screen coordinates
            tx, ty, tz = self.target.pos
            target_screen = (int(CENTER[0] + tx * SCALE), int(CENTER[1] - ty * SCALE))
            pygame.draw.line(screen, "red", sat_screen, target_screen, 1)


        # status text
        # r = math.sqrt(x*x + y*y + z*z)
        # txt = font.render(
        #     f"Status: {manager.status.name}   r={r:.1f}",
        #     True, (220,220,220)
        # )
        # screen.blit(txt, (10, 10))


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
    center = (0, 0, 0), 
    axesParams = (EARTH_SEMI_MAJOR, EARTH_SEMI_MINOR, (EARTH_SEMI_MAJOR + EARTH_SEMI_MINOR)/2)
)

targetStation = Entity(
    mass = 0, 
    pos = (EARTH_SEMI_MAJOR, 0, 0)
)

# pos = planetOne.getLaunchLocation()
CubeSat_One = Satellite(
    mass = 1, 
    target = targetStation,
    primaryBody = planetOne,
    targetAlt = 400_000 # 400 kilometres
)

CubeSat_One_MissionManager = MissionManager(
    satellite = CubeSat_One
)

# game loop and text 
running = True
font = pygame.font.SysFont(None, 20)
t = 0

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    dt = clock.tick(MAX_FPS)
    t += dt

    # Fill the screen with a dark background
    screen.fill((0, 0, 10))

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