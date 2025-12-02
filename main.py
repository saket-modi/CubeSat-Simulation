import pygame
import sys
import math
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
    def __init__(self, status, satellite, J2, DRAG, SOLAR, MAG):
        self.status = status
        self.satellite = satellite

        # define ideal orbit here
        self.keplerOrbit = 1 # placeholder value

        # boolean array defining any anomalies
        self.anomalies = [False, False, False, False]
        if (J2): self.anomalies[MissionAnomalies.J2_EFFECT] = True
        if (DRAG): self.anomalies[MissionAnomalies.DRAG] = True
        if (SOLAR): self.anomalies[MissionAnomalies.SOLAR_RADIATION_PRESSURE] = True
        if (MAG): self.anomalies[MissionAnomalies.MAGNETIC_EFFECTS] = True

    def updateStatus(self, status):
        self.status = status

        # performed predefined functions based on current status
        match status:
            case MissionStatus.IN_ORBIT_AND_MONITORING:
                if not self.satellite.checkRange():
                    self.satellite.toggleMonitoring()
                    self.status = MissionStatus.IN_ORBIT_AND_NOT_MONITORING

            case MissionStatus.IN_ORBIT_AND_NOT_MONITORING:
                if self.satellite.checkRange():
                    self.satellite.toggleMonitoring()
                    self.status = MissionStatus.IN_ORBIT_AND_MONITORING

            case MissionStatus.OUT_OF_ORBIT_AND_MONITORING:
                if not self.satellite.checkRange():
                    self.satellite.toggleMonitoring()
                    self.status = MissionStatus.OUT_OF_ORBIT_AND_NOT_MONITORING
                self.satellite.fixOrbit()

            case MissionStatus.OUT_OF_ORBIT_AND_NOT_MONITORING:
                if self.satellite.checkRange():
                    self.satellite.toggleMonitoring()
                    self.status = MissionStatus.OUT_OF_ORBIT_AND_MONITORING
                self.satellite.fixOrbit()

######## RIGID BODY IMPLEMENTATION ########
class Dynamics:
    def __init__(self, pos = (0, 0), vel = (0, 0), acc = (0, 0)):
        self.pos = pos
        self.vel = vel
        self.acc = acc

class Entity:
    def __init__ (self, mass, type, sprite):
        self.mass = mass
        self.type = type
        self.sprite = sprite

class Satellite(Entity, Dynamics):
    def __init__(self, mass, type, sprite, pos, vel, acc):
        self.mass = mass
        self.type = type
        self.sprite = sprite
        self.pos = pos
        self.vel = vel
        self.acc = acc

    # position in polar coordinates (r, theta)
    def updatePosition(self, radius, theta):
        self.dyn.pos = (
            radius * math.cos(theta), 
            radius * math.sin(theta)
        )
    
    def updateVelocity(self, speed, theta):
        self.dyn.vel = (
            speed * math.cos(theta),
            speed * math.sin(theta)
        )

    def updateAcceleration(self, rate, theta):
        self.dyn.acc = (
            rate * math.cos(theta),
            rate * math.sin(theta)
        )

    def updateDynamics(self, radius, speed, rate, theta):
        self.updatePosition(radius, theta)
        self.updateVelocity(speed, theta)
        self.updateAcceleration(rate, theta)

#################### SIMULATION LOGIC ####################

pygame.init()

SCREEN_WIDTH, SCREEN_HEIGHT, MAX_FPS = 800, 600, 60

clock = pygame.time.Clock()
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
running = True
dt = 0 # delta time in seconds

while (running):
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit()

    dt = clock.tick(MAX_FPS) / 1000