import math

######## MISCELLANEOUS HELPER METHODS ########
def toStr(v):
    return str(v[0]) + str(v[1]) + str(v[2])

def getCoords(v, CENTER):
    return (CENTER[0] + v[0], CENTER[1] - v[1])

# Takes a 3D vector and transforms it to 2D using the formula:
# T(x, y, z) = (x - z/sqrt(2), y - z/sqrt(2))
# This works because the program defines the z-axis at an angle pi/4 from both axes
# it's a pseudo 3D implementation
def transform(v, CENTER):
    x, y, z = v
    return getCoords((x - z/math.sqrt(2), y - z/math.sqrt(2)), CENTER)

######## VECTOR MATH ########
# Performs vector operations on tuples
def vectorAdd(v1, v2): # v1 + v2
    return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2])

def vectorSub(v1, v2): # v1 - v2
    return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])

def vectorScalar(v, k): # v = kv
    return (v[0] * k, v[1] * k, v[2] * k)

def vectorDot(v1, v2): # v1 . v2
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def vectorCross(v1, v2): # v1 x v2
    return (
        v1[1] * v2[2] - v2[1] * v1[2], 
        -1 * (v1[0] * v2[2] - v1[2] * v2[0]),
        v1[0] * v2[1] - v1[1] * v2[0]
    )

def vectorMag(v): # |v|
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

def vectorAngle(v1, v2):
    # a . b = |a||b|cos(x)
    mag = vectorMag(v1) * vectorMag(v2)
    if mag == 0:
        return 0
    
    cosVal = vectorDot(v1, v2) / mag

    if cosVal >= 1: return 0
    if cosVal <= -1: return math.pi
    return math.acos(cosVal)

def vectorUnit(v): # v / |v|
    mag = vectorMag(v)
    if mag == 0:
        return (0, 0, 0)
    return (v[0] / mag, v[1] / mag, v[2] / mag)