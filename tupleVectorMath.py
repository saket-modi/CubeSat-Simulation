import math

######## VECTOR MATH ########
# Performs vector operations on tuples
class Vector:
    @staticmethod
    def vectorAdd(v1, v2): # v1 + v2
        return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2])

    @staticmethod
    def vectorSub(v1, v2): # v1 - v2
        return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])
    
    @staticmethod
    def vectorScalar(v, k): # v = kv
        return (v[0] * k, v[1] * k, v[2] * k)

    @staticmethod
    def vectorDot(v1, v2): # v1 . v2
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    @staticmethod
    def vectorCross(v1, v2): # v1 x v2
        return (
            v1[1] * v2[2] - v2[1] * v1[2], 
            -1 * (v1[0] * v2[2] - v1[2] * v2[0]),
            v1[0] * v2[1] - v1[1] * v2[0]
        )
    
    @staticmethod
    def vectorMag(v): # |v|
        return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    
    @staticmethod
    def vectorAngle(v1, v2):
        # a . b = |a||b|cos(x)
        mag = Vector.vectorMag(v1) * Vector.vectorMag(v2)
        if mag == 0:
            return 0
        
        cosVal = Vector.vectorDot(v1, v2) / mag

        if cosVal >= 1: return 0
        if cosVal <= -1: return math.pi
        return math.acos(cosVal)
    
    @staticmethod
    def vectorUnit(v): # v / |v|
        mag = Vector.vectorMag(v)
        if mag == 0:
            return (0, 0, 0)
        return (v[0] / mag, v[1] / mag, v[2] / mag)