import math
class Rayleigh:
    def getTR(M, gamma):
        return ((gamma + 1) ** 2 * M ** 2) / (1 + gamma * M ** 2) ** 2
    
    def getT0R(M, gamma):
        return ((2 * (gamma + 1) * M ** 2) / (1 + gamma * M ** 2) ** 2) * (1 + (gamma - 1) / 2 * M ** 2)
    
    def getPR(M, gamma):
        return (gamma + 1) / (1 + gamma * M ** 2)
    
    def getP0R(M, gamma):
        return (gamma + 1) / (1 + gamma * M ** 2) * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2)) ** (gamma / (gamma - 1))
    
    def getUR(M, gamma):
        return ((gamma + 1) * M ** 2) / (1 + gamma * M ** 2)
    
    def getMfromTR(TR, gamma):
        a = 2 / (gamma + 1)
        b = (gamma - 1) / 2
        return math.sqrt(1 / (TR * a * b) - 1 / b)