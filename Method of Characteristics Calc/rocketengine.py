# Author: Vaishnav Srivaths (vsrivat@purdue.edu)
class rocketEngine:
    def __init__(self, P_c, T_c, P_e, gamma, R, r_t, numPoints, numLines):
        self.P_c = P_c
        self.T_c = T_c
        self.P_e = P_e
        self.gamma = gamma
        self.R = R
        self.r_t = r_t
        self.numPoints = numPoints
        self.numLines = numLines

    def getMachExit(self, P_c, P_e, gamma):
        M_e = (((P_e / P_c) ** ((1 - gamma) / gamma) - 1) * 2 / (gamma - 1)) ** 0.5
        return M_e
    
    def getTempExit(self, T_c, M_e, gamma):
        T_e = T_c / (1 + (gamma - 1) / 2 * M_e ** 2)
        return T_e
    
    def getSoundSpeedExit(self, gamma, R, T_e):
        a_e = (gamma * R * T_e) ** 0.5
        return a_e
    
    def getVelocityExit(self, a_e, M_e):
        v_e = a_e * M_e
        return v_e
    
    def getAreaRatio(self, M_e, gamma):
        areaRatio = (((gamma + 1) / 2) / (1 + (gamma - 1) / 2 * M_e ** 2)) ** ((gamma + 1) / (2 - 2 * gamma)) / M_e
        return areaRatio
    
    def getExitRadius(self, areaRatio, r_t):
        r_e = r_t * areaRatio ** 0.5
        return r_e
    
    