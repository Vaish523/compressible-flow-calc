import math
class NormalShock:
    def getM2(M1, gamma):
        return math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
    
    def getP2_P1(M1, gamma):
        return (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
    
    def getRho2_Rho1(M1, gamma):
        return ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
    
    def getT2_T1(M1, gamma):
        return ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
    
    def getP02_P01(M1, gamma):
        return ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
    
    def getM1fromM2(M2, gamma):
        return math.sqrt( ((gamma - 1) * M2 ** 2 + 2) / (2 * gamma * M2 ** 2 - (gamma - 1)) )
    
    def getM1fromPR(p2_p1, gamma):
        return math.sqrt((p2_p1 * (gamma + 1) + (gamma - 1)) / (2 * gamma))
    
    def getM1fromDR(rho2_rho1, gamma):
        return math.sqrt((-2 * rho2_rho1) / (rho2_rho1 * (gamma - 1) - (gamma + 1)))
    
    def getM1fromTR(t2_t1, gamma):
        M1_guess = 1.00001
        tol = 0.0001
        while M_guess < 25:
            t2_t1_test = ((2 * gamma * M1_guess ** 2 - (gamma - 1)) * ((gamma - 1) * M1_guess ** 2 + 2)) / ((gamma + 1) ** 2 * M1_guess ** 2)
            diff = t2_t1 - t2_t1_test
            if diff < tol:
                return M_guess
            else:
                M_guess += 0.00001
        return M_guess
    
    def getM1fromP0R(p02_p01, gamma):
        M1_guess = 1.00001
        tol = 0.0001
        while M_guess < 25:
            p02_p01_test = ((((gamma + 1) * M1_guess ** 2) / ((gamma - 1) * M1_guess ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1_guess ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
            diff = p02_p01 - p02_p01_test
            if diff < tol:
                return M_guess
            else:
                M_guess += 0.00001
        return M_guess