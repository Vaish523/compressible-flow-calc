import math
class Fanno:
    def getTR(M, gamma):
        return 1 / ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2))
    
    def getPR(M, gamma):
        return 1 / M * 1 / math.sqrt((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2))
    
    def getP0R(M, gamma):
        return 1 / M * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2)) ** ((gamma + 1) / (2 * gamma - 2))
    
    def getUR(M, gamma):
        return M / math.sqrt((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2))
    
    def getFric(M, gamma):
        return ((1 - gamma * M ** 2) / (gamma * M ** 2)) + ((gamma + 1) / (2 * gamma)) * math.log(M ** 2 / ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M ** 2)))
    
    def getMfromTR(TR, gamma):
        a = 2 / (gamma + 1)
        b = (gamma - 1) / 2
        return math.sqrt(1 / (TR * a * b) - 1 / b)
    
    def getMfromPR(PR, gamma):
        M_guess = 0.001
        tol = 0.0001
        while M_guess < 25:
            PR_test = 1 / M_guess * 1 / math.sqrt((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * M_guess ** 2))
            diff = PR - PR_test
            if diff > 0:
                if diff < tol:
                    return M_guess
                else:
                    M_guess -= 0.00001
            if diff < 0:
                if abs(diff) < tol:
                    return M_guess
                else:
                    M_guess += 0.00001
        return M_guess
    
    def getMfromP0R(P0R, gamma):
        Msub_guess = 0.001
        Msup_guess = 1.001
        tol = 0.0001
        while Msub_guess < 1:
            P0Rsub_test = 1 / Msub_guess * ((2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * Msub_guess ** 2)) ** ((gamma + 1) / (2 * gamma - 2))
            diff_sub = P0R - P0Rsub_test
            if diff_sub > 0:
                if diff_sub < tol:
                    return Msub_guess
                else:
                    Msub_guess -= 0.00001
            if diff_sub < 0:
                if abs(diff_sub) < tol:
                    return Msub_guess
                else:
                    Msub_guess += 0.00001
        return Msub_guess