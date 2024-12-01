import math
class Isentropic:
    def getP_P0(M, gamma):
        p_p0 = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (1 - gamma))
        return p_p0
    
    def getT_T0(M, gamma):
        t_t0 = (1 + (gamma - 1) / 2 * M ** 2) ** (-1)
        return t_t0
    
    def getRho_Rho0(M, gamma):
        rho_rho0 = (1 + (gamma - 1) / 2 * M ** 2) ** (1 / (1 - gamma))
        return rho_rho0
    
    def getAlpha(M):
        if M >= 1:
            alpha = math.asin(1 / M) * 180 / math.pi
            return alpha
        else:
            return "N/A"
    
    def getNu(M, gamma):
        if M >= 1:
            nu = 1 * (math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1)))
            return nu
        else:
            return "N/A"
    
    def getA_Astar(M, gamma):
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        return A_Astar

    def getMfromPR(pr, gamma):
        return math.sqrt((2 * (pr ** (1 / gamma - 1) - 1)) / (gamma - 1))
    
    def getMfromTR(tr, gamma):
        return math.sqrt((2 * (1/tr - 1)) / (gamma - 1))
    
    def getMfromDR(dr, gamma):
        return math.sqrt((2 * (dr ** (1 - gamma) - 1)) / (gamma - 1))
    
    def getMfromAlpha(alpha, gamma):
        M_guess = 1.00001
        tol = 0.0001
        while M_guess < 25:
            alpha_test = math.asin(1 / M_guess) * 180 / math.pi
            diff = alpha - alpha_test
            if diff < tol:
                return M_guess
            else:
                M_guess += 0.00001
        return M_guess
    
    def getMfromNu(nu, gamma):
        M_guess = 1.00001
        tol = 0.0001
        while M_guess < 25:
            nu_test = abs(180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M_guess ** 2 - 1))) + math.atan(math.sqrt(M_guess ** 2 - 1))))
            diff = nu - nu_test
            if diff < tol:
                return M_guess
            else:
                M_guess += 0.00001
        return M_guess