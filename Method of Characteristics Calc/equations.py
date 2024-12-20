# Author: Vaishnav Srivaths (vsrivat@purdue.edu)
import math
class Equations:
    def getPMAngle(M, gamma):
        nu = 180 / math.pi * (math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) - math.atan(math.sqrt(M ** 2 - 1)))
        return abs(nu)
    
    def getMachFromPM(nu, gamma):
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
    
    def getMachAngle(M):
        mu = 180 / math.pi * math.asin(1/M)
        return mu
    
    def interpolateTheta(theta_start, theta, numLines):
        deltaTheta = (theta - theta_start) / (numLines - 2)
        thetas = [0]
        for index in range(numLines - 2):
            thetas.append(thetas[index] + deltaTheta)
        thetas.append(theta)
        return thetas
    
    def intersection(x_left, y_left, theta_right, x_right, y_right, theta_left):
        slope_right = math.tan(theta_right * math.pi / 180)
        slope_left = math.tan(theta_left * math.pi / 180)
        x = (x_left * slope_right - x_right * slope_left + y_right - y_left) / (slope_right - slope_left)
        y = (slope_right * slope_left * (x_left - x_right) + slope_right * y_right - slope_left * y_left) / (slope_right - slope_left)
        return [x, y]

