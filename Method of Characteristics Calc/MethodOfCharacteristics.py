# Author: Vaishnav Srivaths (vsrivat@purdue.edu)
import math
import numpy as np
from equations import Equations
from rocketengine import rocketEngine
from charPoint import charPoint
equations = Equations

## INPUTS ##
numLines = 14
numPoints = int(numLines + numLines * (numLines + 1) / 2)
gamma = 1.4
R = 287 # J/kgK
P_e = 101.325 * 10 ** 3 # kPa
P_c = 1/0.1278 * P_e # kPa
T_c = 1000 # K
r_t = 0.0254 # m
initialFlowAngle = 0.19 # deg

engine = rocketEngine(P_c, T_c, P_e, gamma, R, r_t, numPoints, numLines)
M_e = engine.getMachExit(engine.P_c, engine.P_e, engine.gamma)
nu_e = equations.getPMAngle(M_e, gamma)
max_theta = nu_e / 2
points = []
print(max_theta)
thetaInterpolation = equations.interpolateTheta(initialFlowAngle, max_theta, numLines)
print(thetaInterpolation)

j = 1 + numLines
k = 0
for i in range(numPoints):
    points.append(charPoint(i + 1, False, False))
    currentPoint = points[i]
    if currentPoint.number == j + k:
        currentPoint.setWallPoint(True)
        k += 1
        j = j + numLines - k
    else:
        currentPoint.setWallPoint(False)
for i in range(numPoints):
    currentPoint = points[i]
    if i != 0:
        if points[i-1].isWallPoint:
            currentPoint.setCenterPoint(True)
            currentPoint.setTheta(0)
            currentPoint.setY(0)
        else:
            currentPoint.setCenterPoint(False)
    if i == 0:
        currentPoint.setCenterPoint(True)
        currentPoint.setTheta(0)
        currentPoint.setY(0)

for i in range(numPoints):
    currentPoint = points[i]
    print(currentPoint.number, currentPoint.isWallPoint, currentPoint.isCenterPoint)

# throat point (point a)
x_a = 0
theta_a = initialFlowAngle
nu_a = initialFlowAngle
M_a = equations.getMachFromPM(nu_a, gamma)
nu_1 = 2 * nu_a
M_1 = equations.getMachFromPM(nu_1, gamma)
mu_a = equations.getMachAngle(M_a)
mu_1 = equations.getMachAngle(M_1)
slope_a1 = ((theta_a - mu_a) + (0 - mu_1)) / 2
print(slope_a1)


# point 2 to numLines + 1
for index in range(numLines + 1):
    if(index == 0):
        # point 1 (at center line) properties
        points[index].setX(r_t * math.tan((90 + slope_a1) * math.pi / 180))
        print(points[index].x)
        points[index].setNu(theta_a + nu_a)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        points[index].setMu(equations.getMachAngle(points[index].mach))
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)
        print(M_e)
        print(M_a)
        print(M_1)
        print("point: ", index + 1)
        print("   theta: ", points[index].theta)
        print("   nu: ", points[index].nu)
        print("   mach: ", points[index].mach)
        print("   mu: ", points[index].mu)
        print("   coordinate: ", points[index].x, points[index].y)
        
    elif(points[index].isWallPoint == False):
        #points after point 1 in the same reflection
        print("point: ", index + 1)
        points[index].setTheta(thetaInterpolation[index])
        print("   theta: ", points[index].theta)
        points[index].setNu(thetaInterpolation[index])
        print("   nu: ", points[index].nu)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        print("   mach: ", points[index].mach)
        points[index].setMu(equations.getMachAngle(points[index].mach))
        print("   mu: ", points[index].mu)
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)

        theta_nu_ax = thetaInterpolation[index] + (points[index - 1].theta - points[index - 1].nu) / 2
        M_ax = equations.getMachFromPM(theta_nu_ax, gamma)
        mu_ax = equations.getMachAngle(M_ax)

        slope_left = (points[index - 1].theta + points[index - 1].mu + points[index].theta + points[index].mu) / 2
        slope_right = (theta_nu_ax - mu_ax + points[index].theta - points[index].mu) / 2

        intersectionPoint = equations.intersection(x_a, r_t, slope_right, points[index - 1].x, points[index - 1].y, slope_left)
        points[index].setX(intersectionPoint[0])
        points[index].setY(intersectionPoint[1])
        print("   left running slope: ", slope_left)
        print("   right running slope: ", slope_right)
        print("   coordinate: ", points[index].x, points[index].y)
    elif(points[index].isWallPoint):
        # for first wall point
        points[index].setTheta(points[index - 1].theta)
        points[index].setNu(points[index - 1].nu)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        points[index].setMu(equations.getMachAngle(points[index].mach))
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)
        print("point: ", index + 1)
        print("   theta: ", points[index].theta)
        print("   nu: ", points[index].nu)
        print("   mach: ", points[index].mach)
        print("   mu: ", points[index].mu)

        slope_left = (points[index].theta + points[index].mu + points[index - 1].theta + points[index - 1].mu) / 2
        slope_second_left = points[index - 1].theta

        intersectionPoint = equations.intersection(x_a, r_t, slope_second_left, points[index - 1].x, points[index - 1].y, slope_left)
        print(intersectionPoint)
        points[index].setX(intersectionPoint[0])
        points[index].setY(intersectionPoint[1])
        print("   left running slope: ", slope_left)
        print("   second left running slope: ", slope_second_left)
        print("   coordinate: ", points[index].x, points[index].y)


        

jindex = 0 # how much to subtract from numLines to reach left point
kindex = 1 # delta Theta index, ends up pointing to 0
for index in range(numLines + 1, numPoints):
    if(points[index].isCenterPoint):
        # need nu, mach, mu, left slope, intersection point. (already set theta as 0 and right slope is 0)
        leftPoint = points[index - numLines + jindex]
        previousCenterPoint = points[index - numLines + jindex - 1]
        print("point: ", index + 1)
        print("    reference point: ", index - numLines + jindex + 1)
        points[index].setNu(leftPoint.nu + leftPoint.theta)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        points[index].setMu(equations.getMachAngle(points[index].mach))
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)
        print("   theta: ", points[index].theta)
        print("   nu: ", points[index].nu)
        print("   mach: ", points[index].mach)
        print("   mu: ", points[index].mu)

        slope_left = 0
        slope_right = (points[index].theta - points[index].mu + leftPoint.theta - leftPoint.mu) / 2

        intersectionPoint = equations.intersection(leftPoint.x, leftPoint.y, slope_right, previousCenterPoint.x, previousCenterPoint.y, slope_left)
        points[index].setX(intersectionPoint[0])
        points[index].setY(intersectionPoint[1])
        print("   left running slope: ", slope_left)
        print("   right running slope: ", slope_right)
        print("   coordinate: ", points[index].x, points[index].y)
    elif(points[index].isCenterPoint == False and points[index].isWallPoint == False):
        # need theta, nu, mach, mu, left slope, intersection point
        leftPoint = points[index - numLines + jindex]
        previousPoint = points[index - 1]
        print("point: ", index + 1)
        print("    reference point 1: ", index - numLines + jindex + 1)
        print("    reference point 2: ", index - 1 + 1)
        points[index].setTheta(thetaInterpolation[kindex])
        points[index].setNu((leftPoint.nu + leftPoint.theta + previousPoint.nu - previousPoint.theta) / 2)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        points[index].setMu(equations.getMachAngle(points[index].mach))
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)
        print("   theta: ", points[index].theta)
        print("   nu: ", points[index].nu)
        print("   mach: ", points[index].mach)
        print("   mu: ", points[index].mu)

        slope_left = (previousPoint.theta + previousPoint.mu + points[index].theta + points[index].mu) / 2
        slope_right = (points[index].theta - points[index].mu + leftPoint.theta - leftPoint.mu) / 2

        intersectionPoint = equations.intersection(leftPoint.x, leftPoint.y, slope_right, previousPoint.x, previousPoint.y, slope_left)
        
        points[index].setX(intersectionPoint[0])
        points[index].setY(intersectionPoint[1])
        print("   left running slope: ", slope_left)
        print("   right running slope: ", slope_right)
        print("   coordinate: ", points[index].x, points[index].y)
        
        kindex += 1
    elif(points[index].isWallPoint):
        # need theta, nu, mach, mu, right slope, intersection point (left slope is theta from reference point)
        leftPoint = points[index - numLines + jindex]
        previousPoint = points[index - 1]
        print("point: ", index + 1)
        print("    reference point : ", index - 1 + 1)
        points[index].setTheta(previousPoint.theta)
        points[index].setNu(previousPoint.nu)
        points[index].setMach(equations.getMachFromPM(points[index].nu, gamma))
        points[index].setMu(equations.getMachAngle(points[index].mach))
        points[index].getPosConstant(points[index].theta, points[index].nu)
        points[index].getNegConstant(points[index].theta, points[index].nu)
        print("   theta: ", points[index].theta)
        print("   nu: ", points[index].nu)
        print("   mach: ", points[index].mach)
        print("   mu: ", points[index].mu)
        slope_left = (points[index].theta + points[index].mu + previousPoint.theta + previousPoint.mu) / 2
        slope_second_left = previousPoint.theta

        intersectionPoint = equations.intersection(leftPoint.x, leftPoint.y, slope_second_left, previousPoint.x, previousPoint.y, slope_left)
        points[index].setX(intersectionPoint[0])
        points[index].setY(intersectionPoint[1])
        print("   left running slope: ", slope_left)
        print("   second left running slope: ", slope_second_left)
        print("   coordinate: ", points[index].x, points[index].y)

        kindex = 1 # resets delta theta
        jindex += 1 # moves over one reflection

pointNumbers = [0]
x_points = [0]
y_points = [r_t]
machs = [1]
flow_angles = [0]
mach_angles = [90]
pm_angles = [0]
for index in range(numPoints):
    pointNumbers.append(points[index].number)
    x_points.append(points[index].x)
    y_points.append(points[index].y)
    machs.append(points[index].mach)
    flow_angles.append(points[index].theta)
    mach_angles.append(points[index].mu)
    pm_angles.append(points[index].nu)

data = np.asarray([pointNumbers, 
        x_points, 
        y_points, 
        machs, 
        flow_angles, 
        mach_angles, 
        pm_angles]).transpose()
print(data)

np.savetxt('EnginePoints.csv', data, delimiter=',')