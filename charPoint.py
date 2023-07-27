# Author: Vaishnav Srivaths (vsrivat@purdue.edu)
class charPoint:
    def __init__(self, number, wallPoint, centerPoint):
        self.number = number
        self.wallPoint = wallPoint
        self.centerPoint = centerPoint
    
    def setWallPoint(self, isWallPoint):
        self.isWallPoint = isWallPoint

    def setCenterPoint(self, isCenterPoint):
        self.isCenterPoint = isCenterPoint
        
    def setTheta(self, theta):
        self.theta = theta
    
    def setNu(self, nu):
        self.nu = nu
    
    def setMu(self, mu):
        self.mu = mu
    
    def setMach(self, mach):
        self.mach = mach
    
    def setX(self, x):
        self.x = x

    def setY(self, y):
        self.y = y
    
    def getX(self):
        return self.x_coord
    
    def getY(self):
        return self.y_coord
    
    def getPosConstant(self, theta, nu):
        self.posConstant = theta - nu
        return self.posConstant
    
    def getNegConstant(self, theta, nu):
        self.negConstant = theta + nu
        return self.negConstant
    
    