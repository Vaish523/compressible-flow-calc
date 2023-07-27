# isentropic relations
# t/t0 = (1 + (k-1)/2 * M^2)^-1
import math
import sympy
from sympy.abc import x
import numpy

gamma = float(input("Ratio of specific heats: ")) # ratio of specific heats

print('\nHere are the following relations:\n1. Isentropic Flow Relations\n2. Normal Shock Relations\n3. Oblique Shock Relations\n4. Conical Shock Relations\n5. Fanno Flow Relations\n6. Rayleigh Flow Relations')
relation = int(input("\nWhich relation would you like to calculate properties for? Enter a number: "))

if relation == 1:
    print('Here are possible inputs:\n1. Mach number\n2. T/T0\n3. p/p0\n4. rho/rho0\n5. A/A*(subsonic)\n6. A/A*(supersonic)\n7. Mach angle(degrees)\n8. Prandtl-Meyer Angle(degrees)')
    choice = int(input("\nWhich property do you already know? Enter a number: "))

    if choice == 1:
        M = float(input("Mach number = "))
        p_p0 = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (1 - gamma))
        print("p/p0 = ", p_p0)
        t_t0 = (1 + (gamma - 1) / 2 * M ** 2) ** (-1)
        print("T/T0 = ", t_t0)
        rho_rho0 = (1 + (gamma - 1) / 2 * M ** 2) ** (1 / (1 - gamma))
        print("rho/rho0 = ", rho_rho0)
        alpha = math.asin(1 / M) * 180 / math.pi
        print("Mach angle (degrees) = ", alpha)
        nu = 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) + math.atan(math.sqrt(M ** 2 - 1)))
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)
    if choice == 2:
        t_t0 = float(input("T/T0 = "))
        M = math.sqrt((2 * (1/t_t0 - 1)) / (gamma - 1))
        print("Mach number = ", M)
        p_p0 = t_t0 ** (gamma / (gamma - 1))
        print("p/p0 = ", p_p0)
        rho_rho0 = t_t0 ** (1 / (gamma - 1))
        print("rho/rho0 = ", rho_rho0)
        alpha = math.asin(1 / M) * 180 / math.pi
        print("Mach angle (degrees) = ", alpha)
        nu = 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) + math.atan(math.sqrt(M ** 2 - 1)))
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)
    if choice == 3:
        p_p0 = float(input("p/p0 = "))
        M = math.sqrt((2 * (p_p0 ** (1 / gamma - 1) - 1)) / (gamma - 1))
        print("Mach number = ", M)
        t_t0 = p_p0 ** ((gamma - 1)/ gamma)
        print("T/T0 = ", t_t0)
        rho_rho0 = p_p0 ** (1 / gamma)
        print("rho/rho0 = ", rho_rho0)
        alpha = math.asin(1 / M) * 180 / math.pi
        print("Mach angle (degrees) = ", alpha)
        nu = 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) + math.atan(math.sqrt(M ** 2 - 1)))
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)
    if choice == 4:
        rho_rho0 = float(input("rho/rho0 = "))
        M = math.sqrt((2 * (rho_rho0 ** (1 - gamma) - 1)) / (gamma - 1))
        print("Mach number = ", M)
        p_p0 = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (1 - gamma))
        print("p/p0 = ", p_p0)
        t_t0 = p_p0 ** ((gamma - 1)/ gamma)
        print("T/T0 = ", t_t0)
        alpha = math.asin(1 / M) * 180 / math.pi
        print("Mach angle (degrees) = ", alpha)
        nu = 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) + math.atan(math.sqrt(M ** 2 - 1)))
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)
    if choice == 7: 
        alpha = float(input("Mach angle (degrees) = "))
        M =  1 / math.sin(alpha * math.pi / 180)
        print("Mach number = ", M)
        p_p0 = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (1 - gamma))
        print("p/p0 = ", p_p0)
        t_t0 = (1 + (gamma - 1) / 2 * M ** 2) ** (-1)
        print("T/T0 = ", t_t0)
        rho_rho0 = (1 + (gamma - 1) / 2 * M ** 2) ** (1 / (1 - gamma))
        print("rho/rho0 = ", rho_rho0)
        nu = 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (M ** 2 - 1))) + math.atan(math.sqrt(M ** 2 - 1)))
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)
    if choice == 8:
        nu = float(input("Prandtl-Meyer angle (degrees) = "))
        M_sols = sympy.solve(nu - 180 / math.pi * (-1 * math.sqrt((gamma + 1) / (gamma - 1)) * math.atan(math.sqrt((gamma - 1) / (gamma + 1) * (x ** 2 - 1))) + math.atan(math.sqrt(x ** 2 - 1))), "x")
        M = M_sols[1]
        print("Mach number = ", M)
        p_p0 = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (1 - gamma))
        print("p/p0 = ", p_p0)
        t_t0 = (1 + (gamma - 1) / 2 * M ** 2) ** (-1)
        print("T/T0 = ", t_t0)
        rho_rho0 = (1 + (gamma - 1) / 2 * M ** 2) ** (1 / (1 - gamma))
        print("rho/rho0 = ", rho_rho0)
        alpha = math.asin(1 / M) * 180 / math.pi
        print("Mach angle (degrees) = ", alpha)
        A_Astar = 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / (1 + (gamma - 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))
        print("A/A* = ", A_Astar)

if relation == 2:
    print('Here are possible inputs:\n1. M1\n2. M2\n3. p2/p1\n4. rho2/rho1\n5. T2/T1\n6. p02/p01')
    choice = int(input("\nWhich property do you already know? Enter a number: "))

    if choice == 1:
        M1 = float(input("M1 = "))
        M2 = math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
        print("M2 = ", M2)
        p2_p1 = (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
        print("T2/T1 = ", t2_t1)
        p02_p01 = ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
        print("P02/P01 = ", p02_p01)
    if choice == 2:
        M2 = float(input("M2 = "))
        M1 = math.sqrt( ((gamma - 1) * M2 ** 2 + 2) / (2 * gamma * M2 ** 2 - (gamma - 1)) )
        print("M1 = ", M1)
        p2_p1 = (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
        print("T2/T1 = ", t2_t1)
        p02_p01 = ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
        print("P02/P01 = ", p02_p01)
    if choice == 3:
        p2_p1 = float(input("P2/P1 = "))
        M1 = math.sqrt((p2_p1 * (gamma + 1) + (gamma - 1)) / (2 * gamma))
        print("M1 = ", M1)
        M2 = math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
        print("M2 = ", M2)
        rho2_rho1 = ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
        print("T2/T1 = ", t2_t1)
        p02_p01 = ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
        print("P02/P01 = ", p02_p01)
    if choice == 4:
        rho2_rho1 = float(input("rho2/rho1 = "))
        M1 = math.sqrt((-2 * p2_p1) / (p2_p1 * (gamma - 1) - (gamma + 1)))
        print("M1 = ", M1)
        M2 = math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
        print("M2 = ", M2)
        p2_p1 = (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
        print("P2/P1 = ", p2_p1)
        t2_t1 = ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
        print("T2/T1 = ", t2_t1)
        p02_p01 = ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
        print("P02/P01 = ", p02_p01)
    if choice == 5:
        t2_t1 = float(input("T2/T1 = "))
        M1_sols = sympy.solve(t2_t1 - ((2 * gamma * x ** 2 - (gamma - 1)) * ((gamma - 1) * x ** 2 + 2)) / ((gamma + 1) ** 2 * x ** 2), "x")
        M1 = M1_sols[1]
        print("M1 = ", M1)
        M2 = math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
        print("M2 = ", M2)
        p2_p1 = (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
        print("rho2/rho1 = ", rho2_rho1)
        p02_p01 = ((((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * M1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
        print("P02/P01 = ", p02_p01)
    if choice == 6:
        # FIX
        p02_p01 = float(input("P02/P01 = "))
        M1_sols = sympy.solve(p02_p01 - ((((gamma + 1) * x ** 2) / ((gamma - 1) * x ** 2 + 2)) ** (gamma / (gamma - 1))) * (((gamma + 1) / (2 * gamma * x ** 2 - (gamma - 1))) ** (1 / (gamma - 1))), "x")
        M1 = M1_sols
        print("M1 = ", M1)
        M2 = math.sqrt( ((gamma - 1) * M1 ** 2 + 2) / (2 * gamma * M1 ** 2 - (gamma - 1)) )
        print("M2 = ", M2)
        p2_p1 = (2 * gamma * M1 ** 2 - (gamma - 1)) / (gamma + 1)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = ((gamma + 1) * M1 ** 2) / ((gamma - 1) * M1 ** 2 + 2)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = ((2 * gamma * M1 ** 2 - (gamma - 1)) * ((gamma - 1) * M1 ** 2 + 2)) / ((gamma + 1) ** 2 * M1 ** 2)
        print("T2/T1 = ", t2_t1)
