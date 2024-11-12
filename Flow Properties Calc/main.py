# isentropic relations
# t/t0 = (1 + (k-1)/2 * M^2)^-1
import math
import sympy
from sympy.abc import x
import numpy
from isentropic import Isentropic
from normalShock import NormalShock
from fanno import Fanno
#from rayleigh import Rayleigh

isentropic = Isentropic
normalShock = NormalShock
fanno = Fanno
#rayleigh = Rayleigh

gamma = float(input("Ratio of specific heats: ")) # ratio of specific heats

print('\nHere are the following relations:\n1. Isentropic Flow Relations\n2. Normal Shock Relations\n3. Oblique Shock Relations\n4. Fanno Flow Relations\n5. Rayleigh Flow Relations')
relation = int(input("\nWhich relation would you like to calculate properties for? Enter a number: "))

if relation == 1:
    print('Here are possible inputs:\n1. Mach number\n2. T/T0\n3. p/p0\n4. rho/rho0\n5. A/A*(subsonic)\n6. A/A*(supersonic)\n7. Mach angle(degrees)\n8. Prandtl-Meyer Angle(degrees)')
    choice = int(input("\nWhich property do you already know? Enter a number: "))

    if choice == 1:
        M = float(input("Mach number = "))
        p_p0 = isentropic.getP_P0(M, gamma)
        print("p/p0 = ", p_p0)
        t_t0 = isentropic.getT_T0(M, gamma)
        print("T/T0 = ", t_t0)
        rho_rho0 = isentropic.getRho_Rho0(M, gamma)
        print("rho/rho0 = ", rho_rho0)
        alpha = isentropic.getAlpha(M)
        print("Mach angle (degrees) = ", alpha)
        nu = isentropic.getNu(M, gamma)
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)
    if choice == 2:
        t_t0 = float(input("T/T0 = "))
        M = isentropic.getMfromTR(t_t0, gamma)
        print("Mach number = ", M)
        p_p0 = isentropic.getP_P0(M, gamma)
        print("p/p0 = ", p_p0)
        rho_rho0 = isentropic.getRho_Rho0(M, gamma)
        print("rho/rho0 = ", rho_rho0)
        alpha = isentropic.getAlpha(M)
        print("Mach angle (degrees) = ", alpha)
        nu = isentropic.getNu(M, gamma)
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)
    if choice == 3:
        p_p0 = float(input("p/p0 = "))
        M = isentropic.getMfromPR(p_p0, gamma)
        print("Mach number = ", M)
        t_t0 = isentropic.getT_T0(M, gamma)
        print("T/T0 = ", t_t0)
        rho_rho0 = isentropic.getRho_Rho0(M, gamma)
        print("rho/rho0 = ", rho_rho0)
        alpha = isentropic.getAlpha(M)
        print("Mach angle (degrees) = ", alpha)
        nu = isentropic.getNu(M, gamma)
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)
    if choice == 4:
        rho_rho0 = float(input("rho/rho0 = "))
        M = isentropic.getMfromDR(rho_rho0, gamma)
        print("Mach number = ", M)
        p_p0 = isentropic.getMfromPR(M, gamma)
        print("p/p0 = ", p_p0)
        t_t0 = isentropic.getT_T0(M, gamma)
        print("T/T0 = ", t_t0)
        alpha = isentropic.getAlpha(M)
        print("Mach angle (degrees) = ", alpha)
        nu = isentropic.getNu(M, gamma)
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)
    if choice == 7: 
        alpha = float(input("Mach angle (degrees) = "))
        M =  isentropic.getMfromAlpha(alpha, gamma)
        print("Mach number = ", M)
        p_p0 = isentropic.getMfromPR(M, gamma)
        print("p/p0 = ", p_p0)
        t_t0 = isentropic.getT_T0(M, gamma)
        print("T/T0 = ", t_t0)
        rho_rho0 = isentropic.getRho_Rho0(M, gamma)
        print("rho/rho0 = ", rho_rho0)
        nu = isentropic.getNu(M, gamma)
        print("Prandtl-Meyer angle (degrees) = ", nu)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)
    if choice == 8:
        nu = float(input("Prandtl-Meyer angle (degrees) = "))
        M = isentropic.getMfromNu(nu, gamma)
        print("Mach number = ", M)
        p_p0 = isentropic.getMfromPR(M, gamma)
        print("p/p0 = ", p_p0)
        t_t0 = isentropic.getT_T0(M, gamma)
        print("T/T0 = ", t_t0)
        rho_rho0 = isentropic.getRho_Rho0(M, gamma)
        print("rho/rho0 = ", rho_rho0)
        alpha = isentropic.getAlpha(M)
        print("Mach angle (degrees) = ", alpha)
        A_Astar = isentropic.getA_Astar(M, gamma)
        print("A/A* = ", A_Astar)

if relation == 2:
    print('Here are possible inputs:\n1. M1\n2. M2\n3. p2/p1\n4. rho2/rho1\n5. T2/T1\n6. p02/p01')
    choice = int(input("\nWhich property do you already know? Enter a number: "))

    if choice == 1:
        M1 = float(input("M1 = "))
        M2 = normalShock.getM2(M1, gamma)
        print("M2 = ", M2)
        p2_p1 = normalShock.getP2_P1(M1, gamma)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = normalShock.getRho2_Rho1(M1, gamma)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = normalShock.getT2_T1(M1, gamma)
        print("T2/T1 = ", t2_t1)
        p02_p01 = normalShock.getP02_P01(M1, gamma)
        print("P02/P01 = ", p02_p01)
    if choice == 2:
        M2 = float(input("M2 = "))
        M1 = normalShock.getM1fromM2(M2, gamma)
        print("M1 = ", M1)
        p2_p1 = normalShock.getP2_P1(M1, gamma)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = normalShock.getRho2_Rho1(M1, gamma)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = normalShock.getT2_T1(M1, gamma)
        print("T2/T1 = ", t2_t1)
        p02_p01 = normalShock.getP02_P01(M1, gamma)
        print("P02/P01 = ", p02_p01)
    if choice == 3:
        p2_p1 = float(input("P2/P1 = "))
        M1 = normalShock.getM1fromPR(p2_p1, gamma)
        print("M1 = ", M1)
        M2 = normalShock.getM2(M1, gamma)
        print("M2 = ", M2)
        rho2_rho1 = normalShock.getRho2_Rho1(M1, gamma)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = normalShock.getT2_T1(M1, gamma)
        print("T2/T1 = ", t2_t1)
        p02_p01 = normalShock.getP02_P01(M1, gamma)
        print("P02/P01 = ", p02_p01)
    if choice == 4:
        rho2_rho1 = float(input("rho2/rho1 = "))
        M1 = normalShock.getM1fromDR(rho2_rho1, gamma)
        print("M1 = ", M1)
        M2 = normalShock.getM2(M1, gamma)
        print("M2 = ", M2)
        p2_p1 = normalShock.getP2_P1(M1, gamma)
        print("P2/P1 = ", p2_p1)
        t2_t1 = normalShock.getT2_T1(M1, gamma)
        print("T2/T1 = ", t2_t1)
        p02_p01 = normalShock.getP02_P01(M1, gamma)
        print("P02/P01 = ", p02_p01)
    if choice == 5:
        t2_t1 = float(input("T2/T1 = "))
        M1 = normalShock.getM1fromTR(t2_t1, gamma)
        print("M1 = ", M1)
        M2 = normalShock.getM2(M1, gamma)
        print("M2 = ", M2)
        p2_p1 = normalShock.getP2_P1(M1, gamma)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = normalShock.getRho2_Rho1(M1, gamma)
        print("rho2/rho1 = ", rho2_rho1)
        p02_p01 = normalShock.getP02_P01(M1, gamma)
        print("P02/P01 = ", p02_p01)
    if choice == 6:
        p02_p01 = float(input("P02/P01 = "))
        M1 = normalShock.getM1fromP0R(p02_p01, gamma)
        print("M1 = ", M1)
        M2 = normalShock.getM2(M1, gamma)
        print("M2 = ", M2)
        p2_p1 = normalShock.getP2_P1(M1, gamma)
        print("P2/P1 = ", p2_p1)
        rho2_rho1 = normalShock.getRho2_Rho1(M1, gamma)
        print("rho2/rho1 = ", rho2_rho1)
        t2_t1 = normalShock.getT2_T1(M1, gamma)
        print("T2/T1 = ", t2_t1)

if relation == 3:
    M1 = float(input("M1 = "))
    a = float(input("turn angle (deg) = "))

if relation == 4:
    print('Here are possible inputs:\n1. M\n2. t/t*\n3. p/p*\n4. p0/p0* (sub)\n5. p0/p0* (sup)\n6. u/u*\n7. 4fL*/D (sub)\n8. 4fL*/D (sup)')
    choice = int(input("\nWhich property do you already know? Enter a number: "))

    if choice == 1:
        M = float(input("M = "))
        t_tstar = fanno.getTR(M, gamma)
        print("t0/t0* = ", t_tstar)
        p_pstar = fanno.getPR(M, gamma)
        print("p/p* = ", p_pstar)
        p0_p0star = fanno.getP0R(M, gamma)
        print("p0/p0* = ", p0_p0star)
        u_ustar = fanno.getUR(M, gamma)
        print("u/u* = ", u_ustar)
        fric = fanno.getFric(M, gamma)
        print("4fL*/D = ", fric)

    if choice == 2:
        t_tstar = float(input("t/t* = "))
        if t_tstar >= (gamma + 1) / 2:
            print("error: t/t* must be less than ", (gamma + 1) / 2)
        else:
            M = fanno.getMfromTR(t_tstar, gamma)
            print("M = ", M)
            p_pstar = fanno.getPR(M, gamma)
            print("p/p* = ", p_pstar)
            p0_p0star = fanno.getP0R(M, gamma)
            print("p0/p0* = ", p0_p0star)
            u_ustar = fanno.getUR(M, gamma)
            print("u/u* = ", u_ustar)
            fric = fanno.getFric(M, gamma)
            print("4fL*/D = ", fric)

    if choice == 3:
        p_pstar = float(input("p/p* = "))
        M = fanno.getMfromPR(p_pstar, gamma)
        print("M = ", M)
        t_tstar = fanno.getTR(M, gamma)
        print("t0/t0* = ", t_tstar)
        p0_p0star = fanno.getP0R(M, gamma)
        print("p0/p0* = ", p0_p0star)
        u_ustar = fanno.getUR(M, gamma)
        print("u/u* = ", u_ustar)
        fric = fanno.getFric(M, gamma)
        print("4fL*/D = ", fric)

    if choice == 4:
        p0_p0star = float(input("p0/p0* (sub) = "))
        M = fanno.getMfromP0Rsub(p0_p0star, gamma)
        print("M = ", M)
        p_pstar = fanno.getPR(M, gamma)
        print("p/p* = ", p_pstar)
        t_tstar = fanno.getTR(M, gamma)
        print("t0/t0* = ", t_tstar)
        u_ustar = fanno.getUR(M, gamma)
        print("u/u* = ", u_ustar)
        fric = fanno.getFric(M, gamma)
        print("4fL*/D = ", fric)

    if choice == 5:
        p0_p0star = float(input("p0/p0* (sup) = "))
        M = fanno.getMfromP0Rsup(p0_p0star, gamma)
        print("M = ", M)
        p_pstar = fanno.getPR(M, gamma)
        print("p/p* = ", p_pstar)
        t_tstar = fanno.getTR(M, gamma)
        print("t0/t0* = ", t_tstar)
        u_ustar = fanno.getUR(M, gamma)
        print("u/u* = ", u_ustar)
        fric = fanno.getFric(M, gamma)
        print("4fL*/D = ", fric)