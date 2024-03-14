"""
    This module defines the motion of a damped driven pendulum:
        m l d2teta/dt2 = Fg + fd + Fr
    and integrates it using the 4th Runge-Kutta method with a uniform
    timestep.
    
    MIA STEIN, spring/2013 (based on Mike Zingale's codes)

"""


import numpy as np
from math import sin, cos, pi, floor



class pendulumHistory:
    """ a simple container to store the pendulum history"""
    def __init__(self, t = None, theta = None, omega = None):
        self.t = np.array(t)
        self.theta = np.array(theta)
        self.omega = np.array(omega)



def rhs(b, theta, omega, t, q, omegad):
        """ RHS of the equations of motion """
        thetadot = omega      
        omegadot = - sin(theta) -q*omega + b*cos(omegad*t)

        return thetadot, omegadot
        


def intRK4(b, dt, tmax, rhs, theta0, q, omegad ):
        """ integrate the equations of motion using 4th order R-K """

        # initial conditions
        t = 0.0
        theta = theta0
        omega = 0.0

        # store the history for plotting
        tPoints = [t]
        thetaPoints = [theta]
        omegaPoints = [omega]

        
        while (t < tmax):
                   
            thetanew, omeganew =  RK4_singlestep(b, theta, omega, t, dt, rhs, theta0, q, omegad)     
            t += dt

            # set for the next step
            #if (thetanew < pi):
            #    theta = thetanew + 2*pi
            #if (thetanew > pi):
            #    theta = thetanew -2*pi
            #else:
            theta = thetanew
            omega = omeganew

            # store
            tPoints.append(t)
            thetaPoints.append(theta)
            omegaPoints.append(omega)


        # restrict theta to be [-pi,pi]
        for n in range(len(thetaPoints)):
            thetaPoints[n] += pi
            thetaPoints[n] -= 2.0*pi*floor(thetaPoints[n]/(2.0*pi))
            thetaPoints[n] -= pi

        # return a orbitHistory object with the trajectory
        H = pendulumHistory(tPoints, thetaPoints, omegaPoints)
        
        return H




def RK4_singlestep(b, theta0tpm, omega0tpm, t, dt, rhs, theta0, q, omegad):
    """ take a single RK-4 timestep from t to t+dt for the system 
        ydot = rhs """

    theta = theta0tpm
    omega = omega0tpm
    
    # get the RHS at several points
    thetadot1, omegadot1 = rhs(b, theta, omega,t, q, omegad)
    thetadot2, omegadot2 = rhs(b, theta+0.5*dt*thetadot1, omega+0.5*dt*omegadot1,t, q, omegad)
    thetadot3, omegadot3 = rhs(b, theta+0.5*dt*thetadot2, omega+0.5*dt*omegadot2,t, q, omegad)
    thetadot4, omegadot4 = rhs(b, theta+0.5*dt*thetadot3, omega+0.5*dt*omegadot3,t, q, omegad)
            
    # advance
    thetanew = theta + (dt/6.0)*(thetadot1 + 2.0*thetadot2 + 2.0*thetadot3 + thetadot4)
    omeganew = omega + (dt/6.0)*(omegadot1 + 2.0*omegadot2 + 2.0*omegadot3 + omegadot4)

    return thetanew, omeganew
