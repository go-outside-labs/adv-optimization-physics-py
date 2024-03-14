########################################################################
# This program illustrate many approaches to solve ODEs. In this case
#    we use the simple pendulum as the system. 
#
# This module implements all the classes and functions for the problem.
#
#    MIA STEIN, spring/2013 
#    (based on Mike Zingale's codes)
#
########################################################################


import math 
import numpy as np



class pendulumtHistory:
    """ container to store the integrated history of the pendulum.  """
    
    def __init__(self):
        self.t = None
        self.theta = None
        self.omega = None
        self.E = None          
        self.dE = None
        
    def energyError(self, EN,E0):
        """ absolute energy error given by the starting and ending point """
        return abs( (EN - E0)/ E0)



class pendulumTrajectory:
    """ hold the initial conditions the pendulum and integrate """
    
    def __init__(self, CONST_theta0, CONST_l, CONST_g, CONST_m):   
        self.theta0_degrees = CONST_theta0                 # initial condition -> theta0 in degree
        self.theta0 = self.theta0_degrees*math.pi/180      # convert theta0 to radians
        self.t0 = 0.0                                       # time starts at zero
        self.omega0 = 0.0                                   # initial condition -> omega0 starts at 0
        self.l = CONST_l                                    # length of pendulum
        self.g = CONST_g                                    # acceleration of gravity
        self.m = CONST_m                                    # mass
    
    
    def pendulumPeriod(self):
        """ return the estimative of the period of the pendulum, which 
                only depend on the theta0 (initial displacement). """
        return 2.0*math.pi*math.sqrt(self.l/self.g)*(1.0 + self.theta0**2*(1.0/16.0))


    def pendulumEnergy(self):
        """ return the estimative of the energy of the pendulum, 
                in function of the initial values. The total energy
                should be conserved so this value should be a 
                constant."""
        return abs(0.5*self.m*(self.l)**2*(self.omega0)**2 - self.m*self.g*self.l*(math.cos(self.theta0)))
    
    
    def newEnergy(self, the, ome):
        """ return the energy of the pendulum by substituting the time steps,
                we use this to show that Euler is not good for the pendulum
                system since the total energy seems to monotonically increases
                in time. This wont happen on the other methods."""    
        return abs(0.5*self.m*(self.l)**2*(ome)**2 - self.m*self.g*self.l*(math.cos(the)))


    def intEuler(self, dt, tmax):
        """ integrate the equations of motion using Euler's method. """

        # initial conditions
        t = self.t0
        theta = self.theta0
        omega = self.omega0
        e = self.pendulumEnergy()

        # store the history for plotting
        tpoints = [t]
        thetapoints = [theta]
        omegapoints = [omega]
        energypoints = [e]
        
        # loop over the time steps
        while (t < tmax):
            
            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t+dt > tmax:
                dt = tmax-t            

            # get the RHS
            thetadot, omegadot = self.rhs(theta, omega)

            # advance
            thetanew = theta + dt*thetadot
            omeganew = omega + dt*omegadot
            t += dt
            
            # calculates the step energy 
            energynew = self.newEnergy(theta,omega)
        
            # set for the next step
            theta = thetanew
            omega = omeganew

            # store
            tpoints.append(t)
            thetapoints.append(thetanew*180/(math.pi))
            omegapoints.append(omeganew)
            energypoints.append(energynew)

        
        # return a orbitHistory object with the angular displacement
        H = pendulumtHistory()
        H.t = np.array(tpoints)
        H.theta = np.array(thetapoints) 
        H.omega = np.array(omegapoints)
        H.energy = np.array(energypoints)
        
        return H
    
    
    def intEulerCromer(self, dt, tmax):
        """ integrate the equations of motion using Euler-Cromer's method."""

        # initial conditions
        t = self.t0
        theta = self.theta0
        omega = self.omega0
        e = self.pendulumEnergy()

        # store the history for plotting
        tpoints = [t]
        thetapoints = [theta]
        omegapoints = [omega]
        energypoints = [e]

        # loop over the time steps
        while (t < tmax):
            
            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t+dt > tmax:
                dt = tmax-t            

            # get the RHS
            thetadot, omegadot = self.rhs(theta, omega)

            # advance
            omeganew = omega + dt*omegadot
            
            # these line is the only change from Euler
            thetanew = theta + dt*omeganew
    
            t += dt
            
            # calculates the step energy 
            energynew = self.newEnergy(theta,omega)
            
            # set for the next step
            theta = thetanew
            omega = omeganew

            # store
            tpoints.append(t)
            thetapoints.append(thetanew*180/(math.pi))
            omegapoints.append(omeganew)
            energypoints.append(energynew)


        
        # return a orbitHistory object with the angular displacement
        H = pendulumtHistory()
        H.t = np.array(tpoints)
        H.theta = np.array(thetapoints) 
        H.omega = np.array(omegapoints)
        H.energy = np.array(energypoints)
        
        return H

    
    
    def intVerlet(self, Dt, tmax):
        """ integrate the equations of motion using Euler's method. """

        # initial conditions
        t = self.t0
        theta = self.theta0
        omega = self.omega0
        energy_ver = self.pendulumEnergy()
        dt = Dt
        
        while (t < tmax):
            
            # make sure that the next step doesn't take us past where
            # we want to be, because of roundoff
            if t+dt > tmax:
                dt = tmax-t            
            # Take one backward step to start Verlet

            #alpha = - (self.g/self.l)*math.sin(theta)
            #thetaold = theta - omega*dt + 0.5*dt**2*alpha
            # get the RHS
            thetadot, alpha = self.rhs(theta, omega)
    

            # advance
            thetanew = theta + dt*omega +0.5*dt**2*alpha
            thetadot, alphaplus1 = self.rhs(thetanew, omega)
            omeganew = omega + 0.5*dt*(alpha+alphaplus1) 
            
            # calculate energy  for this time step
            energynew = self.newEnergy(thetanew,omeganew)            
            t += dt
            # set for the next step
            theta = thetanew;
            omega = omeganew
            
        # return a orbitHistory object with the angular displacement
        H = pendulumtHistory()
        H.t = Dt
        H.E = energynew
        H.dE = H.energyError(energynew,energy_ver)
        
        return H
    
    
    def rhs(self, the, ome):
        """ RHS of the equations of motion """
        thetadot = ome      
        omegadot = -(self.g/self.l)*math.sin(the)

        return thetadot, omegadot

