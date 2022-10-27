# simple example of classes -- here we have a projective class
# that returns the height and distance for a projectile with
# specified initial velocity and angle

import math

class projectile:

    # this is called every time a new object is created
    def __init__ (self, v=1.0, theta=45, grav=9.81):

        self.v = v           # velocity m/s
        self.theta = theta   # angle (degrees)
        
        self.thetaRad = theta*math.pi/180.0
        self.vx = v*math.cos(self.thetaRad)
        self.vy = v*math.sin(self.thetaRad)

        self.g = grav


    def height(self):

        # how high does this projectile go?
        # vf_y^2 = 0 = vi_y^2 - 2 g h
        h = self.vy**2/(2.0*self.g)

        return h

    def distance(self):
        
        # time of flight up
        # vf_y = 0 = vi_y - g t
        t = self.vy/self.g

        # total time = up + down
        t = 2.0*t

        d = self.vx*t

        return d

    def __str__(self):
        # a string representation for this class -- so we can print it
        str = " v: %s m/s\n theta: %s (degrees)\n height: %s m\n distance: %s m\n" % \
            (`self.v`, `self.theta`, `self.height()`, `self.distance()`)
        
        return str

# fire a projectile
p1 = projectile()
print p1

# create a list of projectiles
projectiles = []
projectiles.append(p1)

projectiles.append(projectile(v = 100, theta = 70))
projectiles.append(projectile(v = 1000, theta = 30))

print projectiles

print " "
print "heights:"
for p in projectiles:
    print p.height()



        
