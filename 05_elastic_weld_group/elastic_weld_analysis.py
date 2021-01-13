import numpy as np

class Weld_Load(object):
    
    def __init__(self, name, pz=0.0, vx=0.0, vy=0.0, mx=0.0, my=0.0, mz=0.0, x=0.0, y=0.0, z=0.0):
        
        self.name = name

        # Applied load. Provide a diagram showing the sign convention. Units are k and in*k.
        self.load = np.array([vx, vy, pz, mx, my, mz])
        self.coords = np.array([x, y, z])

    def calc_cg(self, cg_coords):

        # Distance from load point to weld group centroid.
        self.coords_prime = self.coords - cg_coords

        # Applied moments and moments form forces away from the weld group centroid.
        self.load_prime = np.array([self.load[0],
                                    self.load[1],
                                    self.load[2],
                                    self.load[3] - self.load[1]*self.coords_prime[2] + self.load[2]*self.coords_prime[1], 
                                    self.load[4] + self.load[0]*self.coords_prime[2] - self.load[2]*self.coords_prime[0],
                                    self.load[5] - self.load[0]*self.coords_prime[1] + self.load[1]*self.coords_prime[0]])

class Weld(object):
    
    def __init__(self, name, xi=0.0, yi=0.0, xj=0.0, yj=0.0):
        
        self.name = name
        
        # Weld end point coordinate relative to origin. Welds are defined as [[xi, yi, xj, yj]]. Units are inches.
        self.i = np.array([xi, yi, 0.0])
        self.j = np.array([xj, yj, 0.0])
        self.coords = np.array([self.i, self.j])

        # Weld length in each orthogonal axis, x and y, and along the longitundianl axis of the weld [in].
        self.vec = self.coords[1] - self.coords[0]
        self.L = np.linalg.norm(self.vec)

        # Location of the x and y centroids of weld line [in].
        self.cg = np.mean(self.coords, axis=0)

        # Second moment of area per weld width about x and y axis [in^4/in]. The length of the weld is treated as an area.
        self.I = self.L*self.vec**2/12

        # Area moments [in^3]
        self.Lc = self.cg*self.L
        
    def calc_prime(self, cg_coords):

        # Coordinates relative to group centroid.
        self.i_prime = self.i - cg_coords
        self.j_prime = self.j - cg_coords
        self.d_prime = self.cg - cg_coords
        
        # Moment of inertia relative to group centroid using parallel axis theorem.
        # I_prime = [Iyy', Ixx', Izz'] = [Iyy, Ixx, Izz] + L*[dx, dy, dz]
        self.I_prime = self.I + self.L*self.d_prime**2
        
        # Product of inertia.
        self.Ixy_prime = self.L*self.d_prime[0]*self.d_prime[1]
        
    def calc_stress(self, load, s_direct, I, Ixy, H, Ip):

        # Stress in x and y from bending about z (torsion).
        # sx = -Mz/Ip*y
        # sy = Mz/Ip*x
        sxi = -1*load[5]/Ip*self.i_prime[1]
        syi = load[5]/Ip*self.i_prime[0]
        sxj = -1*load[5]/Ip*self.j_prime[1]
        syj = load[5]/Ip*self.j_prime[0]
        
        if Ixy != 0:
            # If Ixy = 0, the weld group centroidal axis are is aligned with the principal axis. H is probably
            # also 0.
    
            # Stress in z from bending about x and y. The following equation will give the stress about a 
            # set of arbitrary axes with moments and geometric properties relative to those axes. Look up interia 
            # tensor and how to use.
            # sz = -((My*Ix + Mx*Ixy)/(Ix*Iy - I^2_xy))*x + ((Mx*Iy + My*Ixy)/(Ix*Iy - I^2_xy))*y
            # The bending about the principal axes, using prinicpal axes properties, is
            # sz = Mx'y'/Ix - My'x'/Iy'
            stress_coeff = np.array([-1*(load[4]*I[1] + load[3]*Ixy), load[3]*I[0] + load[4]*Ixy, 0])/H
            szi = np.dot(self.i_prime, stress_coeff)
            szj = np.dot(self.j_prime, stress_coeff)
        
        elif I[0] == 0:
            # Weld group has no resistance to moment about y-axis.
            
            # stress in z direction from moments at start of weld line.
            # szi = -load[3]*self.i_prime[1]/I[1] + load[4]*self.i_prime[0]/I[0]
            # szj = -load[3]*self.j_prime[1]/I[1] + load[4]*self.j_prime[0]/I[0]

            if load[4] != 0: 
                # Note that ignoring the z-direction stress from My will effectively ignore My and may
                # lead the user to beleive the weld configuration effectively resists the applied loads.
                # Set high stress to indicate ineffectiveness.
                szi = szj = 9*10**6
            else:
                szi = -load[3]*self.i_prime[1]/I[1]
                szj = -load[3]*self.j_prime[1]/I[1]
            
        elif I[1] == 0:
            # Weld group has no resistance to moment about x-axis.

            if load[3] != 0: 
                szi = szj = 9*10**6
                
            else:
                szi = load[4]*self.i_prime[0]/I[0]
                szj = load[4]*self.j_prime[0]/I[0]
            
        # Stresses due to moments.
        self.smi = np.array([sxi, syi, szi])
        self.smj = np.array([sxj, syj, szj])
        
        # Combined stress in each direction from moment and direct shear.
        self.si = s_direct + self.smi
        self.sj = s_direct + self.smj
        
        # Total stress resultant.
        self.si_mag = np.linalg.norm(self.si)
        self.sj_mag = np.linalg.norm(self.sj)
        
            
class Weld_Group(object):
    
    def __init__(self, welds, loads, Fexx=70.0, phi=0.75):
        
        self.name = "Weld Group"
        self.welds = {w.name: w for w in welds}
        self.loads = {l.name: l for l in loads}
        self.Fexx = Fexx
        self.phi = phi

    def new_segment(self, weld):
        self.welds[weld.name] = weld

    def del_segment(self, name):
        del self.welds[name]
        
    def new_load(self, load):
        self.loads[load.name] = load

    def del_load(self, name):
        del self.loads[name]
        
    def analyze(self):
        
        # Weld group properties.
        self.A = sum([self.welds[k].L for k in self.welds])
        self.cg = sum([self.welds[k].Lc for k in self.welds])/self.A
        
        # Calculate centroidal properties.
        for k in self.welds:
            self.welds[k].calc_prime(self.cg)
        
        # Sum I about group centroid to get I for the weld group.
        self.I = sum([self.welds[k].I_prime for k in self.welds])
        # Product of inertia. Ixy = J.
        self.Ixy = sum([self.welds[k].Ixy_prime for k in self.welds])
        # Polar moment of inertia.
        self.Ip = self.I[0] + self.I[1]
        # Precalculate Ixx'*Iyy' - Ixy'^2 for use in unsymmetric bending equation.
        self.H = self.I[0]*self.I[1] - self.Ixy**2
        
        # Calculate centroidal properties of the load.
        for k in self.loads:
            self.loads[k].calc_cg(self.cg)
            
        # Aggregate centroidal loads.
        self.load = sum(self.loads[k].load_prime for k in self.loads)

        # Calculate direct stresses for weld group. s = [sx, sy, sz]
        self.s_direct = self.load[:3]/self.A
        
        # Calculate stresses at start and end coordinates for each weld.
        for k in self.welds:
            self.welds[k].calc_stress(self.load, self.s_direct, self.I, self.Ixy, self.H, self.Ip)
            
    def design(self):
        
        # Find max stress.
        self.s_max = max([max(self.welds[k].si_mag, self.welds[k].sj_mag) for k in self.welds])
        
        # phi_Rn = phi*(0.7071*w*F_W)
        # The weld shear strength varies with the load direction. Use F_W = 0.6*F_EXX for the entire weld
        # w = phi_Rn/(phi*0.7071*0.6*Fexx) for a fillet weld.
        self.FW = 0.6*self.Fexx
        self.Rn_req = self.s_max/self.phi
        self.throat_req = self.Rn_req/self.FW
        self.fillet = self.throat_req/0.7071 # sqrt(2)/2
        self.phi_Rn = self.phi*(0.7071*self.throat_req*self.FW)
        self.Rn = 0.7071*self.fillet*self.FW

if __name__ == "__main__": 
    
    # Example 0: Symmetrical bending of box shape.
    # group = Weld_Group([Weld(xj=6.0),
    #                     Weld(yj=6.0),
    #                     Weld(0, 6, 6, 6), 
    #                     Weld(6, 0, 6, 6)], 
    #                    [Weld_Load(vy=100, x=3.0, y=3.0)])
    
    # Example 1: Unsymmetrical bending of singly-symmetric shape ([ shape)
    # group = Weld_Group([Weld(xj=6.0),
    #                     Weld(yj=6.0),
    #                     Weld(0, 6, 6, 6)],
    #                    [Weld_Load(vy=100, x=5.0),
    #                     Weld_Load(pz=10.0, vx=15.0, vy=20.0, mx=16.0, my=15.0, mz=14.0, x=1.0, y=1.0, z=1.0)])
    
    # Example 2: Unsymmetrical bending occurs when a beam is loaded through a transverse plane not passing 
    # through the shear center.
    # group = Weld_Group([Weld(xj=6.0),
    #                     Weld(yj=6.0)], 
    #                    [Weld_Load(my=100.0, x=3.0)])
    
    # Example 3: Unstable Weld Pattern
    # group = Weld_Group([Weld(xj=6.0)], 
    #                    [Weld_Load(mx=100.0, x=5.0)])
    
    # Example 3: H is 0 if Ix = 0 or Iy = 0 and Ixy = 0. For each weld, the product of inertia, Ixy, is 0 when at least one 
    # of the coordinates of the centroid are at the weld group centroid. For each weld, Ix = 0 or Iy = 0 when the weld has no 
    # moment of inertia in the relevant direction.
    group = Weld_Group([Weld(name='1', xj=6.0)],
                       [Weld_Load(name='1', my=100, vy=100.0, x=1.0)])
    
    group.analyze()
    group.design()
    
    print('H', group.H)
    print('I', group.I)
    print('Ixy', group.Ixy)
    print('Iy*Ix - Ixy^2', group.I[0], group.I[1], group.Ixy)
    print('smax', group.s_max)
