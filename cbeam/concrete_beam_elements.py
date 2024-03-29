import numpy as np

class AppliedLoad(object):
    def __init__(self, parent, Mu=0.0, Vu=0.0):
        self.Mu = Mu
        self.Vu = Vu
        
    def calc(self, Mu, Vu):
        self.Mu = Mu
        self.Vu = Vu

class ConcMaterial(object):
    def __init__(self, parent, fc_prime=5000.0, strain_ult=0.003, Ec=3605000.0):
        self.parent = parent
        self.strain_ult = strain_ult
        self.fc_prime = fc_prime
        self.Ec = Ec
        self.beta1 = self.calc_beta1(fc_prime)
        self.fc = self.calc_fc(fc_prime)    # Whitney block stress.
        
    def calc_fc(self, fc_prime):
        fc = fc_prime*0.85
        
        return fc
    
    def calc_beta1(self, fc_prime):
        if fc_prime <= 4000.0:
            beta1 = 0.85
        elif fc_prime >= 8000.0:
            beta1 = 0.65
        else:
            beta1 = 0.85 - 0.00005*(fc_prime - 4000.0)
            
        return beta1
    
    def calc(self, fc_prime, strain_ult, Ec):
        self.fc_prime = fc_prime
        self.strain_ult = strain_ult
        self.Ec = Ec
        self.beta1 = self.calc_beta1(self.fc_prime)
        self.fc = self.calc_fc(fc_prime)

class LongReinf(object):
    def __init__(self, parent, name=1, size='9', qty=1.0, x=0.0, y=0.0, fy=60000.0, Es=29000000.0):
        self.parent = parent
        self.name = name
        self.size = size
        self.qty = qty
        self.D, self.A, self.weight = self.get_bar_props(size, qty)
        self.x = x
        self.y = y
        self.x_rc = None
        self.y_rc = None
        self.fy = fy
        self.Es = Es
        self.strain = None
        self.stress = None
        self.F = None
        self.Mz_rc = None
        self.My_rc = None
        self.Mz = None
        self.My = None
        
    def get_bar_props(self, size, qty):
        bar_lib = {'3': {'D': 0.375, 'A': 0.11, 'weight': 0.376},
                   '4': {'D':0.500, 'A': 0.20, 'weight': 0.668},
                   '5': {'D':0.625, 'A': 0.31, 'weight': 1.043},
                   '6': {'D':0.750, 'A': 0.44, 'weight': 1.502},
                   '7': {'D':0.875, 'A': 0.60, 'weight': 2.044},
                   '8': {'D':1.000, 'A': 0.79, 'weight': 2.67},
                   '9': {'D':1.128, 'A': 1.00, 'weight': 3.4},
                   '10': {'D':1.270, 'A': 1.27, 'weight': 4.303},
                   '11': {'D':1.410, 'A': 1.56, 'weight': 5.313},
                   '14': {'D':1.693, 'A': 2.25, 'weight': 7.65},
                   '18': {'D':2.257, 'A': 4.00, 'weight': 13.6}}
        D = bar_lib[size]['D']
        A = bar_lib[size]['A']
        weight = bar_lib[size]['weight']
        
        return D, A*qty, weight*qty
        
    def calc(self, size, qty, x, y, fy, Es):
        self.size = size
        self.qty = qty
        self.x = x
        self.y = y
        self.fy = fy
        self.Es = Es
        self.D, self.A, self.weight = self.get_bar_props(self.size, self.qty)
        
class XsectShape(object):
    def __init__(self, parent, shape_type="Rectangular", verts=None, b=12.0, h=24.0, bf=None, hf=None, cover=1.5):
        self.parent = parent
        self.cover = cover
        self.verts_nc = self.set_verts(shape_type, verts, b, h, bf, hf)
        self.A, self.z_gc, self.y_gc = self.calc_props_from_vertices(self.verts_nc)

    def set_verts(self, shape_type, verts, b, h, bf, hf):
        self.shape_type = shape_type
        self.verts = verts
        self.b = b
        self.h = h
        self.bf = bf
        self.hf = hf

        if shape_type == "Rectangular":
            verts = np.array([[0, 0], [b, 0], [b, h], [0, h]])
        if shape_type == "T-Beam":
            verts = np.array([[0, 0], 
                              [bf, 0], 
                              [bf, hf], 
                              [(bf+b)/2, hf], 
                              [(bf+b)/2, h], 
                              [(bf-b)/2, h], 
                              [(bf-b)/2, hf], 
                              [0, hf]])
        else: #shape_type == "Custom":
            verts = np.array(verts)

        self.min_y = min(verts[:,1])
        self.max_y = max(verts[:,1])
        self.min_x = min(verts[:,0])
        self.max_x = max(verts[:,0])
        self.max_h = abs(self.max_y - self.min_y)
        
        return verts
    
    def calc_props_from_vertices(self, native_verts):
        # Append the first row as the last row as well. The passed array should
        # be copied so that it is not altered. This copy process should work 
        # with Python or Numpy arrays.
        verts = np.array(native_verts)
        verts = np.vstack((verts, verts[0, :]))
        # Convert to a numpy array.
        # Compute the component cross-sectional properties.
        z0 = verts[:-1,0]       # z-coordinates from the first row to one row before the end.
        z1 = verts[1:,0]        # z-coordinates from the second row to the end.
        y0 = verts[:-1,1]
        y1 = verts[1:,1]
        Ai = np.multiply(z0, y1) - np.multiply(z1, y0)
        zci = np.multiply(z0 + z1, Ai)
        yci = np.multiply(y0 + y1, Ai)
        
        # Calculate the total cross-sectional properties.
        # If the vertices were entered clockwise the area will be negative.
        Ac = (1/2)*abs(sum(Ai))
        if Ac == 0.0:
            print('Area of shape is 0. Shape is likely degenerate.')
        else:
            # zc and yc will be represented as location relative to the 
            # reference coordinates.
            zc = (1/(6*Ac))*sum(zci)
            yc = (1/(6*Ac))*sum(yci)
            
        return (Ac, zc, yc)

    def calc(self, shape_type, verts, b, h, bf, hf):
        self.verts_nc = self.set_verts(shape_type, verts, b, h, bf, hf)
        
        # Calculate the analytical section properties of the concrete cross-section relative to the geometric centroid.
        # Coordinates of the geometric centroid relative to the native centroid.
        self.A, self.z_gc, self.y_gc = self.calc_props_from_vertices(self.verts_nc)
        
class Vertex(object):
    def __init__(self, x=0.0, y=0.0):
        self.x = x
        self.y = y
        