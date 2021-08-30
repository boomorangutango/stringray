from concrete_beam_elements import LongReinf, XsectShape, ConcMaterial
import numpy as np
import matplotlib.pyplot as plt

class ConcreteBeam(object):
    def __init__(self):
        self.long_reinf = []
        
        # Nominal capacities of the concrete beam section.
        self.Pn = 0.0
        self.Mnx = 0.0
        self.Mny = 0.0
        self.Mnx_rc = 0.0
        self.Mny_rc = 0.0
        
        # Analytical properties related to the concrete.
        self.Fc = None
        self.Mz_conc_rc = None
        self.My_conc_rc = None
        self.Mz_conc = None
        self.My_conc = None
        # Compression block area and coordinates.
        self.Ac = None
        
        # Initialize the coordinates of the P.C. of the beam section.
        self.y_pc = None
        self.x_pc = None
    
    def new_long_reinf(self, size='6', qty=1.0, x=0.0, y=0.0, fy=60000.0, Es=29000000.0):
        bar = LongReinf(self, size, qty, x, y, fy, Es)
        self.long_reinf.append(bar)
        return bar
    
    def new_xsect_shape(self, verts=None, b=12.0, h=24.0, bf=0.0, hf=0.0, cover=1.5):
        kwargs = {'verts': verts, 'b': b, 'h': h, 'hf': hf, 'bf': bf, 'cover': cover}
        self.xsect = XsectShape(self, **kwargs)
        return self.xsect
    
    def new_material(self, fc_prime=5000.0, strain_ult=0.003, Ec=3605000.0):
        self.material = ConcMaterial(self, fc_prime, strain_ult, Ec)
        return self.material

    def calc_phi(self, long_reinf):
        # Calculate the strength reduction factor. This factor is depedent on
        # the strain in the tension steel.
        min_strain = min(bar.strain for bar in long_reinf)
        if min_strain > -0.002:
            phi = 0.0
        elif min_strain < -0.005:
            phi = 0.9
        else: 
            phi = 0.65 + (abs(min_strain) - 0.002)*(250/3)
        return phi

    # def calc_PC_coords(self):
    #     # ! This method doesn't make sense to me. The location of the plastic centroid should 
    #     # be relative to the cracked section. This method seems to assume every bar is 
    #     # yeilding in compression and the whole section is in compression.
        
    #     '''
    #     Calculate the plastic centroid.
        
    #     The plastic centroid is the centriod of the section when the entire portion of the 
    #     concrete cross-section experiencing compression is at maximum compressive capacity.
    #     The plastic nuetral axis is located at the geometric centroid of the transformed 
    #     section. The transformed section may be found by equating the area of steel 
    #     reinforcement to an equivalent area of concrete. Since the strain in the steel 
    #     and concrete in the section is assumed to be equal, the area of concrete equivalent 
    #     to an area of steel may be found by the modular ratio eta = Es/Ec or by the stress
    #     ratio fy/fc_prime. The equivalent area is as follows.
    #     As_eq = As*(fy/fc_prime - 1)
    #     '''
    #     # Calculate the analytical section properties relative to the geometric
    #     # centroid.
        
    #     # Calculate the first moments of area (statical moments) relative to 
    #     # the native centroid. The native centroid is or origin of the coordinate system
    #     # in which the cross-section was defined. The first moment of area is calculated as 
    #     # the area of a region multiplied by the distance from the geometric centroid of the 
    #     # region to an axis.
    #     Qz_nc = self.xsect.A*self.y_gc
    #     Qy_nc = self.xsect.A*self.z_gc
        
    #     # Calculate the analytical section properties of the cross-section relative to 
    #     # the P.C. by including the steel reinforcement as an equivalent concrete area.
    #     for bar in self.long_reinf:
    #         # Equivalent concrete area.
    #         # ! This doesn't seem right. Doesn't this assume the bar always yeilds?
    #         bar.A_eq = bar.A*(bar.fy/self.fc_prime - 1)
            
    #         # Calculate the first moments of area (statical moments) of the 
    #         # transformed reinforcement.
    #         bar.Qy_eq_nc = bar.coords[0]*bar.A_eq
    #         bar.Qz_eq_nc = bar.coords[1]*bar.A_eq
            
    #     # Get the total first moment of area of the steel relative to the 
    #     # native coordinate system and get the total area of the steel.
    #     tot_Qz_eq_steel_nc = sum([b.Qz_eq_nc for b in self.long_reinf])
    #     tot_Qy_eq_steel_nc = sum([b.Qy_eq_nc for b in self.long_reinf])
    #     tot_A_eq_steel = sum([b.A_eq for b in self.long_reinf])
        
    #     # Find the coordinates of the P.C. of the cross-section relative to the native 
    #     # coordinate system.
    #     self.y_pc = (Qz_nc + tot_Qz_eq_steel_nc)/(self.xsect.A + tot_A_eq_steel)
    #     self.z_pc = (Qy_nc + tot_Qy_eq_steel_nc)/(self.xsect.A + tot_A_eq_steel)
    
    #     # Translate the concrete cross-section vertices to be relative to the P.C.
    #     self.conc_verts_pc = self.xsect.verts_nc - np.array([self.z_pc, self.y_pc])
        
    #     # Translate the bar coordinates.
    #     for bar in self.long_reinf:
    #         bar.coords_pc = bar.coords - np.array([self.z_pc, self.y_pc])    

    def transform_to_rotated_coords(self, rotation=0.0):

        # Form the rotation matrix about z-axis from an angle given in degrees.
        # This rotation matrix pre-multiplies a column vector to give vector 
        # components in the rotated coordinate system. Positive rotation is 
        # assumed to occur in the counterclockwise direction when looking 
        # towards the origin (right hand rule). This means that positive 
        # bending about the y-axis requires flipping the sign of the moment 
        # about the y-axis.
        cos_x = np.cos(rotation)
        sin_x = np.sin(rotation)
        # Rotation matrix about x-axis for 2D.
        self.Rx_to_rc = np.array([[cos_x, -sin_x],
                                  [sin_x, cos_x]])

        # Using the symmetry of the cosine function and the antisymmetry of the 
        # sine function, the rotatioin matrix to return to the non-rotated axis
        # is as follows.
        self.Rx_from_rc = np.array([[cos_x, sin_x],
                                    [-sin_x, cos_x]])
                                    
        # Rotate the vertices of the concrete section with the rotation matrix.
        # The vertices are listed as row vectors so the rotation matrix must 
        # post-multiply with the transpose.
        # conc_verts_rc = np.dot(self.conc_verts_pc, self.Rx_to_rc.T)
        conc_verts_rc = np.dot(self.xsect.verts_nc, self.Rx_to_rc.T)
        
        # Rotate the bar coordinates.
        for bar in self.long_reinf:
            # bar.coords_rc = np.dot(bar.coords_pc, self.Rx_to_rc.T)
            bar.x_rc, bar.y_rc = np.dot(np.array([bar.x, bar.y]).T, self.Rx_to_rc.T)
        
        return conc_verts_rc

    def get_whitney_block(self, conc_verts_rc, pna):
        
        # The area of concrete in compression is assumed to be equal to 
        # the area of the Whitney compression block. The Whitney compression 
        # block is assumed to start at the most extreme compression fiber 
        # and extend a distance of a = beta1*c  toward the most extreme tension fiber.
        # self.a = self.max_y_rc - self.beta1*(self.max_y_rc - c)
        a = self.material.beta1*pna
        
        # The area of concrete in compression is assumed to be equal to 
        # the area of the Whitney compression block. The Whitney compression 
        # block is assumed to start at the most extreme compression fiber 
        # and extend a distance of beta1*c  toward the most extreme tension fiber.
        # self.a = self.max_y_rc - self.beta1*(self.max_y_rc - c)

        # Solve for the location of the N.A., given the depth of the Whitney 
        # stress block.
        # c = (a - max_y_rc)/self.beta1 + max_y_rc
        # self.c[ri].append(c)
                
        # Find the intersections of the lower boundary of the stress 
        # block with the concrete cross-section. We may want to implement
        # line to line segment intersection detection techniques.
        whitney_stress_region = []
        
        # Create a set of vertices for iteration with the first vertex 
        # appended to the end.
        verts = np.append(conc_verts_rc, [conc_verts_rc[0,:]], axis=0)
     
        for i in range(len(verts)-1):
            # Everytime the lower bound of the stress block intersects
            # the boundary of the cross-section, the boundary is 
            # either above, below, or at the 
            # N.A. If the i+1 coordinate is below the bound, the next 
            # boundary vertices will be below the bound until another 
            # intersection occurs. If the i+1 coordinate is above the bound,
            # the next boundary vertices will be above the bound until 
            # another intersection occurs. If the i+1 coordinate is at the
            # bound, the next coordinates may be above or below the bound; in
            # this case the i+2 coordinate is needed to determine.
            # Determine if the bound intersects the boundary of the 
            # cross-section at between the current set of coordinates and the 
            # next set.
            x1, y1 = verts[i]
            x2, y2 = verts[i+1]
            # y1 = verts[i, 1]
            # y2 = verts[i+1, 1]
            # Note that this procedure assumes the coordinate system of the shape is defined with the 
            # origin at the compression set of the beam.
            if y1 <= a:
                # Current vertex is at or above bound.
                whitney_stress_region.append([x1, y1])
            # if a == y1:
            #     # The bound is at the current vertex. Add the current 
            #     # vertex to the set.
            #     whitney_stress_region.append([x1, y1])
            if (a > y1 and a < y2) or (a < y1 and a > y2):
                # The bound is also between the current and the next vertex.
                # Add the point at the intersection of the 
                # N.A. and the line connecting the current and next 
                # vertices. We know the y-coordinate because the nuetral axis
                # depth is parallel to the y-axis. The x-coordinate of the 
                # intersection of the nuetral axis and the line segment is
                # as follows.

                # Catch the case of the intersection with a vertical line.
                # This will be a very common case.
                if x1 == x2:
                    whitney_stress_region.append([x1, a])
                else:
                    x = (a - (y1 - (y2 - y1)/(x2 - x1)*x1))/((y2 - y1)/(x2 - x1))
                    # Record the new point in the array of vertices definining 
                    # the bounded region.
                    whitney_stress_region.append([x, a])
            # elif a > y1 and a > y2:
            #     # The bound is below both the current and next vertex.
            #     # Add the current vertex to the set.
            #     whitney_stress_region.append([x1, y1])
            # elif a < y1 and a < y2:
            #     # The bound is higher both the current and next vertex. Neither vertex is in the set.
            #     continue
            # elif a < y1 and a > y2:
            #     # The bound is between the current and the next vertex with 
            #     # the current vertex above the N.A. and the next vertex 
            #     # below the bound Add the current vertex and the point at 
            #     # the intersection of the bound and the line connecting the 
            #     # current and next vertices.
                
            #     # Add the current vertex to the set.
            #     whitney_stress_region.append([x1, y1])
                
            #     # Add the point at the intersection of the current vertex
            #     # and the next vertex.
            #     if x1 == x2:
            #         # Catch the case of the intersection with a vertical line.
            #         # This will be a very common case.
            #         whitney_stress_region.append([x1, a])
            #     else:
            #         # Find the x-coordinate of the intersection.
            #         x = (a - (y1 - (y2 - y1)/(x2 - x1)*x1))/((y2 - y1)/(x2 - x1))
            #         # Record the new point in the array of vertices definining 
            #         # the bounded region.
            #         whitney_stress_region.append([x, a])
                    
        return whitney_stress_region

    def calc_moments(self):
        
        # Accumulate the equivalent concrete for the area of the bars in compression. The force 
        # and moments are subtracted from the overall force in the concrete in compression.
        accum_Mz = 0.0
        accum_My = 0.0
        
        for bar in self.long_reinf:
            
            # Subtract the force of concrete occupied by the area of steel in compression.
            if bar.strain > 0:
                bar.Fc_eq = bar.A*0.85*self.material.fc_prime
                accum_Mz += bar.y_rc*bar.Fc_eq
                accum_My += -bar.x_rc*bar.Fc_eq
            
            # The moment contribution of each bar is the force in the bar multiplied by the eccentricity.
            bar.Mz_rc = bar.F*bar.y_rc
            bar.My_rc = -bar.F*bar.x_rc
            self.Mnx_rc += bar.Mz_rc
            self.Mny_rc += bar.My_rc

            # Rotate the steel moments back to the original coordinate system with the rotation matrix.
            Ms = np.dot(self.Rx_from_rc, [[bar.Mz_rc], [bar.My_rc]])
            bar.Mz = Ms[0,0]
            bar.My = Ms[1,0]
            # Add the steel moment to the section moments.
            self.Mnx += Ms[0,0]
            self.Mny += Ms[1,0]
        
        # Moment of the concrete in compression about the origin.
        Mcx_rc = self.Fc*self.ecy - accum_Mz
        Mcy_rc = -self.Fc*self.ecz - accum_My
        self.Mz_conc_rc = Mcx_rc
        self.My_conc_rc = Mcy_rc
        self.Mnx_rc += Mcx_rc
        self.Mny_rc += Mcy_rc

        # Rotate the moments back to native coordinates.
        Mc = np.dot(self.Rx_from_rc, [[Mcx_rc], [Mcy_rc]])
        self.Mz_conc = Mc[0,0]
        self.My_conc = Mc[1,0]
        self.Mnx += Mc[0,0]
        self.Mny += Mc[1,0]
        
    def calc_axial(self, pna):
        # Calculate the properties of the Whitney stress block. The center of force in the concrete 
        # is assumed to be at the centroid of the Whitney compression block and the moment arm of 
        # that force is measured to the P.C.
        self.conc_stress_region = self.get_whitney_block(self.conc_verts_rc, pna)
        self.Ac, self.ecz, self.ecy = self.xsect.calc_props_from_vertices(self.conc_stress_region)
        
        # Accumulate the equivalent concrete for the area of the bars in compression. The force 
        # and moments are subtracted from the overall force in the concrete in compression.
        Pn = 0.0
        Asc = 0.0
        
        for bar in self.long_reinf:
            # Find the strain in the bar and append the strain curve. From 
            # linear strain curve, the proportion of the distance from the 
            # N.A. to the extreme compression fiber to the strain at the 
            # extreme compression fiber is equal to to the proportion of 
            # the distance to the extreme tension fiber to the strain in 
            # the extreme tension fiber.
            # Note that the spreadsheet that accompanies SP-17 sets the bar
            # strain rather than the depth of the N.A. This is a convenient 
            # method for determining the critical points for zero stress and 0.5*fy.
            bar.strain = -self.material.strain_ult*(bar.y_rc/pna - 1)
            
            # The stress in the steel is found from the strain. 
            # If the strain in the bar is more than the yield strain, the 
            # stress in the bar is equal to the yield stress.
            # If the strain is less than the yield strain, the stress 
            # is proportional to the strain through the elasic modulus. 
            bar.stress = max(min(bar.fy, bar.Es*bar.strain), -bar.fy)

            # The force in each bar is equal to the stress in the bar 
            # multiplied by the bar area.
            bar.F = bar.A*bar.stress
            Pn += bar.F

            # Test if the bar strain is compressive (assumed to be positive).
            # If so, account for the missing area of concrete occupied by the
            # bar in both the axial force and the moment calculations. Note
            # that many programs substract the force that would have occurred 
            # in the concrete from the bar force to make this adjustment. 
            # This program keeps the force that would have occured in the 
            # concrete separate from the bar force to ensure the contribution 
            # to moment and axial force from the individual bar is recorded 
            # accurately. The method employed here creates some error because 
            # the concrete missing at the location of the bar is not accounted
            # for in the calculation of the moment arm for the resultant force 
            # in the concrete.
            if bar.strain > 0:
                Asc += bar.A

        # Calculate the compressive force in the concrete.
        self.Fc = 0.85*self.material.fc_prime*(self.Ac - Asc)
        Pn += self.Fc

        return Pn
        
    def find_pna(self, rotation=0.0):
        
        # Calculate the section properties for the rotated section.
        self.conc_verts_rc = self.transform_to_rotated_coords(rotation)
        
        # Initialize the location of the lower bound of the Whitney stress block at y = 0.
        top = 0
        bot = max(self.conc_verts_rc[:, 1])
        pna = (top + bot)/2
        Pn = self.calc_axial(pna)
        self.log = []
        i = 0
        while i < 1000 and abs(Pn) > 0.0001:
            # Check for force equalibrium in the section. The maximum moment capacity occurs where 
            # the axial capacity is zero.
            if Pn > 0:    # Compression is greater than tension. Move c up.
                bot = pna
                pna = (top + pna)/2
            else:    # Tension is greater than compression. Move c down.
                top = pna
                pna = (pna + bot)/2
            Pn = self.calc_axial(pna)
            i += 1
        self.Pn = Pn
        self.pna = pna
        self.a = self.material.beta1*pna
    
    def run_calc(self):
        self.find_pna()
        self.calc_moments()
        self.phi = self.calc_phi(self.long_reinf)
        self.phi_Mnx = self.phi*self.Mnx

    def plot(self):
        
        plt.close("all")
        
        # Plot section.
        fig1, ax1 = plt.subplots()
        ax1.set_title('Section')
        ax1.set_xlabel('x [in]')
        ax1.set_ylabel('y [in]')
        ax1.invert_yaxis()

        # Plot section shape.
        shape_verts = np.vstack((self.xsect.verts_nc, self.xsect.verts_nc[0, :]))
        ax1.plot(shape_verts[:,0], shape_verts[:,1], marker='.', color='green', label='shape')
        
        # Plot Whitney block.
        region = np.array(self.conc_stress_region)
        ax1.fill(region[:,0], region[:,1], color='red', alpha=0.5, label='Whitney block')
        # Plot and annotate centroid.
        ax1.plot(self.ecz, self.ecy, marker='+', color='black')
        ax1.annotate("("+str(round(self.ecz,2))+", "+str(round(self.ecy,2))+")", 
                     (self.ecz, self.ecy), color='black', ha='center', va='bottom')
        
        # Plot and annotate bars.
        for b in self.long_reinf:
            bar_plot = plt.Circle((b.x, b.y), b.D/2, color='blue', alpha=0.5)
            ax1.add_patch(bar_plot)
            ax1.plot(b.x, b.y, marker='+', color='black')
            ax1.annotate(str(b.qty)+"#"+b.size+" ("+str(b.x)+", "+str(b.y)+")", (b.x, b.y), 
                         textcoords='offset points', xytext=(0,5),
                         color='black', ha='center', va='bottom')

        # Plot PNA.
        ax1.hlines(self.pna, xmin=min(shape_verts[:,0]) , xmax=max(shape_verts[:,0]),
                   colors='grey', linestyles='dashed', label="N.A.")

        # Plot and label points of interest.
        plot_points = np.unique(np.vstack((self.xsect.verts_nc, self.conc_stress_region)), axis=0)
        ax1.scatter(plot_points[:,0], plot_points[:,1], marker='.', color='black')
        lbls = ["("+str(round(plot_points[i,0],2))+", "+str(round(plot_points[i,1],2))+")" 
                for i in range(len(plot_points))]
        for i, l in enumerate(lbls):
            ax1.annotate(l, plot_points[i], textcoords='offset points', xytext=(0,5),
                         color='black', ha='center', va='bottom')

        # Set section plot properties.
        ax1.legend(loc='lower left')
        ax1.set_aspect(aspect='equal', adjustable='datalim', share=True)

        # Plot strain.
        fig2, ax2 = plt.subplots()
        ax2.set_title('Strain')
        ax2.set_xlabel('eps [in/in]')
        ax2.set_ylabel('y [in]')
        ax2.invert_yaxis()
        
        # Plot strain curve.
        ax2.plot([self.material.strain_ult, -self.material.strain_ult*(self.xsect.max_y/self.pna - 1)],
                 [self.xsect.min_y, self.xsect.max_y], linestyle='dotted', color='grey', label='strain curve')

        # Plot and label ultimate concrete strain. 
        ax2.scatter(self.material.strain_ult, self.xsect.min_y, marker='s', color='red', label='eps_cu')
        ax2.annotate("("+str(self.material.strain_ult)+", "+str(self.xsect.min_y)+")", 
                     (self.material.strain_ult, self.xsect.min_y), 
                     textcoords='offset points', xytext=(0,5),
                     color='red', ha='center', va='bottom')
        # Plot and label bar stain.
        for b in self.long_reinf:
            ax2.scatter(b.strain, b.y, marker='o', color='blue')
            ax2.annotate("("+str(round(b.strain,4))+", "+str(b.y)+")", 
                         (b.strain, b.y), 
                         textcoords='offset points', xytext=(0,5),
                         color='blue', ha='center', va='bottom')
            
        ax2.axvline(x=0, color='k')
        ax2.grid(True, which='both')
        ax2.legend(loc='upper left')

        # Plot stress.
        fig3, ax3 = plt.subplots()
        ax3.set_title('Stress')
        ax3.set_xlabel('eps [lb/in^2]')
        ax3.set_ylabel('y [in]')
        ax3.invert_yaxis()
        
        # Plot and label concrete stress. 
        fc = self.material.fc_prime*0.85
        sigma_c_curve = np.array([[fc, self.xsect.min_y], [fc, self.a]])
        ax3.plot(sigma_c_curve[:,0], sigma_c_curve[:,1], marker='s', color='red', label='sigma_c')
        ax3.annotate("("+str(self.material.strain_ult)+", "+str(round(self.ecy,2))+")", 
                     (self.material.strain_ult, self.ecy), 
                     textcoords='offset points', xytext=(5,0),
                     color='red', ha='right', va='center')        
        
        # Plot and label bar stress.
        for b in self.long_reinf:
            sigma_s_curve = np.array([[b.stress, b.y+b.D/2], [b.stress, b.y-b.D/2]])
            ax3.plot(sigma_s_curve[:,0], sigma_s_curve[:,1], marker='o', color='blue')
            ax3.annotate("("+str(round(b.stress,2))+", "+str(b.y)+")", 
                         (b.stress, b.y), 
                         textcoords='offset points', xytext=(5,0),
                         color='blue', ha='left', va='center')
            
        ax3.axvline(x=0, color='k')
        ax3.grid(True, which='both')
        ax3.legend(loc='upper left')
        
        # Plot force.
        fig4, ax4 = plt.subplots()
        ax4.set_title('Force')
        ax4.set_xlabel('F [lb]')
        ax4.set_ylabel('y [in]')
        ax4.invert_yaxis()
        
        # Plot and label concrete stress.
        ax4.arrow(0, self.ecy, self.Fc, 0, color='red', width=0.1, 
                  head_width=self.xsect.max_h/48, head_length=self.Fc/15, 
                  length_includes_head=True, label='Fc')
        ax4.annotate("("+str(round(self.Fc,2))+", "+str(round(self.ecy,2))+")", 
                     (self.Fc, self.ecy), 
                     textcoords='offset points', xytext=(-10,-5),
                     color='red', ha='right', va='top')        

        # Plot and label bar stress.
        for b in self.long_reinf:
            ax4.arrow(0, b.y, b.F, 0, color='blue', width=0.1, 
                      head_width=self.xsect.max_h/48, head_length=self.Fc/15, 
                      length_includes_head=True)
            ax4.annotate("("+str(round(b.F,2))+", "+str(b.y)+")", 
                          (b.F, b.y), 
                          textcoords='offset points', xytext=(10,5),
                          color='blue', ha='left', va='bottom')

        ax4.axvline(x=0, color='k')
        ax4.grid(True, which='both')
        
        plt.show()

        return
        
#%% __main__

if __name__ == "__main__":
    beam = ConcreteBeam()
    beam.new_xsect_shape(b=16.0, h=24.0, bf=0.0, hf=0)
    beam.new_long_reinf(size='10', qty=3.0, x=10, y=10.0)
    beam.new_long_reinf(size='10', qty=3.0, x=20, y=21.0)
    # beam.new_xsect_shape(verts=[[3,0],
    #                             [9,0],
    #                             [12, 24],
    #                             [0, 24]])
    beam.new_material()
    beam.run_calc()
    print(beam.phi_Mnx/1000/12)
    beam.plot()
