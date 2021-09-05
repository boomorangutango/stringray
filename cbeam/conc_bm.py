#%%
from concrete_beam_elements import LongReinf, XsectShape, ConcMaterial
import numpy as np
import plotly.graph_objects as go
#%%
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
    
    def new_long_reinf(self, size='9', qty=1.0, x=6.0, y=21.375, fy=60000.0, Es=29000000.0):
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
            # Everytime the lower bound of the stress block intersects the boundary of the 
            # cross-section, the boundary is either above, below, or at the N.A. If the i+1 
            # coordinate is below the bound, the next boundary vertices will be below the bound 
            # until another intersection occurs. If the i+1 coordinate is above the bound,
            # the next boundary vertices will be above the bound until another intersection 
            # occurs. If the i+1 coordinate is at the bound, the next coordinates may be above 
            # or below the bound; in this case the i+2 coordinate is needed to determine.
            # Determine if the bound intersects the boundary of the 
            # cross-section at between the current set of coordinates and the next set.
            x1, y1 = verts[i]
            x2, y2 = verts[i+1]
            # Note that this procedure assumes the coordinate system of the shape is defined with the 
            # origin at the compression set of the beam.
            if y1 <= a:
                # Current vertex is at or above bound.
                whitney_stress_region.append([x1, y1])
            if (a > y1 and a < y2) or (a < y1 and a > y2):
                # The bound is also between the current and the next vertex.
                # Add the point at the intersection of the 
                # N.A. and the line connecting the current and next 
                # vertices. We know the y-coordinate because the nuetral axis
                # depth is parallel to the y-axis. The x-coordinate of the 
                # intersection of the nuetral axis and the line segment is
                # as follows.

                # Catch the case of the intersection with a vertical line. This will be a very common case.
                if x1 == x2:
                    whitney_stress_region.append([x1, a])
                else:
                    x = (a - (y1 - (y2 - y1)/(x2 - x1)*x1))/((y2 - y1)/(x2 - x1))
                    # Record the new point in the array of vertices definining the bounded region.
                    whitney_stress_region.append([x, a])
                    
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

    def plot_ly(self):
        # Plot section, Whitney block, PNA.
        section_fig = go.Figure(layout=go.Layout(title="Section",
                                                 height=600, width=600,
                                                 yaxis=dict(title="y [in]", scaleanchor="x", scaleratio=1, 
                                                            autorange='reversed', dtick=2),
                                                 xaxis=dict(title="x [in]", dtick=2)))
        shape_verts = np.vstack((self.xsect.verts_nc, self.xsect.verts_nc[0, :]))
        block_verts = np.vstack((np.array(self.conc_stress_region), self.xsect.verts_nc[0, :]))
        section_fig.add_trace(go.Scatter(x=shape_verts[:,0], y=shape_verts[:,1], name="Shape"))
        section_fig.add_trace(go.Scatter(x=block_verts[:,0], y=block_verts[:,1], name="Whitney Block", fill='toself'))
        section_fig.add_trace(go.Scatter(x=[self.ecz], y=[self.ecy], name="Whitney Block CG", mode='markers', 
                                         marker_symbol='cross'))
        section_fig.add_trace(go.Scatter(x=[min(shape_verts[:,0]), max(shape_verts[:,0])], 
                                          y=[self.pna, self.pna], name="P.N.A.", mode='lines', 
                                          line=dict(dash='dash')))
        # Plot bars.
        for b in self.long_reinf:
            section_fig.add_shape(type="circle", xref="x", yref="y", x0=b.x-b.D/2, y0=b.y-b.D/2, 
                                  x1=b.x+b.D/2, y1=b.y+b.D/2)
            section_fig.add_trace(go.Scatter(x=[b.x], y=[b.y], mode='markers', marker_symbol='cross', 
                                             name=str(b.qty)+"#"+b.size))

        # Plot strain.
        strain_fig = go.Figure(layout=go.Layout(title="Strain", 
                                                height=600, width=600,
                                                yaxis=dict(title="y [in]", autorange='reversed', dtick=2),
                                                xaxis=dict(title="strain [in/in]")))
        # Plot vertical line representing the concrete section.
        strain_fig.add_trace(go.Scatter(x=[0, 0], y=[self.xsect.max_y, self.xsect.min_y], name="section"))
        # Plot PNA.
        strain_fig.add_trace(go.Scatter(x=[0], y=[self.pna], name="P.N.A.", mode='markers', marker_symbol='line-ew',
                                        marker_line_width=2, marker_size=20))
        # Plot strain curve.
        strain_fig.add_trace(go.Scatter(x=[self.material.strain_ult, -self.material.strain_ult*(self.xsect.max_y/self.pna - 1)],
                                        y=[self.xsect.min_y, self.xsect.max_y], name="strain curve", mode='lines', 
                                        line=dict(dash='dash')))
        # Plot and label ultimate concrete strain.
        strain_fig.add_trace(go.Scatter(x=[self.material.strain_ult], y=[self.xsect.min_y], 
                                        mode='markers', name="eps_cu", marker_symbol='cross'))
        # Plot bar stain.
        for b in self.long_reinf:
            strain_fig.add_trace(go.Scatter(x=[b.strain], y=[b.y], mode='markers', marker_symbol='cross', name=str(b.qty)+"#"+b.size))
            
        # Plot stress.
        stress_fig = go.Figure(layout=go.Layout(title="Stress", 
                                                height=600, width=600,
                                                yaxis=dict(title="y [in]", autorange='reversed', dtick=2),
                                                xaxis=dict(title="stress [lb/in^2]")))
        # Plot vertical line representing the concrete section.
        stress_fig.add_trace(go.Scatter(x=[0, 0], y=[self.xsect.max_y, self.xsect.min_y], name="section"))
        # Plot PNA.
        stress_fig.add_trace(go.Scatter(x=[0], y=[self.pna], name="P.N.A.", mode='markers', marker_symbol='line-ew',
                                        marker_line_width=2, marker_size=20))
        # Plot concrete stress.
        sigma_c_curve = np.array([[self.material.fc, self.xsect.min_y], [self.material.fc, self.a]])
        stress_fig.add_trace(go.Scatter(x=sigma_c_curve[:,0], y=sigma_c_curve[:,1], 
                                        fill='tozerox', name="concrete"))
        # Plot bar stress.
        for b in self.long_reinf:
            sigma_s_curve = np.array([[b.stress, b.y+b.D/2], [b.stress, b.y-b.D/2]])
            stress_fig.add_trace(go.Scatter(x=sigma_s_curve[:,0], y=sigma_s_curve[:,1], 
                                            fill='tozerox', name=str(b.qty)+"#"+b.size))
        
        # Plot force.
        force_fig = go.Figure(layout=go.Layout(title="Force", 
                                               height=600, width=600,
                                               yaxis=dict(title="y [in]", autorange='reversed', dtick=2),
                                               xaxis=dict(title="force [lb]")))
        # Plot vertical line representing the concrete section.
        force_fig.add_trace(go.Scatter(x=[0, 0], y=[self.xsect.max_y, self.xsect.min_y], name="section"))
        # Plot concrete force.
        force_fig.add_trace(go.Scatter(x=[0, self.Fc], y=[self.ecy, self.ecy], name="concrete force"))
        # Plot bar force.
        for b in self.long_reinf:
            force_fig.add_traces(go.Scatter(x=[0, b.F], y=[b.y, b.y], name=str(b.qty)+"#"+b.size))
        
        return [section_fig, strain_fig, stress_fig, force_fig]
#%% __main__

if __name__ == "__main__":
    beam = ConcreteBeam()
    beam.new_xsect_shape(b=12.0, h=24.0, bf=0.0, hf=0)
    beam.new_long_reinf(size='9', qty=1.0, x=6, y=21.0)
    # beam.new_long_reinf(size='10', qty=3.0, x=20, y=21.0)
    # beam.new_xsect_shape(verts=[[3,0],
    #                             [9,0],
    #                             [12, 24],
    #                             [0, 24]])
    beam.new_material()
    beam.run_calc()
    print(beam.phi_Mnx/1000/12)    