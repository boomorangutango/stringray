#%% Imports
import ipywidgets as wg
import pandas as pd

#%% ConcBeamLayout
class ConcBeamLayout(object):
    
    def __init__(self, model=None):
        self.model = model
        self.long_reinf = []
        self.output_box = self.build_output_box()
        self.plot_box = self.build_plot_box()
        self.input_box = self.build_input_box()

    def build_plot_box(self):
        # Create tab layout for plotting. Setup plot output widgets.
        self.section_out = wg.Output()
        self.strain_out = wg.Output()
        self.stress_out = wg.Output()
        self.force_out = wg.Output()
        plot = wg.Tab(children=[self.section_out, self.strain_out, self.stress_out, self.force_out])
        plot.set_title(0, 'Section')
        plot.set_title(1, 'Strain')
        plot.set_title(2, 'Stress')
        plot.set_title(3, 'Force')
        return plot
    
    def build_output_box(self):
        simple = wg.Output()
        detailed = wg.Output()
        output = wg.Tab(children=[simple, detailed])
        output.set_title(0, 'Simple')
        output.set_title(1, 'Detailed')
        return output

    def build_input_box(self):
        
        self.applied_load = self.new_applied_load()
        self.xsect_shape = self.new_xsect_shape()
        self.material = self.new_material()

        calc_btn = wg.Button(description='Run Calc', layout=wg.Layout(width='auto'))
        calc_btn.on_click(self.calc_plot)
        
        # Load input box. 
        applied_load_input = wg.VBox([self.applied_load.Mu, self.applied_load.Vu])
        
        self.shape_input = wg.VBox([self.xsect_shape.shape_type, self.xsect_shape.b, self.xsect_shape.h])
#         # Shape inputs box.
#         if self.xsect_shape.shape_type.value == "Rectangular":
            
#         elif self.xsect_shape.shape_type.value == "T-Beam":
#             self.shape_input = wg.VBox([self.xsect_shape.shape_type, self.xsect_shape.b, self.xsect_shape.h, 
#                                         self.xsect_shape.bf, self.xsect_shape.hf])
#         else: # Assumed "Custom"
#             self.shape_input = wg.VBox([self.xsect_shape.shape_type, self.xsect_shape.new_vert_btn, self.xsect_shape.verts_tbl])
            
        # Material inputs box.
        material_input = wg.VBox([self.material.fc_prime, self.material.strain_ult, self.material.Ec])
        
        # Create longitudinal steel new layer button and table.
        new_long_reinf_btn = wg.Button(description='New Bar Layer')
        new_long_reinf_btn.on_click(self.new_long_reinf)
        long_reinf_tbl = self.make_long_reinf_tbl()
        steel_input = wg.VBox([new_long_reinf_btn, long_reinf_tbl])
        
        # Build the input layout.
        applied_load_input_box = wg.Accordion(children=[applied_load_input])
        applied_load_input_box.set_title(0, 'Loads')
        shape_input_box = wg.Accordion(children=[self.shape_input])
        shape_input_box.set_title(0, 'Shape')
        material_input_box = wg.Accordion(children=[material_input])
        material_input_box.set_title(0, 'Material')
        steel_input_box = wg.Accordion(children=[steel_input])
        steel_input_box.set_title(0, 'Reinforcement')
        # Consolidate all input into one box.
        input_box = wg.VBox([applied_load_input_box, shape_input_box, material_input_box, 
                             steel_input_box, calc_btn])

        # The first plot request is run when this first bar is created.
        self.new_long_reinf(None)
    
        return input_box
        
    def calc_plot(self, b):
        
        # Calc all.
        self.xsect_shape.calc()
        self.material.calc()
        for b in self.long_reinf:
            b.calc()
        self.model.calc()
        # Generate plots.
        self.plot_all()
        # Generaate output.
        self.create_results()        
        
    def plot_all(self):
        # Clear existing plots.
        for plot in self.plot_box.children:
            plot.clear_output()
        
        # Create new plots.
        section_plot = self.model.plot_xsect()
        self.model.plot_section_analysis()
        strain_plot = self.model.plot_strain()
        stress_plot = self.model.plot_stress()
        force_plot = self.model.plot_force()
        
        # Show new plots.
        with self.section_out:
            section_plot.show()
        with self.strain_out:
            strain_plot.show()
        with self.stress_out:
            stress_plot.show()
        with self.force_out:
            force_plot.show()

    def create_results(self):
        
        for box in self.output_box.children:
            box.clear_output()
        
        # Create output.
        bar_cols = ['size', 'qty', 'x', 'y', 'strain', 'stress', 'F' ]
        bar_data = {}
        for b in self.model.long_reinf:
            bar_data[b.name] = [b.size, b.qty, b.x, b.y, b.strain, b.stress, b.F]
        bars_df = pd.DataFrame.from_dict(bar_data, orient='index', columns=bar_cols)
        
        with self.output_box.children[1]:
            display(wg.HTMLMath("<h3> Bars </h3>"))
            display(bars_df)
            
            print('c =', self.model.c, 'in')
            print('a =', self.model.a, 'in')
            print('beta1 =', self.model.material.beta1)
            print('phi =', self.model.phi)
            print('Mn =', self.model.Mnx, 'lb*in')
            print('phi*Mn =', self.model.phi_Mnx, 'lb*in')
            print('phi*Mn =', self.model.phi_Mnx/12/1000, 'k*ft')

    def new_applied_load(self):
        model = self.model.new_applied_load()
        layout = AppliedLoadLayout(self, model)
        return layout
        
    def new_xsect_shape(self):
        model = self.model.new_xsect_shape()
        layout = XsectShapeLayout(self, model)
        return layout
    
    def new_material(self):
        model = self.model.new_material()
        layout = ConcMaterialLayout(self, model)
        return layout
    
    def new_long_reinf(self, b):
        # Create the new bar layer model and layout. Give the model to the new bar layer layout 
        # and add the layout to the concrete beam layout. 
        model = self.model.new_long_reinf()
        layout = LongReinfLayout(self, model)
        self.long_reinf.append(layout)
        self.rebuild_long_reinf_tbl()
        self.calc_plot(None)
        
    def delete_long_reinf(self, long_reinf):
        # Remove the bar from the list of bars in both the model and layout and rebuild the bar layer table.
        self.model.long_reinf.remove(long_reinf.model)
        self.long_reinf.remove(long_reinf)
        self.rebuild_long_reinf_tbl()
        
    def rebuild_long_reinf_tbl(self):
        # Distribute bar info into columns.
        self.name_box.children = tuple([self.name_box.children[0]] + [b.name for b in self.long_reinf])
        self.size_box.children = tuple([self.size_box.children[0]] + [b.size for b in self.long_reinf])
        self.qty_box.children = tuple([self.qty_box.children[0]] + [b.qty for b in self.long_reinf])
        self.x_box.children = tuple([self.x_box.children[0]] + [b.x for b in self.long_reinf])
        self.y_box.children = tuple([ self.y_box.children[0]] + [b.y for b in self.long_reinf])
        self.fy_box.children = tuple([self.fy_box.children[0]] + [b.fy for b in self.long_reinf])
        self.Es_box.children = tuple([self.Es_box.children[0]] + [b.Es for b in self.long_reinf])
        self.delete_box.children = tuple([self.delete_box.children[0]] + [b.delete_btn for b in self.long_reinf])

    def make_long_reinf_tbl(self):
        # Define table layout properties.
        col_layout = wg.Layout(width='100%')
        
        # Create longitudinal steel columns.
        self.name_box = wg.VBox([wg.Label(value='name')], layout=col_layout)
        self.size_box = wg.VBox([wg.Label(value='size')], layout=col_layout)
        self.qty_box = wg.VBox([wg.Label(value='qty')], layout=col_layout)
        self.x_box = wg.VBox([wg.Label(value='x')], layout=col_layout)
        self.y_box = wg.VBox([wg.Label(value='y')], layout=col_layout)
        self.fy_box =  wg.VBox([wg.Label(value="fy")], layout=col_layout)
        self.Es_box =  wg.VBox([wg.Label(value="Es")], layout=col_layout)
        self.delete_box = wg.VBox([wg.Label(value='Delete')], layout=col_layout)
        
        # Create longitudinal table.
        long_reinf_tbl = wg.HBox([self.name_box, self.size_box, self.qty_box, 
                                  self.x_box, self.y_box, self.fy_box, self.Es_box, self.delete_box])
        return long_reinf_tbl

#%% Other Layouts
class AppliedLoadLayout(object):
    
    def __init__(self, parent, model):
        self.parent = parent
        self.model = model
        
        cell_layout = wg.Layout(width='auto')

        self.Mu = wg.FloatText(value=model.Mu, description=r'$M_{u}$ &nbsp <i>[lb*in]</i>', layout=cell_layout)
        self.Vu = wg.FloatText(value=model.Vu, description=r'$V_{u}$ &nbsp <i>[lb]</i>', layout=cell_layout)
        
    def calc(self):
        self.model.calc(Mu=self.Mu.value, Vu=self.Vu.value)

class ConcMaterialLayout(object):
    
    def __init__(self, parent, model):
        self.parent = parent
        self.model = model
        
        cell_layout = wg.Layout(width='auto')

        self.fc_prime = wg.FloatText(value=model.fc_prime, description=r'$f^{\prime}_{c}$ &nbsp <i>[psi]</i>', layout=cell_layout)
        self.strain_ult = wg.FloatText(value=model.strain_ult, description=r'$\epsilon_{u}$ &nbsp <i>[in/in]</i>', layout=cell_layout)
        self.Ec = wg.FloatText(value=model.Ec, description=r'$E_c$ &nbsp <i>[psi]</i>', layout=cell_layout)
        
    def calc(self):
        self.model.calc(fc_prime=self.fc_prime.value, strain_ult=self.strain_ult.value, Ec=self.Ec.value)
        
class VertexLayout(object):
    def __init__(self, parent):
        self.parent = parent
        
        cell_layout = wg.Layout(width='auto')

        self.x = wg.FloatText(value=0.0, layout=cell_layout)
        self.y = wg.FloatText(value=0.0, layout=cell_layout)
        
        self.delete_btn = wg.Button(description='Delete', layout=cell_layout)
        self.delete_btn.on_click(self.delete)

        self.insert_btn = wg.Button(description='Insert Vertex Below', layout=cell_layout)
        self.insert_btn.on_click(self.insert_vert_below)
        
    def delete(self, b):
        self.parent.delete_vert(self)
        
    def insert_vert_below(self, b):
        self.parent.insert_vert_below(self)


class XsectShapeLayout(object):
    
    def __init__(self, parent, model):
        self.parent = parent
        self.model = model
        self.verts = []
        
        cell_layout = wg.Layout(width='auto')

        self.shape_type = wg.Dropdown(description="Shape", options=['Rectangular', 'T-Beam', 'Custom'], value=model.shape_type, layout=cell_layout)
        self.b = wg.FloatText(value=model.b, description=r'$b$ &nbsp <i>[in]</i>', layout=cell_layout)
        self.h = wg.FloatText(value=model.h, description=r'$h$ &nbsp <i>[in]</i>', layout=cell_layout)
        self.bf = wg.FloatText(value=model.bf, description=r'$b_{f}$ &nbsp <i>[in]</i>', layout=cell_layout)
        self.hf = wg.FloatText(value=model.hf, description=r'$h_{f}$ &nbsp <i>[in]</i>', layout=cell_layout)
        
        self.b.observe(self.plot_update, names='value')
        self.h.observe(self.plot_update, names='value')
        self.bf.observe(self.plot_update, names='value')
        self.hf.observe(self.plot_update, names='value')
        self.shape_type.observe(self.change_shape, names='value')

        # Create new vertices button and table.
        self.new_vert_btn = wg.Button(description='New Vertex')
        self.new_vert_btn.on_click(self.new_vert)
        
    def initialize_verts_tbl(self):
        # Add the initial set of vertices.
        self.verts.clear()
        for v in self.model.verts_nc:
            layout = VertexLayout(self)
            self.verts.append(layout)
            layout.x.value = v[0]
            layout.y.value = v[1]

        # Build the initial table of vertices.
        col_layout = wg.Layout(width='100%')    # Define table layout properties.
        
        # Create columns for vertices.
        self.x_box = wg.VBox([wg.Label(value='x')] + [v.x for v in self.verts], layout=col_layout)
        self.y_box = wg.VBox([wg.Label(value='y')] + [v.y for v in self.verts], layout=col_layout)
        self.delete_box = wg.VBox([wg.Label(value='Delete')] + [v.delete_btn for v in self.verts], layout=col_layout)
        self.insert_box = wg.VBox([wg.Label(value='Insert Vertex Below')] + [v.insert_btn for v in self.verts], layout=col_layout)

        # Create longitudinal table.
        verts_tbl = wg.HBox([self.x_box, self.y_box, self.delete_box, self.insert_box])
        
        return verts_tbl
        
    def new_vert(self, b):
        layout = VertexLayout(self)
        self.verts.append(layout)
        self.rebuild_verts_tbl()
        self.parent.calc_plot(None)

    def rebuild_verts_tbl(self):
        # Distribute verts info into columns.
        self.x_box.children = tuple([self.x_box.children[0]] + [v.x for v in self.verts])
        self.y_box.children = tuple([ self.y_box.children[0]] + [v.y for v in self.verts])
        self.delete_box.children = tuple([self.delete_box.children[0]] + [v.delete_btn for v in self.verts])
        self.insert_box.children = tuple([self.insert_box.children[0]] + [v.insert_btn for v in self.verts])

    def delete_vert(self, vert):
        # Remove the vertex from the list of vertices in both the model and layout and rebuild the table.
        self.verts.remove(vert)
        self.rebuild_verts_tbl()
#         self.parent.calc_plot(None)
        
    def insert_vert_below(self, current_vert):
        new_vert = VertexLayout(self)
        self.verts.insert(self.verts.index(current_vert)+1, new_vert)
        self.rebuild_verts_tbl()
        
    def change_shape(self, change):
        if change['new'] == "Rectangular":
            self.parent.shape_input.children = (self.shape_type, self.b, self.h)
        elif change['new'] == "T-Beam":
            self.parent.shape_input.children = (self.shape_type, self.b, self.h, self.bf, self.hf)
        else: # Assume custom layout.
            verts_tbl = self.initialize_verts_tbl()
            self.parent.shape_input.children = (self.shape_type, self.new_vert_btn, verts_tbl)
        
    def plot_update(self, change):
        self.parent.calc_plot(None)
        
    def calc(self):
        verts = [[v.x.value, v.y.value] for v in self.verts]
        self.model.calc(shape_type=self.shape_type.value, verts=verts, b=self.b.value, h=self.h.value, bf=self.bf.value, hf=self.hf.value)
        
        
class LongReinfLayout(object):
    
    def __init__(self, parent, model):
        self.parent = parent
        self.model = model
        
        cell_layout = wg.Layout(width='auto')
        
        self.name = wg.IntText(value=model.name, layout=cell_layout, disabled=True)
        self.size = wg.Dropdown(options=['3', '4', '5', '6', '7', '8', '9', '10', '11', '14', '18'], 
                                value=model.size, layout=cell_layout)
        self.qty = wg.FloatText(value=model.qty, layout=cell_layout)
        self.x = wg.FloatText(value=model.x, layout=cell_layout)
        self.y = wg.FloatText(value=model.y, layout=cell_layout)
        self.fy = wg.FloatText(value=model.fy, layout=cell_layout)
        self.Es = wg.FloatText(value=model.Es, layout=cell_layout)
        
        self.name.observe(self.plot_update, names='value')
        self.size.observe(self.plot_update, names='value')
        self.qty.observe(self.plot_update, names='value')
        self.x.observe(self.plot_update, names='value')
        self.y.observe(self.plot_update, names='value')
        
        self.delete_btn = wg.Button(description='Delete', layout=cell_layout)
        self.delete_btn.on_click(self.delete)
    
    def plot_update(self, change):
        self.parent.calc_plot(None)
        
    def delete(self, b):
        self.parent.delete_long_reinf(self)
        
    def calc(self):
        self.model.calc(size=self.size.value, qty=self.qty.value, x=self.x.value, y=self.y.value, fy=self.fy.value, Es=self.Es.value)

#%% __main__

if __name__ == "__main__":
    from general_app_layout import AppLayout, StackedLayout
    from conc_bm import ConcreteBeam
    beam = ConcBeamLayout(ConcreteBeam())
    app = StackedLayout(beam.input_box, beam.plot_box, beam.output_box)
    app.main_view        
        
