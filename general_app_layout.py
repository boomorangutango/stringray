import ipywidgets as wg

class AppLayout(object):
    
    def __init__(self, input_box, plot_box, output_box):
        
        # File layout.
        self.file_box = wg.Box()

        # Help layout.
        self.help_box = wg.Box()

        # Input area.
        narrow_input_btn = self.autosize_button('🡐')
        narrow_input_btn.on_click(self.on_narrow_input)
        widen_input_btn = self.autosize_button('🡒')
        widen_input_btn.on_click(self.on_widen_input)
        input_header = wg.HBox([wg.Label(value='Input'), narrow_input_btn, widen_input_btn])
        input_layout = wg.VBox([input_header, input_box], layout=wg.Layout(width='auto', grid_area='input'))

        # Plot area.
        narrow_plot = self.autosize_button('🡑')
        narrow_plot.on_click(self.on_narrow_plot)
        widen_plot = self.autosize_button('🡓')
        widen_plot.on_click(self.on_widen_plot)
        plot_header = wg.HBox([wg.Label(value='Plots'), narrow_plot, widen_plot])
        plot_layout = wg.VBox([plot_header, plot_box], layout=wg.Layout(width='auto', grid_area='plot'))

        # Output area.
        output_header = wg.HBox([wg.Label(value='Output')])
        output_layout = wg.VBox([output_header, output_box], layout=wg.Layout(width='auto', grid_area='output'))
        
        self.calc_layout = wg.GridBox(children=[input_layout, plot_layout, output_layout], 
                                      layout=wg.Layout(grid_template_rows='50% auto',
                                                       grid_template_columns='50% auto',
                                                       grid_template_areas='''
                                                       "input plot"
                                                       "input output"
                                                       '''))

        # Tabbed layout for overall app.
        self.app_layout = wg.Tab()
        self.app_layout.children = [self.calc_layout, self.help_box, self.file_box]
        self.app_layout.set_title(0, 'App')
        self.app_layout.set_title(1, 'Help')
        self.app_layout.set_title(2, 'File')

    def autosize_button(self, description):
        return wg.Button(description=description, layout=wg.Layout(height='auto', width='auto'))
    def on_narrow_input(self, b):
        current_width = self.calc_layout.layout.grid_template_columns.split()[0]
        width = max(float(current_width.split('%')[0]) - 10, 10)
        self.calc_layout.layout.grid_template_columns = f'{width}% auto'
    def on_widen_input(self, b):
        current_width = self.calc_layout.layout.grid_template_columns.split()[0]
        width = min(float(current_width.split('%')[0]) + 10, 90)
        self.calc_layout.layout.grid_template_columns = f'{width}% auto'
    def on_narrow_plot(self, b):
        current_width = self.calc_layout.layout.grid_template_rows.split()[0]
        width = max(float(current_width.split('%')[0]) - 10, 10)
        self.calc_layout.layout.grid_template_rows = f'{width}% auto'
    def on_widen_plot(self, b):
        current_width = self.calc_layout.layout.grid_template_rows.split()[0]
        width = min(float(current_width.split('%')[0]) + 10, 90)
        self.calc_layout.layout.grid_template_rows = f'{width}% auto'


class LinearLayout(object):
    
    def __init__(self, input_box, plot_box, output_box):
        
        # File layout.
        # self.file_box = wg.Box()

        # Help layout.
        # self.help_box = wg.Box()

        # Input area.
        # narrow_input_btn = self.autosize_button('🡐')
        # narrow_input_btn.on_click(self.on_narrow_input)
        # widen_input_btn = self.autosize_button('🡒')
        # widen_input_btn.on_click(self.on_widen_input)
        # input_header = wg.HBox([wg.Label(value='Input'), narrow_input_btn, widen_input_btn])
        # input_layout = wg.VBox([input_box], layout=wg.Layout(width='auto', grid_area='input'))

        # Plot area.
        # narrow_plot = self.autosize_button('🡑')
        # narrow_plot.on_click(self.on_narrow_plot)
        # widen_plot = self.autosize_button('🡓')
        # widen_plot.on_click(self.on_widen_plot)
        # plot_header = wg.HBox([wg.Label(value='Plots'), narrow_plot, widen_plot])
        # plot_layout = wg.VBox([plot_box], layout=wg.Layout(width='auto', grid_area='plot'))

        # Output area.
        # output_header = wg.HBox([wg.Label(value='Output')])
        # output_layout = wg.VBox([output_box], layout=wg.Layout(width='auto', grid_area='output'))
        
        # self.calc_layout = wg.GridBox(children=[input_layout, plot_layout, output_layout], 
        #                               layout=wg.Layout(grid_template_rows='50% auto',
        #                                                grid_template_columns='50% auto',
        #                                                grid_template_areas='''
        #                                                "input plot"
        #                                                "input output"
        #                                                '''))


        # Build the T-beam input layout.
        input_accordion = wg.Accordion(children=[input_box])
        input_accordion.set_title(0, 'Input')
        plot_accordion = wg.Accordion(children=[plot_box])
        plot_accordion.set_title(0, 'Plot')
        output_accordion = wg.Accordion(children=[output_box])
        output_accordion.set_title(0, 'Output')
        
        
        # Tabbed layout for overall app.
        self.app_layout = wg.VBox([input_accordion, plot_accordion, output_accordion])
        # self.app_layout.children = [, self.help_box, self.file_box]
        # self.app_layout.set_title(0, 'App')
        # self.app_layout.set_title(1, 'Help')
        # self.app_layout.set_title(2, 'File')

    # def autosize_button(self, description):
    #     return wg.Button(description=description, layout=wg.Layout(height='auto', width='auto'))
    # def on_narrow_input(self, b):
    #     current_width = self.calc_layout.layout.grid_template_columns.split()[0]
    #     width = max(float(current_width.split('%')[0]) - 10, 10)
    #     self.calc_layout.layout.grid_template_columns = f'{width}% auto'
    # def on_widen_input(self, b):
    #     current_width = self.calc_layout.layout.grid_template_columns.split()[0]
    #     width = min(float(current_width.split('%')[0]) + 10, 90)
    #     self.calc_layout.layout.grid_template_columns = f'{width}% auto'
    # def on_narrow_plot(self, b):
    #     current_width = self.calc_layout.layout.grid_template_rows.split()[0]
    #     width = max(float(current_width.split('%')[0]) - 10, 10)
    #     self.calc_layout.layout.grid_template_rows = f'{width}% auto'
    # def on_widen_plot(self, b):
    #     current_width = self.calc_layout.layout.grid_template_rows.split()[0]
    #     width = min(float(current_width.split('%')[0]) + 10, 90)
    #     self.calc_layout.layout.grid_template_rows = f'{width}% auto'
