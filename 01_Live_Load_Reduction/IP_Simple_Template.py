from tkinter import Tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
import pickle
import ipywidgets as wg

class IP_Simple_Template(object):
    """
    Create New, Save, and Load.
    """
    
    def __init__(self, calcs, inputs):
        self.calcs = calcs
        self.inputs = inputs
        self.create_update_calc_title_btn()
        self.create_calc_select()
        self.create_new_btn()
        self.create_save_btn()
        self.create_load_btn()
        self.create_delete_btn()
        self.file_box = wg.HBox([self.new_btn, self.save_btn, self.load_btn])
        self.out_box = wg.Output()
        self.create_info_accordion()
        self.create_input_accordion()
        self.create_output_accordion()

    def create_info_accordion(self):
        # Create the title box.
        self.calc_title = wg.Text(value=self.calcs.current.title)
        self.calc_title_box = wg.HBox([wg.Label('Title (cannot be blank) = '), self.calc_title])
        
        # Create the notes box.
        self.calc_notes = wg.Textarea(value='', placeholder='Put your notes here.')
        self.calc_notes.observe(self.change_notes, names='value')
        self.calc_notes_box = wg.HBox([wg.Label('Notes = '), self.calc_notes])
        
        # Put labels next to the title and notes widgets.
        self.calc_title_hbox = wg.HBox([self.calc_title_box, self.update_calc_title_btn])
        self.info_vbox = wg.VBox([self.delete_btn, self.calc_title_hbox, self.calc_notes_box])
        
        # Create the accordion widget.
        self.info_accordion = wg.Accordion(children=[self.info_vbox])
        self.info_accordion.set_title(0, 'Info')

    def change_notes(self, change):
        self.calcs.current.notes = change['new']
        
    def create_input_accordion(self):
        input_boxes = []
        for k,v in self.inputs.items():
            input_boxes.append(wg.HBox([wg.Label(v['label']), v['widget']]))
        input_vbox = wg.VBox(input_boxes)
        self.input_accordion = wg.Accordion(children=[input_vbox])
        self.input_accordion.set_title(0, 'Input')

    def create_output_accordion(self):
        self.out_accordion = wg.Accordion(children=[self.out_box])
        self.out_accordion.set_title(0, 'Output')

    def clear_out(self):
        # Clear the calc log.
        self.calcs.current.calc_log = []
        # Clear the output widget.
        self.out_box.clear_output()
        
    def write_out(self):
        # Write the calc log to the output widget.
        with self.out_box:
            for line in self.calcs.current.calc_log:
                print(line)

    def create_update_calc_title_btn(self):
        # Update the calc info. The calc title should not be continuously observed for updated because it will 
        # also update the calc selector and the title in the calcs collection.
        self.update_calc_title_btn = wg.Button(description='Update Title')
        self.update_calc_title_btn.on_click(self.update_calc_title)
        
    def update_calc_title(self, b):
        # Change the title of the calc in the calcs collection. Do not allow a blank title and do not change the title if the 
        # user didn't change it.
        self.calcs.update_title(old=self.calcs.current.title, new=self.calc_title.value)
        # Update the calculation selector with the new calc.
        self.calc_select.options = self.calcs.get_titles()
    
    def create_calc_select(self):
        self.calc_select = wg.Select(options=self.calcs.get_titles())
        self.calc_select.observe(self.on_selection_change, names='value')
        
    def on_selection_change(self, change):
        # Make the selected calc the current calc.
        self.calcs.change_current(change['new'])
        self.update_input()
        
    def update_input(self):
        """
        The calculation input widgets need to be updated whenever a new calculation
        is create (which makes the calculation the current calc) or when a different
        calculation is selected in the selection box. If a new calc is created, it
        will automatically be made the current calc in the calcs collection object.
        If a different calculation is selected in the selection box, it will have
        been made the current calc before this method is called.        
        """
        # The calc title and notes are not included in the input list because they
        # are required for every calc object. Update these two parameters separately.
        self.calc_title.value = self.calcs.current.title
        self.calc_notes.value = self.calcs.current.notes
        for k,v in self.inputs.items():
            # Update the widgets with the loaded object values.
            v['widget'].value = getattr(self.calcs.current, k)

    def create_new_btn(self):
        self.new_btn = wg.Button(description='New')
        self.new_btn.on_click(self.new)
        
    def new(self, b):
        # Create a new calculation object.
        self.calcs.new_calc()
        # Update the calc selector.
        self.calc_select.options = self.calcs.get_titles()
        self.calc_select.value = self.calcs.current.title
        # Update the input widgets to the values from the new calc.
        self.update_input()
    
    def create_save_btn(self):
        self.save_btn = wg.Button(description='Save')
        self.save_btn.on_click(self.save)
        
    def save(self, b):
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Get the file path with a file browser.
        file_path = asksaveasfilename(defaultextension='.pkl')
        # Launch the GUI.
        #%gui tk
        
        if file_path != "":
            # Pickle the object.
            pickle_out = open(file_path, "wb")
            pickle.dump(self.calcs, pickle_out)
            pickle_out.close()
                
    def create_load_btn(self):
        self.load_btn = wg.Button(description='Open')
        self.load_btn.on_click(self.load)
        
    def load(self, b):
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Get the file path with a file browser.
        file_path = askopenfilename()
        # Launch the GUI.
        #%gui tk
    
        if file_path != "":
            # Open the file.
            pickle_in = open(file_path, "rb")
            # Unpickle and overwrite the previous object.
            self.calcs = pickle.load(pickle_in)
            self.calc_select.options = self.calcs.get_titles()
            self.calc_select.value = next(iter(self.calc_select.options))

    def create_delete_btn(self):
        self.delete_btn = wg.Button(description='Delete')
        self.delete_btn.on_click(self.delete_calc)

    def delete_calc(self, b):
        # Delete the current calculation 
        self.calcs.del_calc(self.calcs.current.title)
        # Update the calc selector.
        self.calc_select.options = self.calcs.get_titles()
        self.calc_select.value = self.calcs.current.title
        # Update the input widgets to the values from the next calc.
        self.update_input()

"""
The following is the equivalent code that would need to be posted in the IPython
notebook to get similar functionality to the template object above. This code will
not be updated synchronously with the object, so it is likely to be out-of-date.

# Create input widgets and tie each to a function to change the object parameters when the widget value changes.

def update_input(calc_title):
    calcs.change_current(calc_title)
    # Update the widgets with the loaded object values.
    title.value = calcs.current.title
    notes.value = calcs.current.notes
    TA.value = calcs.current.TA
    KLL_cat.value = calcs.current.KLL_cat
    num_floor.value = calcs.current.num_floor
    passenger_vehicle_garage.value = calcs.current.passenger_vehicle_garage
    assembly_area.value = calcs.current.assembly_area
calc_select = wg.Select(options=calcs.get_titles())

def calc_selection_change(change):
    # Make the selected calc the current calc.
    update_input(change['new'])
calc_select.observe(calc_selection_change, names='value')

# Update the calc info. The calc title should not be continuously observed for updated because it will 
# also update the calc selector and the title in the calcs collection.
update_calc_title_btn = wg.Button(description='Update Title')
def update_calc_title(b):
    # Change the title of the calc in the calcs collection. Do not allow a blank title and do not change the title if the 
    # user didn't change it.
    calcs.update_title(old=calcs.current.title, new=title.value)
    # Update the calculation selector with the new calc.        
    calc_select.options = calcs.get_titles()
update_calc_title_btn.on_click(update_calc_title)

# Create save and laod buttons.
from tkinter import ttk, Tk
from tkinter.filedialog import askopenfilename, asksaveasfilename
import pickle

new_btn = wg.Button(description='New')
def new(b):
    # Create a new calculation object and make it the current object.
    # https://www.geeksforgeeks.org/global-local-variables-python/
    calcs.new_calc()
    # Update the calc selector.
    calc_select.options = calcs.get_titles()
    calc_select.value = calcs.current.title
    # Update the input widgets to the values from the new calc.
    update_input(calcs.current.title)
new_btn.on_click(new)

save_btn = wg.Button(description='Save')
def save(b):
    # Create Tk root
    root = Tk()
    # Hide the main window
    root.withdraw()
    # Get the file path with a file browser.
    file_path = asksaveasfilename(defaultextension='.pkl')
    # Launch the GUI.
    %gui tk
    
    if file_path is not "":
        # Pickle the object.
        pickle_out = open(file_path, "wb")
        pickle.dump(calcs, pickle_out)
        pickle_out.close()
save_btn.on_click(save)

# Create save and laod buttons.
load_btn = wg.Button(description='Load')
def load(b):
    # Create Tk root
    root = Tk()
    # Hide the main window
    root.withdraw()
    # Get the file path with a file browser.
    file_path = askopenfilename()
    # Launch the GUI.
    %gui tk

    if file_path is not "":
        # Open the file.
        pickle_in = open(file_path, "rb")
        # Unpickle and overwrite the previous object.
        global calcs
        calcs = pickle.load(pickle_in)
        calc_select.options = calcs.get_titles()
        calc_select.value = next(iter(calc_select.options))
load_btn.on_click(load)

del_calc_btn = wg.Button(description='Delete')
def del_calc(b):
    # Delete the current calculation 
    calcs.del_calc(calcs.current.title)
    # Update the calc selector.
    calc_select.options = calcs.get_titles()
    calc_select.value = calcs.current.title
    # Update the input widgets to the values from the next calc.
    update_input(calcs.current.title)
del_calc_btn.on_click(del_calc)
"""