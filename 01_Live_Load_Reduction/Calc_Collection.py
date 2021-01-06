class Calc_Collection(object):
    """
    Collect live load reduction calculations and manage the collection.
    """
    def __init__(self, base_obj):
        self.base_obj = base_obj
        self.calcs = {}
        self.new_calc()
        
    def new_calc(self):
        self.current = self.base_obj()
        self.calcs[self.current.title] = self.current
    
    def del_calc(self, title):
        # If the deleted calc is self.current, re-establish current.
        # Don't allow the last calc to be deleted.
        if len(self.calcs) < 2:
            pass
        else:
            del self.calcs[title]
            # This fails for the last calc.
            if self.current.title == title:
                self.current = self.calcs[next(iter(self.calcs))]
        
    def get_calc(self, title):
        return self.calcs[title]
    
    def get_titles(self):
        return self.calcs.keys()
    
    def update_title(self, old, new):
        if new != "" and new != old:
            # Add the old calc to the calc collection under the new name.
            self.calcs[new] = self.calcs[old]
            # Update the title of the calc.
            self.current.title = new
            # Delete the calc under the old title.
            del self.calcs[old]
    
    def change_current(self, title):
        self.current = self.calcs[title]
