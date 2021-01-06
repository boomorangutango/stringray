import math

class LL_Reduction_Factor(object):
    """
    Calculate reduced live loads.

    Legend:
    KLL = live load element factor
    TA = tributary area (ft^2)
    L = reduced design live load per ft^2 of area supported by the member

    Resources:
    (1) ASCE 7-10
    (2) IBC

    Notes:
    (1) Live loads may only be reduced if there are required minimum loads.
    (2) The live load reduction factor may be calculated by setting L = 1.
    (3) Special rules apply for one-way spans and are not considered.
    (4) Roof live loads from assemblies, gardens, or other special purposes may
        may be reduced as floor live loads.
    (5) Heavy loads are not considered. If the live load exceeds 100 PSF for any 
        portion of the area considered, that load should not be included in the 
        reduced loads calculated by this program.

    Legend:
    L =  reduced design live load per ft^2 (m^2) of area supported by the member
    Lo =  unreduced design live load per ft^2 (m^2) of area supported by the member (see Table 4-1)
    KLL = live load element factor (see Table 4-2)
    TA = tributary area in ft^2 (m^2). This is the total tributary area of the element from all floors.
    """
    
    def __init__(self, title='**Title Me**', notes='', TA=1.0, KLL_cat='Interior columns', 
                 num_floor=1, passenger_vehicle_garage=False, assembly_area=False):
        self.title = title
        self.notes = notes
        self.TA = TA
        self.KLL_cat = KLL_cat
        self.num_floor = num_floor
        self.passenger_vehicle_garage = passenger_vehicle_garage
        self.assembly_area = assembly_area
        self.calc_log = []
        
    def calc_KLL(self):
        #self.calc_log.append('Looking up KLL.')
        KLL_tbl = {'Interior columns': 4.0,
                   'Exterior columns without canti. slab': 4.0,
                   'Corner columns with canti. slabs': 3.0, 
                   'Edge beams without canti. slabs': 2.0, 
                   'Interior beams': 2.0, 
                   'All other members not identified': 1.0}
        self.KLL = KLL_tbl[self.KLL_cat]
        self.calc_log.append('KLL = {}'.format(self.KLL))
        
    def calc_LL_reduction_factor(self):
        self.calc_KLL()
        # Influence area.
        infl_area = self.KLL*self.TA
        #print(self.KLL*self.TA)
        #self.calc_log.append('Influence area = {}'.format(infl_area))
        
        # Minimum portions of the original load considering limitations on the type of load.
        if (infl_area < 400.0 or self.assembly_area == True):
            # Does not qualify for live load reductions.
            self.calc_log.append('Live load reduction cannot be used for influence area < 400 ft^2 or for assembly areas')
            self.reduction_factor_limit = 1
        elif (self.passenger_vehicle_garage == True and self.num_floor > 1):
            # Passenger vehicle garage live load reduction.
            self.calc_log.append('Max allowed live load reduction for members supporting more than one floor ' +
                                 'in passenger vehicle garages = 20%.')
            self.reduction_factor_limit = 0.8
        else:
            # Live loads may be reduced subject to normal limitations.
            if self.num_floor == 1:
                self.calc_log.append('Max allowed live load reduction = 50% for 1 floor.')
                self.reduction_factor_limit = 0.5
            else: # More than one floor considered.
                self.reduction_factor_limit = 0.4
                self.calc_log.append('Max allowed live load reduction = 60%.')
            
        # Calculated minimum portion of the original load.
        self.max_reduction_factor = min(1, 0.25 + 15.0/math.sqrt(infl_area))
        # The load reduction factor (the portion of the original load).
        self.reduction_factor = max(self.max_reduction_factor, self.reduction_factor_limit)
        self.calc_log.append('Reduction factor = {0:.2f}'.format(self.reduction_factor))
        
if __name__ == "__main__":
    l = LL_Reduction_Factor(title='', 
                            notes='', 
                            TA=3000, 
                            KLL_cat='Interior columns', 
                            num_floor=1, 
                            passenger_vehicle_garage=False, 
                            assembly_area=False)
    
    l.calc_LL_reduction_factor()