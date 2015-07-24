#!/usr/bin/env python

"""Classes representing parameters for GeoClaw runs"""
"""change by Hui Tang"""

import os

import clawpack.clawutil.data as clawdata



class SedimentData(clawdata.ClawData):

    def __init__(self):
    
        super(SedimentData,self).__init__()
    
        # Sediments data
        self.add_attribute('Sediment_density',2650)
        self.add_attribute('Water_density',1000)
        self.add_attribute('Porosity',0.45)
        self.add_attribute('Maximum_sediment_concentration',0.2)
        self.add_attribute('Control_sediment_diffusion_coefficient',-200.0)
        self.add_attribute('Sediment_thickness_for_each_layer',0.05)
        self.add_attribute('Initial_active_sediment_layer',350e3)
        self.add_attribute('Number_Grain_size_classes',10)
        self.add_attribute('Number_sediment_layers',100)
        self.add_attribute('Number_ghost_cell',2)
        self.add_attribute('Water_depth_consider_sediment',0.01)
        self.add_attribute('Grainsize',[0.005,0.005])
        #Physics parameter
        self.add_attribute('gravity',9.81)
        self.add_attribute('kinematic_viscosity',1e-7)
        self.add_attribute('Water_temperature',20)
        self.add_attribute('Representative_wave_period',10.0)
        self.add_attribute('Threshold_water_depth',0.01)
        self.add_attribute('von_kaman_coefficient',0.41)
        self.add_attribute('mining_coeffcient',0.3)
        self.add_attribute('Coefficient_source_term',0.1)
        self.add_attribute('Minimum_adaptation_time',0.5)
        self.add_attribute('maximum_Shields_parameter',-1.0)
        self.add_attribute('Friction_coefficient_flow',0.003)
        self.add_attribute('horizontal_background_viscosity',0.1)
        self.add_attribute('turbulent_horizontal_viscosity',1.0)
        self.add_attribute('Water_depth_swtich',0.005)
        self.add_attribute('Critical_avalanching_slope_under_water',0.03)
        self.add_attribute('Critical_avalanching_slope_above_water',0.01)
        self.add_attribute('toler_for_sediment_flux_limitor',1e-6)
        #processes control paramter
        self.add_attribute('sourcesink_terms',0)
        self.add_attribute('have_avalanched','F')
        self.add_attribute('Include_avalanching',1)
        self.add_attribute('Switch_for_hard_structures',1)
        self.add_attribute('morphological_acceleration_factor',1)
        self.add_attribute('Coefficient_determining_scheme',1.0)
        self.add_attribute('short_wave',0.0)
        self.add_attribute('Start_time',0.0)
        self.add_attribute('Split_threshold',1.01)
        self.add_attribute('Merge_threshold',0.01)
        # algorithm parameters
        self.add_attribute('order_accuracy',1.0)
        self.add_attribute('fully_upwind',-1.0)
        self.add_attribute('flux_limiter_method','VanAlbada')
        self.add_attribute('sediment_concentration_method','soulsby_vanrijn')
        self.add_attribute('sediment_flux_method','SVL')
    
    def write(self,out_file='./sediment.data',data_source='setrun.py'):
    
        self.open_data_file(out_file,data_source)
        
        self.data_write('Sediment_density')
        self.data_write('Water_density')
        self.data_write('Porosity')
        self.data_write('Maximum_sediment_concentration')
        self.data_write('Control_sediment_diffusion_coefficient')
        self.data_write('Sediment_thickness_for_each_layer')
        self.data_write('Initial_active_sediment_layer')
        self.data_write('Number_Grain_size_classes')
        self.data_write('Number_sediment_layers')
        self.data_write('Number_ghost_cell')
        self.data_write('Water_depth_consider_sediment')
        #Physics parameter
        self.data_write('gravity')
        self.data_write('kinematic_viscosity')
        self.data_write('Water_temperature')
        self.data_write('Representative_wave_period')
        self.data_write('Threshold_water_depth')
        self.data_write('von_kaman_coefficient')
        self.data_write('mining_coeffcient')
        self.data_write('Coefficient_source_term')
        self.data_write('Minimum_adaptation_time')
        self.data_write('maximum_Shields_parameter')
        self.data_write('Friction_coefficient_flow')
        self.data_write('horizontal_background_viscosity')
        self.data_write('turbulent_horizontal_viscosity')
        self.data_write('Water_depth_swtich')
        self.data_write('Critical_avalanching_slope_under_water')
        self.data_write('Critical_avalanching_slope_above_water')
        self.data_write('toler_for_sediment_flux_limitor')
        #processes control paramter
        self.data_write('sourcesink_terms')
        self.data_write('have_avalanched')
        self.data_write('Include_avalanching')
        self.data_write('Switch_for_hard_structures')
        self.data_write('morphological_acceleration_factor')
        self.data_write('Coefficient_determining_scheme')
        self.data_write('short_wave')
        self.data_write('Start_time')
        self.data_write('Split_threshold')
        self.data_write('Merge_threshold')
        # algorithm parameters
        self.data_write('order_accuracy')
        self.data_write('fully_upwind')
        #Method control
        self.data_write('flux_limiter_method')
        self.data_write('sediment_concentration_method')
        self.data_write('sediment_flux_method')
        
        self.close_data_file()


























