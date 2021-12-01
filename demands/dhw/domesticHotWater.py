#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

@author: Tobias
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import pycity_base.classes.demand.domestic_hot_water as DomesticHotWater

import pycity_base.classes.Timer
import pycity_base.classes.Weather
import pycity_base.classes.Environment
import pycity_base.classes.Prices
import pycity_base.classes.demand.Occupancy


def run_test(nr_occ=3,do_plot=True):
    timeDiscretization = 3600
    total_nb_timesteps = 365 * 24 * 60 * 60 / timeDiscretization
    timer = pycity_base.classes.Timer.Timer(timeDiscretization=timeDiscretization,
                                       timestepstotal=total_nb_timesteps)
    weather = pycity_base.classes.Weather.Weather(timer, useTRY=True)
    prices = pycity_base.classes.Prices.Prices()

    environment = pycity_base.classes.Environment.Environment(timer, weather,
                                                         prices)

    dhw_annex42 = DomesticHotWater.DomesticHotWater(environment,
                                                    tFlow=60,
                                                    thermal=True,
                                                    method=1,  # Annex 42
                                                    dailyConsumption=100,
                                                    supplyTemperature=25)

    results = dhw_annex42.get_power(currentValues=False)

    print('Results for Annex42 profile:')
    print()
    print("Thermal demand: " + str(results[0]))
    print("Required flow temperature: " + str(results[1]))
    print()
    
    #  Convert into energy values in kWh
    dhw_energy_curve = results[0] * timeDiscretization / (3600*1000)
    annual_energy_demand = np.sum(dhw_energy_curve)

    print('Annual dhw energy demand in kWh: ', annual_energy_demand)

    print('#----------------------------------------------------------------')
    #  #---------------------------------------------------------------------
    #  Compute active occupants for one year
    #  Max. occupancy is 5 people simultaneously
#    occupancy = np.random.geometric(p=0.8, size=6 * 24 * 365) - 1
#    occupancy = np.minimum(5, occupancy)
    occup_obj = pycity_base.classes.demand.Occupancy.Occupancy(environment,
                                                          number_occupants=nr_occ)
    occupancy = occup_obj.occupancy

    dhw_stochastical = DomesticHotWater.DomesticHotWater(environment,
                                                         tFlow=60,
                                                         thermal=True,
                                                         method=2,
                                                         supplyTemperature=20,
                                                         occupancy=occupancy)

    dhw_power_curve = dhw_stochastical.get_power(currentValues=False,
                                                 returnTemperature=False)
    #  Convert into energy values in kWh
    dhw_energy_curve = dhw_power_curve * timeDiscretization / (3600*1000)
    annual_energy_demand = np.sum(dhw_energy_curve)
    #  DHW volume flow curve in liters/hour
    volume_flow_curve = dhw_stochastical.water
    #  Recalc into water volume in liters
    water_volume_per_timestep = volume_flow_curve / 3600 * timeDiscretization
    # Average daily dhw consumption in liters
    av_daily_dhw_volume = np.sum(water_volume_per_timestep) / 365

    print('Results for stochastic DHW profile:\n')
    print('Max number of occupants:', max(occupancy))
    print('Annual dhw energy demand in kWh: ', annual_energy_demand)
    print('Average daily domestic hot water volume in liters:',
          av_daily_dhw_volume)

    if do_plot:
        ax1 = plt.subplot(2, 1, 1)
        plt.step(np.arange(8760) + 1, dhw_stochastical.loadcurve, linewidth=2)
        plt.ylabel("Heat demand in Watt")
        plt.xlim((0, 8760))

        plt.subplot(2, 1, 2, sharex=ax1)
        plt.step((np.arange(len(occupancy)) * 10 + 10) / 60, occupancy,
                 linewidth=2)
        plt.ylabel("Active occupants")
        offset = 0.2
        plt.ylim((-offset, max(occupancy) + offset))
        plt.yticks(list(range(int(max(occupancy) + 1))))

        plt.show()


    # create .txt for simulation input
    dhw_power_df = pd.DataFrame(dhw_power_curve)
    dhw_power_df.index = dhw_power_df.index * 3600
    dhw_power_df.to_csv('dhwDemand.txt', header=False)

    for i in range(101,133+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
    
    for i in range(201,219+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
        
    for i in range(301,319+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
        
    for i in range(401,402+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
        
    for i in range(501,509+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
            
    for i in range(601,608+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
        
    for i in range(701,716+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt', "w") as f:
            f.writelines(lines)
    
        f.close()
        
    for i in range(801,807+1):
        with open('dhwDemand.txt') as f:
            lines = f.readlines()
    
        lines.insert(0, "#1\n")
        lines.insert(1,'double dhwDemand(8760,2)\n')
    
        with open('dhwDemandClusterC' + str(i) +'.txt',"w") as f:
            f.writelines(lines)
    
        f.close()        


if __name__ == '__main__':
    #  Run program
    run_test(do_plot=True)
