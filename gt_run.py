#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on  10 2021

@author: Tim Adamek
"""
# import python packages
import numpy as np

# import own functions
#import python.parameters as parameters
#import python.building_handler as handler

import python.handler_scenario_distribution as verteilt
import python.load_params as parameters
import python.load_demands as demand
import python.load_devices as devices
import pandas as pd
from pandas import DataFrame


def get_opti_inputs(options):

    # Load Params
    devs_temp = devices.load_devices()
    nodes_temp, building_params = demand.load_demands(options, devs_temp)
    nodes = devices.map_devices(nodes_temp, building_params, devs_temp)
    params = parameters.load_eco()
    net_data = parameters.load_net()


    return nodes, net_data, params


# %% GENERAL PARAMETERS
if __name__ == "__main__":



    options = {"time_hor": 24,
               "tweeks": 4,
               "buildings": 113,
               "pv_share": 0.5,
               "ev_status": True,
               "grid": True,  # True -> consider grid constraints, False -> dont
               "rho": 1,  # start value of rho
               "Allocation": True, # True -> Allocation due to Power, False -> Allocation due to Power^2*Resistance
               "admm_startval_dual": (parameters.load_eco()["eco"]["elec_price"] + parameters.load_eco()["eco"]["eeg_pv"]) / 2, #start value of the P2P price
               "admm_iteration_limit": 700,
               "admm_threshold": 0.01, # epsilon
               "sequence": True,  # True -> optimization simultaneously, False -> optimization after another
               "path_input": "C:/Users/Tim/Desktop/Ma EBC/Python/ADMM/input_data/bedburg/"}



    nodes, net_data, params = get_opti_inputs(options)

    # start optimization
    result_file = {}
    for week in range(options["tweeks"]):

        if week > 1:


            counter, status, costs, res, optiparams, rounds, result = verteilt.admm(nodes, options, week, net_data, params,result_file)

            if week == 0:
                df = DataFrame(
                    {'PtoGrid': result["20"], 'PfromGrid': result["21"], 'PfromP2P': result["22"],
                     'PtoP2P': result["23"],
                     'PreiseP2P': result["24"], 'PVleistung': result["25"]})
                df.to_excel('Woche1.xlsx', sheet_name='sheet1', index=False)

                kl = DataFrame({'Rechenzeit': [result["Berechnungszeit"]], 'Iterationen': [result["Iterationsschritt"]],
                                'Kosten': [result["Gesamtkosten"]], 'Haushalte': [result["Haushalte"]],
                                'MaxZeit': result["26"], 'NetzZeit': result["NetzZeit1"], 'Verluste': result["1Verl"]})
                kl.to_excel('Woche1Werte.xlsx', sheet_name='sheet1', index=False)

                xl = DataFrame({'ZFWERT': [result["30"]]})
                xl.to_excel('Woche1WerteZF.xlsx', sheet_name='sheet1', index=False)
            elif week == 1:
                df = DataFrame(
                    {'PtoGrid': result["20"], 'PfromGrid': result["21"], 'PfromP2P': result["22"],
                     'PtoP2P': result["23"],
                     'PreiseP2P': result["24"], 'PVleistung': result["25"]})
                df.to_excel('Woche2.xlsx', sheet_name='sheet1', index=False)
                kl = DataFrame({'Rechenzeit': [result["Berechnungszeit"]], 'Iterationen': [result["Iterationsschritt"]],
                                'Kosten': [result["Gesamtkosten"]], 'Haushalte': [result["Haushalte"]],
                                'MaxZeit': result["27"], 'NetzZeit': result["NetzZeit2"], 'Verluste': result["2Verl"]})
                kl.to_excel('Woche2Werte.xlsx', sheet_name='sheet1', index=False)
                xl = DataFrame({'ZFWERT': [result["31"]]})
                xl.to_excel('Woche2WerteZF.xlsx', sheet_name='sheet1', index=False)
            elif week == 2:
                df = DataFrame(
                    {'PtoGrid': result["20"], 'PfromGrid': result["21"], 'PfromP2P': result["22"],
                     'PtoP2P': result["23"],
                     'PreiseP2P': result["24"], 'PVleistung': result["25"]})
                df.to_excel('Woche3.xlsx', sheet_name='sheet1', index=False)
                kl = DataFrame({'Rechenzeit': [result["Berechnungszeit"]], 'Iterationen': [result["Iterationsschritt"]],
                                'Kosten': [result["Gesamtkosten"]], 'Haushalte': [result["Haushalte"]],
                                'MaxZeit': result["28"], 'NetzZeit': result["NetzZeit3"], 'Verluste': result["3Verl"]})
                kl.to_excel('Woche3Werte.xlsx', sheet_name='sheet1', index=False)
                xl = DataFrame({'ZFWERT': [result["32"]]})
                xl.to_excel('Woche3WerteZF.xlsx', sheet_name='sheet1', index=False)

            elif week==3:
                df = DataFrame(
                {'PtoGrid': result["20"], 'PfromGrid': result["21"], 'PfromP2P': result["22"], 'PtoP2P': result["23"],
                 'PreiseP2P': result["24"], 'PVleistung': result["25"]})
                df.to_excel('Woche4.xlsx', sheet_name='sheet1', index=False)
                kl = DataFrame({'Rechenzeit': [result["Berechnungszeit"]], 'Iterationen': [result["Iterationsschritt"]], 'Kosten': [result["Gesamtkosten"]], 'Haushalte': [result["Haushalte"]],'MaxZeit': result["29"],'NetzZeit': result["NetzZeit4"], 'Verluste': result["4Verl"]})
                kl.to_excel('Woche4Werte.xlsx', sheet_name='sheet1', index=False)

                xl = DataFrame({'ZFWERT': [result["33"]]})
                xl.to_excel('Woche4WerteZF.xlsx', sheet_name='sheet1', index=False)



