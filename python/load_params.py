#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 14:32:01 2020

@author: timowenzel, Sarah Henn
"""

# import python packages
import python.network as nw
import gurobipy as gp
import pandas as pd
import numpy as np

def load_eco():
    params = {"eco": {"co2_el": 0.516,
                      "eeg_pv": 0.08,
                      "elec_price": 0.275,
                      "peak_price": 1.25
                      },
              "phy": {"rho_w": 1000,  # kg/m³,            density of water
                      "c_w": 4.18  # kJ/(kg*K),        heat capacity of water
                      }
              }

    return params


def load_net():
    # %% CREATE NETWORK BASED ON KERBER NETWORKS

    # empty dict to store net data
    net_data = {}
    # create network
    net = nw.create_kerbernet_bedburg(std_type="NAYY 4x185", p_load_mw=.002, q_load_mvar=0.,
                                      trafotype="0.63 MVA 10/0.4 kV", v_os=10.)
    # set trafo bounds due to technichal limits
    net_data["trafo_max"] = float(net.trafo.sn_mva * 1000.)

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Setting up the network 
    # specify grid nodes for whole grid and trafo; choose and allocate load, injection and battery nodes
    net_data["net_nodes"] = {}
    net_data["net_nodes"]["grid"] = net.bus.index.values
    net_data["net_nodes"]["trafo"] = net.trafo['lv_bus'].values
    net_data["net_nodes"]["load"] = net.load['bus'].values
    # nodeInj = net.load['bus'].to_numpy()
    net_data["net_nodes"]["bat"] = net.load['bus'].values
    net_data["gridnodes"] = list(net_data["net_nodes"]["grid"])

    #print('Bedburg', net_data["gridnodes"])
    # Bedburg grid
    section = {0: [52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69],
               1: [15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 73, 74, 75, 76, 77, 78, 79, 80, 81],
               2: [90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 7, 8, 9, 10, 11, 12, 13],
               3: [14, 0, 1, 2, 3, 4, 5, 6, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 106, 107, 108, 109, 110, 111, 112],
               4: [43, 44, 45, 46, 47, 48, 82, 83, 84, 85, 86, 87, 88, 89],
               5: [31, 32, 49, 50, 51, 70, 71, 72]}


    # create allocation pandapowernetwork and opt grid
    node_list = section[0] + section[1] + section[2] + section[3] + section[4] + section[5]

    net_data["grid_allo"] = pd.DataFrame(node_list, columns=["nodes"])
    #print('xd',net_data["grid_allo"] )
    net_data["grid_allo"]["gridnodes"] = net_data["net_nodes"]["load"]

    # extract existing lines 
    net_data["nodeLines"] = []
    for i in range(len(net.line['from_bus'])):
        net_data["nodeLines"].append((net.line['from_bus'][i], net.line['to_bus'][i]))
        net_data["nodeLines"] = gp.tuplelist(net_data["nodeLines"])



    # extract maximal current for lines
    # multiply with 400 V to get maximal power in kW          
    net_data["powerLine_max"] = {}
    net_data["res"] = {}
    net_data["length"] = {}
    net_data["resQua"] = {}
    for [n, m] in net_data["nodeLines"]:
        net_data["powerLine_max"][n, m] = (net.line['max_i_ka'][net_data["nodeLines"].index((n, m))]) * 400

        net_data["res"][n, m] = (net.line['r_ohm_per_km'][net_data["nodeLines"].index((n, m))]) # resistance for route [n, m]
        net_data["length"][n, m] = (net.line['length_km'][net_data["nodeLines"].index((n, m))]) # length for route [n, m]
        net_data["resQua"][n, m] = (net_data["res"][n, m]*net_data["length"][n, m])**2  # resistance^2 for route [n, m]

    # Berechnung der Widerstände zwischen je 2 Haushalten



    Resi = np.zeros((113, 113))
    for n in range(0, 113):
        for j in range(0, 113):
            if n in section[0] and j in section[0]:
                Resi[n][j] = abs(section[0].index(n) - section[0].index(j)) * 0.01 * 0.165 + 2 * 0.011 * 0.582
            elif n in section[0] and j in section[1]:
                Resi[n][j] = section[0].index(n) * 0.01 * 0.165 + 2 * 0.011 * 0.582 + section[1].index(j) * 0.014 * 0.165
            elif n in section[0] and j in section[2]:
                Resi[n][j] = section[0].index(n) * 0.01 * 0.165 + 2 * 0.011 * 0.582 + section[2].index(j) * 0.023 * 0.165
            elif n in section[0] and j in section[3]:
                Resi[n][j] = section[0].index(n) * 0.01 * 0.165 + 2 * 0.011 * 0.582 + section[3].index(j)* 0.025 * 0.165
            elif n in section[0] and j in section[4]:
                Resi[n][j] = section[0].index(n) * 0.01 * 0.165 + 2 * 0.011 * 0.582 + section[4].index(j) * 0.014 * 0.165
            elif n in section[0] and j in section[5]:
                Resi[n][j] = section[0].index(n) * 0.01 * 0.165 + 2 * 0.011 * 0.582 + section[5].index(j) * 0.05 * 0.165

            elif n in section[1] and j in section[1]:
                Resi[n][j] = abs(section[1].index(n) - section[1].index(j)) * 0.014 * 0.165 + 2 * 0.011 * 0.582

            elif n in section[1] and j in section[2]:
                Resi[n][j] = section[1].index(n) * 0.014 * 0.165 + 2 * 0.011 * 0.582 + section[2].index(j) * 0.023 * 0.165
            elif n in section[1] and j in section[3]:
                Resi[n][j] = section[1].index(n) * 0.014 * 0.165 + 2 * 0.011 * 0.582 + section[3].index(j) * 0.025 * 0.165
            elif n in section[1] and j in section[4]:
                Resi[n][j] = section[1].index(n) * 0.014 * 0.165 + 2 * 0.011 * 0.582 + section[4].index(j) * 0.014 * 0.165
            elif n in section[1] and j in section[5]:
                Resi[n][j] = section[1].index(n) * 0.014 * 0.165 + 2 * 0.011 * 0.582 + section[5].index(j) * 0.05 * 0.165


            elif n in section[2] and j in section[2]:
                Resi[n][j] = abs(section[2].index(n) - section[2].index(j)) * 0.023 * 0.165 + 2 * 0.011 * 0.582

            elif n in section[2] and j in section[3]:
                Resi[n][j] = section[2].index(n) * 0.023 * 0.165 + 2 * 0.011 * 0.582 + section[3].index(j) * 0.025 * 0.165
            elif n in section[2] and j in section[4]:
                Resi[n][j] = section[2].index(n)* 0.023 * 0.165 + 2 * 0.011 * 0.582 + section[4].index(j) * 0.014 * 0.165
            elif n in section[2] and j in section[5]:
                Resi[n][j] = section[2].index(n) * 0.023 * 0.165 + 2 * 0.011 * 0.582 + section[5].index(j) * 0.05 * 0.165


            elif n in section[3] and j in section[3]:
                Resi[n][j] = abs(section[3].index(n) - section[3].index(j)) * 0.025 * 0.165 + 2 * 0.011 * 0.582

            elif n in section[3] and j in section[4]:
                Resi[n][j] = section[3].index(n) * 0.025 * 0.165 + 2 * 0.011 * 0.582 + section[4].index(j) * 0.014 * 0.165
            elif n in section[3] and j in section[5]:
                Resi[n][j] = section[3].index(n) * 0.025 * 0.165 + 2 * 0.011 * 0.582 + section[5].index(j) * 0.05 * 0.165


            elif n in section[4] and j in section[4]:
                Resi[n][j] = abs(section[4].index(n) - section[4].index(j)) * 0.014 * 0.165 + 2 * 0.011 * 0.582

            elif n in section[4] and j in section[5]:
                Resi[n][j] = section[4].index(n) * 0.014 * 0.165 + 2 * 0.011 * 0.582 + section[5].index(j) * 0.05 * 0.165


            elif n in section[5] and j in section[5]:
                Resi[n][j] = abs(section[5].index(n) - section[5].index(j)) * 0.05 * 0.165 + 2 * 0.011 * 0.582

    for n in range(0, 113):
        for j in range(0, 113):

            if n in section[1] and j in section[0]:
                Resi[n][j] = Resi[j][n]
            elif n in section[2] and j in section[0]:
                Resi[n][j] = Resi[j][n]
            elif n in section[2] and j in section[1]:
                Resi[n][j] = Resi[j][n]
            elif n in section[3] and j in section[0]:
                Resi[n][j] = Resi[j][n]
            elif n in section[3] and j in section[1]:
                Resi[n][j] = Resi[j][n]
            elif n in section[3] and j in section[2]:
                Resi[n][j] = Resi[j][n]
            elif n in section[4] and j in section[0]:
                Resi[n][j] = Resi[j][n]
            elif n in section[4] and j in section[1]:
                Resi[n][j] = Resi[j][n]
            elif n in section[4] and j in section[2]:
                Resi[n][j] = Resi[j][n]
            elif n in section[4] and j in section[3]:
                Resi[n][j] = Resi[j][n]
            elif n in section[5] and j in section[0]:
                Resi[n][j] = Resi[j][n]
            elif n in section[5] and j in section[1]:
                Resi[n][j] = Resi[j][n]
            elif n in section[5] and j in section[2]:
                Resi[n][j] = Resi[j][n]
            elif n in section[5] and j in section[3]:
                Resi[n][j] = Resi[j][n]
            elif n in section[5] and j in section[4]:
                Resi[n][j] = Resi[j][n]


    net_data["om"] = Resi

    # Berechnung der Widerstände vom Netz zu Haushalten
    Resist = np.zeros((113))
    for n in range(0, 113):
        if n in section[0]:
            Resist[n] = section[0].index(n) * 0.01 * 0.165 + 0.011 * 0.582
        elif n in section[1]:
            Resist[n] = section[1].index(n) * 0.014 * 0.165 + 0.011 * 0.582
        elif n in section[2]:
            Resist[n] = section[2].index(n) * 0.023 * 0.165 + 0.011 * 0.582
        elif n in section[3]:
            Resist[n] = section[3].index(n) * 0.025 * 0.165 + 0.011 * 0.582
        elif n in section[4]:
            Resist[n] = section[4].index(n) * 0.014 * 0.165 + 0.011 * 0.582
        elif n in section[5]:
            Resist[n] = section[5].index(n) * 0.05 * 0.165 + 0.011 * 0.582

    net_data["omGrid"] = Resist



    return net_data
