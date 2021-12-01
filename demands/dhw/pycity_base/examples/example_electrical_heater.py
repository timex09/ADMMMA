#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example of the electrical heater class.
"""

from __future__ import division

import numpy as np

import pycity_base.classes.supply.electrical_heater as eh

import pycity_base.classes.timer
import pycity_base.classes.weather
import pycity_base.classes.prices
import pycity_base.classes.environment


def run_example():
    # Create environment
    timer = pycity_base.classes.timer.Timer()
    weather = pycity_base.classes.weather.Weather(timer, use_TRY=True)
    prices = pycity_base.classes.prices.Prices()
    environment = pycity_base.classes.environment.Environment(timer, weather, prices)

    # Create Electrical Heater
    lower_activation_limit = 0.5
    q_nominal = 10000
    t_max = 90
    eta = 0.9
    heater = eh.ElectricalHeater(environment, q_nominal, eta, t_max, lower_activation_limit)

    # Print results
    print()
    print(("Type: " + heater.kind))
    print()
    print(("Maximum electricity input: " + str(heater.p_nominal)))
    print(("Maximum heat output: " + str(heater.q_nominal)))
    print(("Efficiency: " + str(heater.eta)))
    print(("Maximum flow temperature: " + str(heater.t_max)))
    print(("Lower activation limit: " + str(heater.lower_activation_limit)))

    print()
    print(("Nominals: " + str(heater.getNominalValues())))

    np.random.seed(0)
    result_p = np.random.rand(timer.timesteps_used_horizon) * heater.p_nominal
    result_q = np.random.rand(timer.timesteps_used_horizon) * q_nominal
    result_schedule = np.random.randint(2, size=timer.timesteps_used_horizon)
    heater.setResults(result_p, result_q, result_schedule)

    results = heater.getResults(True)
    print()
    print("Electricity input: " + str(results[0]))
    print()
    print("Heat output: " + str(results[1]))
    print()
    print("Schedule: " + str(results[2]))


if __name__ == '__main__':
    #  Run program
    run_example()
