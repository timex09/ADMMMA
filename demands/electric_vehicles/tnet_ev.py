import numpy as np
import scipy.io as spio
import pandas as pd
import matplotlib.pyplot as plt

# Root path for EV base profiles (single charger)
ev_profiles_base_path = (r"D:\GIT\bedburg\electric_vehicles\ev_profiles\EV_promoter_")
#r"D:\GIT\bedburg\electric_vehicles\ev_profiles\EV_promoter_"

# Time and periods input data
start = '2017-01-01'
periods = 365*96
freq = '15Min'

# Number of chargers in district by usage pattern and charger type
# Number per charger must not be above 550 --> Can be extended if needed
charger_struct = {'residential': {'3.3': (0, '3.3 kW @ 220V/16A'),
                                  '7.0': (20, '7.0 kW @ 220V/32A'),
                                  '43.0': (0, '43.0 kW @ 400V/63A')},
                  'commercial': {'3.3': (0, '3.3 kW @ 220V/16A'),
                                 '7.0': (0, '7.0 kW @ 220V/32A'),
                                 '43.0': (3, '43.0 kW @ 400V/63A')}}

# Considered vehicle type --> Currently only small, could be extended
vehicle_type = "small"
# usage_pattern = ["residential", "commercial"]

# What to draw from base charger profile --> Others can be e.g. EV occurrence
var_name = 'EV_uncontrolled_charging_profile'

# Create variable and result containers for loop
res_profile = np.expand_dims(np.array(np.arange(periods)), axis=1)
count_chargers = 0
charger_names = []

# Loop over different usage pattern
for _u, _usage in enumerate(charger_struct.keys()):
    # Loop over different charger types
    for _c, _charger in enumerate(charger_struct[_usage].items()):
        if _charger[1][0] > 0:
            _ev_load = spio.loadmat(
                ev_profiles_base_path + "{}_{}_{}_2017.mat".format(
                    _usage, vehicle_type, _charger[0])
            )[var_name][:, :_charger[1][0]]

            count_chargers += _charger[1][0]
            res_profile = np.hstack((res_profile, np.expand_dims(_ev_load.sum(axis=1), axis=1)))

            charger_names.append(_usage[:3] + ' ' + _charger[1][1])

# Put arrays into DataFrame structure
district_profile = pd.DataFrame(data=res_profile[:, 1:],
                                columns=charger_names,
                                index=pd.DatetimeIndex(start=start,
                                                       periods=periods,
                                                       freq=freq))
#district_profile['Total'] = district_profile.sum(axis=1)

# Create plot
ax1 = district_profile['2017-01-01':'2017-01-03'].plot()
ax1.set_ylabel('Charging power in kW')
plt.tight_layout()
plt.savefig('EVProfile.png', dpi=1200, transparent=True)
plt.savefig('EVProfile.pdf', dpi=1200)
plt.show()


# Save as csv
district_profile.to_csv('ev_profile_distrct')
