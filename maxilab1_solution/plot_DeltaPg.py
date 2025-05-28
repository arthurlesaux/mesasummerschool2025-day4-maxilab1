import matplotlib.pyplot as plt
import numpy as np

from tomso import fgong, mesa

project_dir = "/Users/alesaux/MESA_summerschool_2025/"

fgong_filename1 = project_dir+"RC_nuclear_CBM/LOGS/profile41.data.fgong"
profile_fgong1 = fgong.load_fgong(fgong_filename1)

print(profile_fgong1)

profile_file1 = project_dir+"RC_nuclear_CBM/LOGS/profile40.data"
profile_data1 = mesa.load_profile(profile_file1)

print(profile_data1)

history_file1 = project_dir+"RC_nuclear_CBM/LOGS_OVf1_Ledoux/history.data"
history_data1 = mesa.load_history(history_file1)

history_file2 = project_dir+"RC_nuclear_CBM/LOGS_PCf1_Ledoux/history.data"
history_data2 = mesa.load_history(history_file2)

history_file3 = project_dir+"RC_nuclear_CBM/LOGS_maxOV_Ledoux/history.data"
history_data3 = mesa.load_history(history_file3)

history_file4 = project_dir+"RC_nuclear_CBM/LOGS_SC_Ledoux/history.data"
history_data4 = mesa.load_history(history_file4)

print(history_data1)

delta_Pg1 = history_data1["delta_Pg"]
Yc1 = history_data1["center_he4"]
star_age1 = history_data1["star_age"]

delta_Pg2 = history_data2["delta_Pg"]
Yc2 = history_data2["center_he4"]
star_age2 = history_data2["star_age"]

delta_Pg3 = history_data3["delta_Pg"]
Yc3 = history_data3["center_he4"]
star_age3 = history_data3["star_age"]

delta_Pg4 = history_data4["delta_Pg"]
Yc4 = history_data4["center_he4"]
star_age4 = history_data4["star_age"]

fig, (ax1, ax2) = plt.subplots(2,1, sharey= True)
ax1.set_title(r'Period spacing')
ax1.plot(star_age1/1e9,delta_Pg1, color='navy',label='OV')
ax1.plot(star_age2/1e9,delta_Pg2, color='red',label='PC')
ax1.plot(star_age3/1e9,delta_Pg3, color='green',label='MaxOV')
ax1.plot(star_age4/1e9,delta_Pg4, color='goldenrod',label='SC')
ax1.set_xlabel('Age (Gyr)')
ax1.set_ylabel(r'$\Delta \Pi$ (s)')
ax1.legend()

ax2.plot(Yc1,delta_Pg1, color='navy')
ax2.plot(Yc2,delta_Pg2, color='red')
ax2.plot(Yc3,delta_Pg3, color='green')
ax2.plot(Yc4,delta_Pg4, color='goldenrod')
ax2.set_xlabel(r'$Y_c$')
ax2.set_ylabel(r'$\Delta \Pi$ (s)')
ax2.invert_xaxis()

plt.subplots_adjust(hspace = 0.4)

plt.show()
