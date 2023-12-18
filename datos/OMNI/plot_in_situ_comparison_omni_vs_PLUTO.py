"""
Comparing in-situ data to model output
======================================

In this example, we'll see how to compare in-situ data taken by a spacecraft
to equivalent observables in the model along the 1D spacecraft trajectory.

This consists of three steps:
1. Load the spacecraft data
2. Generate the spacecraft trajectory at the in-situ data timestamps
3. Use this trajectory to take samples in the 3D model output

The first two steps are accomplished using ``heliopy.data``
(to get the in-situ data) and ``heliopy.spice`` (to get the trajectory). This
is then fed into `Variable.sample_at_coords` to get the model values, which
we then compare to the in-situ data.
"""
###############################################################################
# First, load the required modules.
import matplotlib.pyplot as plt
import numpy as np
import math

from heliopy.data import omni
import heliopy.data.spice as spicedata
import heliopy.spice as spice
from astropy import units as u

from datetime import datetime
from datetime import timedelta

from psipy.model import PLUTOOutput

from tools.MASweb import get_mas_path
from carrington_dates import get_time_interval
from psipy.model import MASOutput
from psipy.data import sample_data

import ai.cs

from ai import cs

###############################################################################
# Load a set of PLUTO output files.
pluto_path_rotating = '/home/juan/output/'
model_pluto_rotating = PLUTOOutput(pluto_path_rotating)
print(model_pluto_rotating.variables)

vr_model_pluto_rotating = model_pluto_rotating['vx1']
Br_model_pluto_rotating = model_pluto_rotating['Bx1']
rho_model_pluto_rotating = model_pluto_rotating['rho']



# set up the time interval and carrington rotation. 
case_study = "cr2165"
starttime, endtime, cr = get_time_interval(case_study)


###############################################################################
# Load omni data.
#
merged_data = omni.h0_mrg1hr(starttime, endtime)
merged_data.to_dataframe()

merged_vars = merged_data.columns

print(merged_vars)

bx_gse, by_gse, bz_gse = merged_data.quantity('BX_GSE'), merged_data.quantity('BY_GSE'), merged_data.quantity('BZ_GSE')
p, Np, Vp = merged_data.quantity('T'), merged_data.quantity('N'), merged_data.quantity('V')
time = merged_data.index
print(time)

#Convert from datatime.index to datetime - needed for ai.cs
times = time.to_pydatetime()
#print(times)

#merged_data.data
bx_values = np.array(merged_data.data['BX_GSE'])
by_values = np.array(merged_data.data['BY_GSE'])
bz_values = np.array(merged_data.data['BZ_GSE'])
t_values = np.array(merged_data.data['T'])
n_values = np.array(merged_data.data['N'])
v_values = np.array(merged_data.data['V'])

#need to convert to RTN
br, bt, bn = ai.cs.cxform(cs_from = 'GSE', cs_to = 'RTN', dt = times, x = bx_values, y = by_values, z = bz_values)

#convert RTN to HEEQ
br_HEEQ, bt_HEEQ, bn_HEEQ = ai.cs.cxform(cs_from = 'RTN', cs_to = 'HEEQ', dt = times, x = br, y = bt, z = bn)

#convert HEEQ to Stonyhurst heliographic coordinates 

br_SHC = np.sqrt(br_HEEQ**2 + bt_HEEQ**2 + bn_HEEQ **2)
bt_SHC = np.arctan(bn_HEEQ/np.sqrt(br_HEEQ**2 + bt_HEEQ**2))
bp_SHC = np.arctan(br_HEEQ/bt_HEEQ)

#convert RTN to RTP
br_RTP = br
bt_RTP = -bn
bp_RTP = bt

vr_RTP = v_values
rho_RTP = n_values

#do the average on the first 12 hours 
idx = 0
vr_mean = np.nanmean(v_values[(idx):(idx+12)])
rho_mean = np.nanmean(n_values[(idx):(idx+12)])
b_mean = np.nanmean(br[(idx):(idx+12)])

print('vr [km/s] mean/sd', vr_mean)
print('N [cm^-3] mean/sd', rho_mean)
print('Br [nT] mean/sd', b_mean)

# "12 hr avg", 6 behind + current + 6 ahead

vr_mean_12hrs = []
rho_mean_12hrs = []
br_mean_12hrs = []
avg_times = []

for i, obs in enumerate(v_values):
    if 5 < i <= len(v_values)-6:
         vr_mean_12hrs.append(sum(v_values[i-6:i+6])/12)
         rho_mean_12hrs.append(sum(n_values[i-6:i+6])/12)
         br_mean_12hrs.append(sum(br[i-6:i+6])/12)
         avg_times.append(time[i])

###############################################################################
# Generate the trajectory.
#
# We take the timestamps from the previously loaded data, and use `heliopy.spice`
# to generate the trajectory at these times.
earth_traj = spice.Trajectory('Earth')
earth_traj.generate_positions(times=merged_data.index, observing_body='Sun', frame='IAU_SUN')
earth_coords = earth_traj.coords
#print(earth_coords.radius)


fac_km_to_AU = 6.68459e-9

#Pete's solution
# quick fix is to hardwire the min/max of the model domain, but 
# really need to get the min/max in r of the model. 
if (np.min(earth_coords.radius.value) < 0.139600):
    print("Trajectory below 30 Rs....clipping trajectory...")
    r_omni   = (np.clip(earth_coords.radius.value,0.139600,0.99985))*u.km

###############################################################################
#Take a sample of the radial velocity for PLUTO 
 
npoints = 671
vr_model_sampled_pluto_rotating = vr_model_pluto_rotating.sample_at_coords(earth_coords.lon,
                                        earth_coords.lat,
                                        earth_coords.radius,
				         t=np.ones(npoints)*110)                         
###############################################################################
# We can now plot a comparison between the models and in-situ measurements.

fig, axs = plt.subplots(nrows=3, sharex=True, figsize=(20,10))#
vrlabels = ['300', '400', '500', '600', '700']
Nplabels = ['5', '10', '20', '30']
Brlabels = ['-7.5','-5.0', '-2.5', '0.0', '2.5', '5.0', '7.5', '8.0']
timelabels = ['2017-04-29', '2017-05-01', '2017-05-05', '2017-05-09', '2017-05-13', '2017-05-17', '2017-05-21', '2017-05-25']

axs[0].plot(avg_times,  vr_mean_12hrs, color="red", linewidth=3, label='OMNI 12-horas')
axs[0].plot(merged_data.index, vr_model_sampled_pluto_inertial, color = 'green', linewidth=3, label='modelo MHD')
axs[0].axvline(x=starttime, color='orange', linewidth=2)
axs[0].axvline(x=endtime, color='orange', linewidth=2)

axs[0].set_ylabel(r'$v_{r}$ (km/s)', fontsize='20')
#plt.ylim(200, 800)
axs[0].set_ylim(200,800)
axs[0].tick_params(axis='y', labelsize=20)
#axs[0].set_yticklabels(vrlabels, fontsize='20')
axs[0].legend(loc='best', shadow=True, ncol=2, bbox_to_anchor=(0.64, 0.53), fontsize=16)
axs[0].set_title("CR2165", fontsize='26')


###############################################################################
#Take a sample of the mass density for PLUTO output
npoints = 671

rho_model_sampled_pluto_rotating = rho_model_pluto_rotating.sample_at_coords(earth_coords.lon,
                                       earth_coords.lat,
                                       earth_coords.radius,
				        t=np.ones(npoints)*110)
                                  
###############################################################################
# We can now plot a comparison between the models and in-situ measurements.
#fig, ax = plt.subplots(figsize=(20, 8))#
axs[1].plot(avg_times, rho_mean_12hrs, color="red", linewidth=3, label='OMNI 12-horas')
axs[1].plot(merged_data.index, rho_model_sampled_pluto_inertial, color = 'green', linewidth=3, label='modelo MHD')
axs[1].axvline(x=starttime, color='orange', linewidth=3)
axs[1].axvline(x=endtime, color='orange', linewidth=3)

axs[1].set_ylabel(r'$N/cm^{3}$',  fontsize='20')
axs[1].tick_params(axis='y', labelsize=18)

###############################################################################
#Take a sample of the radial magnetic field for PLUTO outoput
npoints = 671
		                       
Br_model_sampled_pluto_rotating = Br_model_pluto_rotating.sample_at_coords(earth_coords.lon,
                                       earth_coords.lat,
                                       earth_coords.radius,
				         t=np.ones(npoints)*110)                             
###############################################################################
# We can now plot a comparison between the models and in-situ measurements.
b_fac_pluto = 0.0458505 #nT
axs[2].plot(avg_times, br_mean_12hrs, color="red", linewidth=3, label='OMNI 12-horas')
axs[2].plot(merged_data.index, Br_model_sampled_pluto_inertial*b_fac_pluto, color = 'green', linewidth=3, label='modelo MHD')
axs[2].axvline(x=starttime, color='orange', linewidth=3)
axs[2].axvline(x=endtime, color='orange', linewidth=3)


axs[2].set_ylabel(r'$B_{r}$ (nT)',  fontsize='20')
plt.ylim(-6, 6)
axs[2].tick_params(axis='y', labelsize=20)
axs[2].tick_params(axis='x', labelsize=20)
fig.autofmt_xdate()  
filename='/home/juan/Downloads/titulo/datos/OMNI/plots/cr2165_modelo_MHD_vs_OMNI.png'
plt.savefig(filename, dpi=150)    
plt.show()
