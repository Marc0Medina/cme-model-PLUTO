from datetime import datetime, date, timedelta 
import matplotlib.pyplot as plt
plt.style.use('seaborn')
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import numpy as np
import matplotlib.ticker as ticker
#def read_file(filename):
year_1 = [] # list time
doy = [] # ...
hora = [] # ...
LAT = []
LONG =[] 
nT = []
KMS =[] 
DP = []
TK = [] 


infile = open("omni2165.txt", "r") # reads file
for line in infile: # each line in txt file
    numbers = line.split() # removes the " "
    anyo = numbers[0]  # año
    di = numbers[1]   # dia
    hr = numbers[2]   # hr
    la = numbers[3]   #4 Lat
    lo = numbers[4]   #5 HGI Long. of the Earth
    nt = numbers[5]   #6 IMF BR, nT (RTN)
    ve = numbers[6]   #7 Speed, km/s
    de = numbers[7]   #8 Density, N/cm^3
    te = numbers[8]   #9 Temperature, K
    year_1.append(anyo) # Inserts 1st column(tt) into list(time) 
    doy.append(di) # ...
    hora.append(hr) # 
    LAT.append(float(la))
    LONG.append(float(lo))
    KMS.append(float(ve))
    nT.append(float(nt))
    TK.append(float(te))
    DP.append(float(de))

#Datos con fechas: dia y hora
concat_func = lambda x,y,z: str(x) + " " + str(y)+ " " + str(z)
concatenate_list = list(map(concat_func,year_1,doy, hora)) # list the map function
string_concatenate_list=' '.join([str(item) for item in concatenate_list])
dates = (datetime.strptime(string_concatenate_list, "%Y %j %H") for string_concatenate_list in concatenate_list)
date_strings = [datetime.strftime(d, '%Y/%m/%d-%H:00') for d in dates]



osdats = list(map(concat_func,date_strings,KMS,DP))
datos='\n'.join([str(item) for item in osdats])

file1 = open('Fluxes.txt', 'w')

file1.write(datos)






fig, axs = plt.subplots(4, 1, figsize=(6.4, 7), layout='constrained')
# common to all three:
for ax in axs:
    # Major ticks every half year, minor ticks every month,
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.grid(True)
    

# different formats:
ax = axs[0]
ax.plot(date_strings, KMS)
ax.set_title('Velocidad radial (km/s)', loc='left', fontsize='medium')
for label in ax.get_xticklabels(which='major'):
    label.set(rotation=8, horizontalalignment='right')

ax = axs[1]
ax.plot(date_strings, DP)
ax.set_title('Densidad (N/cm^3)', loc='left', fontsize='medium')
#ax.set_ylabel='Densidad (N/cm^3)'
for label in ax.get_xticklabels(which='major'):
    label.set(rotation=8, horizontalalignment='right')

ax = axs[2]
ax.plot(date_strings, nT)
ax.set_title('Campo magnético (nT)', loc='left', fontsize='medium')
#ax.ticklabel_format(axis='y',style='sci',scilimits=(7,0))
for label in ax.get_xticklabels(which='major'):
    label.set(rotation=8, horizontalalignment='right')

ax = axs[3]
ax.plot(date_strings, TK)
ax.set_title('Temperatura (K)', loc='left', fontsize='medium')
#ax.ticklabel_format(axis='y',style='sci',scilimits=(7,0))
for label in ax.get_xticklabels(which='major'):
    label.set(rotation=8, horizontalalignment='right')











plt.savefig('OMNI1.svg')
plt.show()

