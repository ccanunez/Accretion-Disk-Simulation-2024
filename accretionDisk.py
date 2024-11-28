import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

# Plotting Style

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    'font.size': 12,         #Tamaño de fuente general
    'axes.titlesize': 12,    #Tamaño de fuente para títulos de ejes
    'axes.labelsize': 12,    #Tamaño de fuente para etiquetas de ejes
    'xtick.labelsize': 12,   #Tamaño de fuente para etiquetas del eje x
    'ytick.labelsize': 12,   #Tamaño de fuente para etiquetas del eje y
    'legend.fontsize': 11,   #Tamaño de fuente para la leyenda
    'figure.titlesize': 12   #Tamaño de fuente para el título de la figura
})
        
'''Conversión de Unidades'''

kb = 1.380650424e-16                   #erg/K    | erg = gr cm2 / s2
mp = 1.672623099e-24                   #gr
mu = 2.35                              #mol
gamma = 1.667                          #
G = 6.67259e-8                         #cm3/gr/s2
Rgas = 8.314472e7                      #erg /K/mol
unit_mass = 1.9891e33                  #gr
unit_length = 1.49598e13               #cm

unit_surf_density = (unit_mass/unit_length**2)     #gr/cm2
unit_volm_density = (unit_mass/unit_length**3)     #gr/cm3
unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length))      #K
unit_time = np.sqrt( pow(unit_length,3.) / G / unit_mass) / 3.154e7 #yr
#unit_disktime = float((float(P["NINTERM"])*float(P["DT"])))*unit_time   #yr/snapshot    [Las unidades de tiempo del disco se toman en R=1]
unit_energy = unit_mass/(unit_time*3.154e7)**2/ unit_length          #erg
#unit_period   = unit_time*2*np.pi*np.sqrt(Rp**3)        #yr
#unit_orbits   = unit_disktime/unit_period                 #orbita/snapshot        [Fraccion de la orbita del planeta]
unit_velocity = unit_length/float(unit_time*3.154e7)                #cm/s

'''Beginning of Code'''

outputdir = '/Users/ccanunez/Documents/USM/8vo Semestre/Discos de Acreción/Simulations/fargo3d/outputs/accretion/'

# Necessary Parameters

NX = 1
NY = 128
Rmin = 0.4
Rmax = 2.5
dr = (Rmax - Rmin)/NY
r = np.arange(Rmin, Rmax, dr)

'''Plotting superfitial density vs radii'''

outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs)))

plt.figure(figsize=(8, 6))
for i, output in enumerate(outputs):
    dens=np.fromfile(path+"gasdens"+str(output)+".dat").reshape(nr,nphi)
    plt.plot(r,(dens*unit_surf_density), label=f"{output} años", color=colors[i])
plt.xlabel("$r$ [AU]")
plt.ylabel("$\log \Sigma$ [g/cm$^2$]") 
#plt.title("Evolución de densidad superficial")
plt.legend(loc='best', title='Densidad')
plt.xlim(0.4, 2.4)
plt.savefig('dens vs r.png', dpi=300)
plt.show()

'''Plotting radial velocity vs radii'''

outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs)))

plt.figure(figsize=(8, 6))
for i, output in enumerate(outputs):
    gas_vy=np.fromfile(path+"gasvy"+str(output)+".dat").reshape(nr,nphi)
    plt.plot(r, gas_vy*unit_velocity, label=f"{output} años", color=colors[i])
plt.xlabel("$r$ [AU]")
plt.ylabel("$v_r$ [cm/s]")
#plt.title("Evolución de la velocidad radial")
plt.legend(loc='best', title='Velocidad Radial')
plt.xlim(0.4, 2.4)
plt.savefig('vr vs r.png', dpi=300)
plt.show()

'''Plotting azimutal velocity vs radii'''

outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs)))

plt.figure(figsize=(8, 6))
for i, output in enumerate(outputs):
    gas_vx=np.fromfile(path+"gasvx"+str(output)+".dat").reshape(nr,nphi)
    plt.plot(r, gas_vx*unit_velocity, label=f"{output} años", color=colors[i])
plt.xlabel("$r$ [AU]")
plt.ylabel("$v_\phi$ [cm/s]")
#plt.title("Evolución de la velocidad azimutal")
plt.legend(loc='best', title='Velocidad Azimutal')
plt.xlim(0.4, 2.4)
plt.savefig('vphi vs r.png', dpi=300)
plt.show()

'''Accreted Mass'''

# Plotting accreted mass vs radii

outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs)))  
plt.figure(figsize=(8, 6))
for j, output in enumerate(outputs):
    radio = []
    accretion_rate_r1 = []
    density = np.fromfile(path + f"gasdens{output}.dat").reshape(nr, nphi) * unit_surf_density
    vr = np.fromfile(path + f"gasvy{output}.dat").reshape(nr, nphi) * unit_velocity


    for i in range(len(r)):
        r_value = r[i]
        acc_rate_r = -np.sum(density[i, :] * vr[i, :] * 2 * np.pi * (r_value * unit_length))
        accretion_rate_r1.append(acc_rate_r)
        radio.append(r_value)


    radio = np.array(radio)
    accretion_rate_r1 = np.array(accretion_rate_r1)


    plt.plot(
        radio,
        accretion_rate_r1 / (1.898 * (10**30)),
        label=f"{output} años",
        color=colors[j]
    )

# Configurar el gráfico
plt.xlabel("$r$ [AU]")
plt.ylabel("$\dot{M}$ [$M_{\odot}$/año]")
#plt.title("Tasa de acreción vs radio en diferentes tiempos")
plt.xlim(0.4, 2.5)  # Límite del radio
plt.legend(loc="upper right", fontsize="small", title='Tasa de Acreción')
plt.savefig('m vs r.png', dpi=300)
plt.show()


# Plotting accreted mass vs time

r_med=0.4
outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs))) 
time=[]
acc_rate=[]

for i in range(0,outputs_totales):
    density=np.fromfile(path+f"gasdens{i}.dat").reshape(nr, nphi)*unit_surf_density
    vr=np.fromfile(path+f"gasvy{i}.dat").reshape(nr, nphi)*unit_velocity
    acc_rate_i=-(density[0]*vr[0]*2*np.pi*(r_med*unit_length))
    time.append(i)
    acc_rate.append(acc_rate_i)
time = np.array(time)
acc_rate = np.array(acc_rate)

plt.figure(figsize=(8, 6))
plt.plot(
    time,
    acc_rate / (1.898 * (10**30)),
    label=f'Tasa de acreción a $r={r_med}$ [UA]',
    color='b'
)
plt.xlabel('$t$ [año]')
plt.ylabel('$\dot{M}$ [$M_{\odot}$/año]')
#plt.title(f'Tasa de acreción en el tiempo (outputs seleccionados) a r={r_med} [UA]')
plt.xlim(0, outputs_totales)
plt.legend(loc='lower right')
plt.savefig('m vs t.png', dpi=300)
plt.show()

'''Plotting luminostity vs time'''

# Plotting for Luminosity for Stellar Disk

time = []
luminosity_over_time = []

# Parámetros
alpha = 1e-2  # viscosidad alfa
radii = r * unit_length  # Radios físicos en cm
H = 0.05 *radii # H/r constante
dr = (r[1]-r[0]) * unit_length  # Diferencia radial constante en cm

for i in range(outputs_totales):
    # Leer densidad superficial y velocidad angular
    density = np.fromfile(path + f"gasdens{i}.dat").reshape(nr) * unit_surf_density
    v_theta = np.fromfile(path + f"gasvx{i}.dat").reshape(nr) * unit_velocity
    cs = np.fromfile(path + f"gasenergy{i}.dat").reshape(nr) * unit_velocity
    # Calcular frecuencia angular Omega
    Omega = v_theta / radii

    # Calcular viscosidad cinemática ν = H * c_s * α
    nu = alpha * H * cs  # Aquí v_theta ~ c_s, suposición en discos delgados

    # Calcular luminosidad L = ∫ (9/4) * ν * Σ * Ω^2 * (2πr) dr
    luminosity = np.sum((9/4) * nu * density * Omega**2 * 2 * np.pi * radii * dr)
    luminosity_over_time.append(luminosity)

    # Almacenar tiempo
    time.append(i)

# Convertir listas a arrays de numpy
time = np.array(time)
luminosity_over_time = np.array(luminosity_over_time)

# Graficar luminosidad vs. tiempo
plt.figure(figsize=(10, 6))
plt.plot(time, luminosity_over_time, color='blue', label='Luminosidad del Disco Protoestelar')
plt.xlabel('$t$ [años]')
plt.ylabel('$L$ [erg/s]')
#plt.title('$L_{disk} vs t$')
plt.legend(loc='upper right')
plt.xlim(0, outputs_totales-1)
plt.savefig('Lproto vs t.png', dpi=300)
plt.show()

# Plotting Luminostity for AGN

r_med=0.4
outputs = range(0, outputs_totales+1, int(outputs_totales/10)) 
colors = plt.cm.viridis(np.linspace(0, 1, len(outputs))) 
time=[]
acc_rate=[]
luminosidad=[]
for i in range(0,outputs_totales):
    density=np.fromfile(path+f"gasdens{i}.dat").reshape(nr, nphi)*unit_surf_density
    vr=np.fromfile(path+f"gasvy{i}.dat").reshape(nr, nphi)*unit_velocity
    acc_rate_i=-(density[0]*vr[0]*2*np.pi*(r_med*unit_length))
    lum=acc_rate_i*unit_mass*G/(r_med*unit_length)
    time.append(i)
    acc_rate.append(acc_rate_i)
    luminosidad.append(lum)
time = np.array(time)
acc_rate = np.array(acc_rate)
luminosidad = np.array(luminosidad)
plt.figure(figsize=(8, 6))
plt.plot(
    time,
    luminosidad,
    label=f'Luminosidad del Disco AGN',
    color='b'
)
plt.xlabel('$t$ [años]')
plt.ylabel('$L$ [erg/s]')
#plt.title('$L_{acc} vs t$')
plt.xlim(0, outputs_totales)
plt.legend(loc='lower right')
plt.savefig('Lagn vs t.png', dpi=300)
plt.show()
