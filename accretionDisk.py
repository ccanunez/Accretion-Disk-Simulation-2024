import numpy as np
import matplotlib.pyplot as plt

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

kb = 1.380650424e-16                   #erg/K    | erg = gr cm2 / s2
mp = 1.672623099e-24                   #gr
mu = 2.35                              #mol
gamma = 1.667              #
G = 6.67259e-8                         #cm3/gr/s2
Rgas = 8.314472e7                      #erg /K/mol
unit_mass = 1.9891e33                #gr
unit_length = 1.49598e13                #cm
        
#Unidades ---------------------------------------------------------#
unit_surf_density = (unit_mass/unit_length**2)     #gr/cm2
unit_volm_density = (unit_mass/unit_length**3)     #gr/cm3
unit_temperature  = ((G*mp*mu)/(kb))*(unit_mass/(unit_length))      #K
unit_time = np.sqrt( pow(unit_length,3.) / G / unit_mass) / 3.154e7 #yr
#unit_disktime = float((float(P["NINTERM"])*float(P["DT"])))*unit_time   #yr/snapshot    [Las unidades de tiempo del disco se toman en R=1]
unit_energy = unit_mass/(unit_time*3.154e7)**2/ unit_length          #erg
#unit_period   = unit_time*2*np.pi*np.sqrt(Rp**3)        #yr
#unit_orbits   = unit_disktime/unit_period                 #orbita/snapshot        [Fraccion de la orbita del planeta]
unit_velocity = unit_length/float(unit_time*3.154e7)                #cm/s

path = "/Users/ccanunez/Documents/USM/8vo Semestre/Discos de Acreción/Simulations/fargo3d/outputs/fargo/"
Rmin = 0.4 
Rmax = 2.5
#NX = 384 # resolución azimutal
#NY = 128 # resolución radial

# Cargar la grilla

domain_x = np.genfromtxt(path + "domain_x.dat")
domain_y = np.genfromtxt(path + "domain_y.dat")[3:-3] # extrae celdas fantasmas

NX = len(domain_x) - 1
NY = len(domain_y) - 1

output = 30
dens_out = np.fromfile(path + f"gasdens{output}.dat").reshape(NY,NX)

R = 0.5*(domain_y[1:]+domain_y[:-1])
phi = 0.5*(domain_x[1:]+domain_x[:-1])

P, R = np.meshgrid(phi, R)
X, Y = R*np.cos(P), R*np.sin(P)

plt.figure(figsize=(8,8))
plt.pcolormesh(X,Y, np.log10(dens_out))
plt.xlabel("$r$ [AU]")
plt.ylabel("$r$ [AU]")
plt.show()
