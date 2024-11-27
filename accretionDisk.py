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

'''Functions'''

'''Plotting superfitial density vs radii'''

# radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3] #(Excluimos las celdas fantasma)
# rmed = 0.5 * (radii[1:] + radii[:-1]) #centramomos las celdas
# output_number = 1
# dens_out = np.fromfile(outputdir + f"gasdens{output_number}.dat").reshape(NY, NX)
# #plt.plot(r, dens_out, 'k--', label='Condición inicial') # Graficamos abs(vr) para usar lo
# for i in range(6,10,1):
#     dens = np.fromfile(f"{outputdir}/gasdens{i}.dat")
#     plt.plot(r, dens*unit_surf_density, label=i)
#     plt.xlim(0.4,2.5)
#     plt.xlabel('$r$ [UA]')
#     plt.ylabel('$\log_{10}\Sigma$ [g/cm$^3$]')
#     #plt.xscale('log')
#     #splt.yscale('log')
#     plt.legend(ncol=3, title='Outputs')
# plt.tight_layout()
# plt.savefig('sigma vs r.png', dpi=300)
# plt.show()

'''Plotting radial velocity vs radii'''

# radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3] #(Excluimos las celdas fantasma)
# rmed = 0.5 * (radii[1:] + radii[:-1]) #centramomos las celdas
# output_number = 50
# v_out = np.fromfile(outputdir + f"gasvy{output_number}.dat").reshape(NY, NX) #
# #plt.plot(r, abs(v_out), 'k--', label='Condición inicial') # Graficamos abs(vr) para usar lo
# for i in range(6,10,3):
#     vr = np.fromfile(f"{outputdir}/gasvy{i}.dat")
#     plt.plot(r, vr*unit_velocity, label=i)
#     plt.xlabel('$r$ [UA]')
#     plt.ylabel('$v_r$ [cm/s]')
#     #plt.xscale('log')
#     #plt.yscale('log')
#     plt.legend(ncol=3, title='Outputs')
# plt.tight_layout()
# plt.savefig('vr vs r.png', dpi=300)
# plt.show()

# Example function to fit: Polynomial (quadratic as an example)
def fit_function(r, a, b, c):
    return a * r**2 + b * r + c

# Load the radii and velocity data
radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3]  # Exclude ghost cells
rmed = 0.5 * (radii[1:] + radii[:-1])  # Center the cells

# Prepare the data for fitting
data = {}
for i in range(2,10):
    vr = np.fromfile(f"{outputdir}/gasvy{i}.dat")
    data[i] = (rmed, vr * unit_velocity)

# Fit the data and plot
plt.figure(figsize=(8, 6))
for i, (r, vr) in data.items():
    # Fit the data
    popt, pcov = curve_fit(fit_function, r, vr)
    vr_fit = fit_function(r, *popt)
    
    # Plot original data
    #plt.plot(r, vr, label=f'Output {i} (data)')
    # Plot fit
    plt.plot(r, vr_fit, linestyle='--', label=f'Output {i} (fit)')
    
    # Add text label at the end of each fit
    plt.text(r[-1], vr_fit[-1], f'{i}', fontsize=10, ha='left', va='center')

# Configure plot
plt.xlabel('$r$ [UA]')
plt.ylabel('$v_r$ [cm/s]')
plt.legend(title='Outputs')
plt.tight_layout()
plt.savefig('vr_vs_r_fit_with_labels.png', dpi=300)
plt.show()


'''Accreted Mass'''

# Plotting accreted mass vs radii

radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3] #(Excluimos las celdas fantasma)
rmed = 0.5 * (radii[1:] + radii[:-1]) #centramomos las celdas
output_number = 5
dens_out = np.fromfile(outputdir + f"gasdens{output_number}.dat").reshape(NY, NX)
v_out = np.fromfile(outputdir + f"gasvy{output_number}.dat").reshape(NY, NX)
#plt.plot(r, abs(v_out), 'k--', label='Condición inicial') # Graficamos abs(vr) para usar lo
for i in range(3,10,3):
    vr = np.fromfile(f"{outputdir}/gasvy{i}.dat")
    dens = np.fromfile(f"{outputdir}/gasdens{i}.dat")
    M = -2*np.pi*r*unit_length*dens*unit_surf_density*vr*unit_velocity/(2*10**(30))
    plt.plot(r, M, label=i)
    plt.xlim(0.4,2.5)
    plt.xlabel('$r$ [UA]')
    plt.ylabel('$\dot{M}$ [M$_\odot$/yr]')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend(ncol=3, title='Outputs')
plt.tight_layout()
plt.savefig('m vs r.png', dpi=300)
plt.show()

# Plotting accreted mass vs time

radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3] #(Excluimos las celdas fantasma)
rmed = 0.5 * (radii[1:] + radii[:-1]) #centramomos las celdas
output_number = 5
#dens_out = np.fromfile(outputdir + f"gasdens{output_number}.dat").reshape(NY, NX)
#v_out = np.fromfile(outputdir + f"gasvy{output_number}.dat").reshape(NY, NX)
#plt.plot(r, abs(v_out), 'k--', label='Condición inicial') # Graficamos abs(vr) para usar lo
for i in range(9,10):
    vr = np.fromfile(f"{outputdir}/gasvy{i}.dat")
    dens = np.fromfile(f"{outputdir}/gasdens{i}.dat")
    M = -2*np.pi*Rmin*unit_length*dens*unit_surf_density*vr*unit_velocity/(2*10**(30))
    plt.plot(np.array(np.linspace(1,128,128)), M, label=i)
    plt.xlim(1,2.5)
    plt.ylim(1.0e-12,1.0e-11)
    plt.xlabel('$t$ [yr]')
    plt.ylabel('$\dot{M}$ [M$_\odot$/yr]')
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(ncol=3, title='Outputs')
plt.tight_layout()
plt.savefig('m vs t.png', dpi=300)
plt.show()

# Plotting for Luminosity

alpha = 1.0e-2
radii = np.loadtxt(f"{outputdir}/domain_y.dat")[3:-3] #(Excluimos las celdas fantasma)
rmed = 0.5 * (radii[1:] + radii[:-1]) #centramomos las celdas
output_number = 5
dens_out = np.fromfile(outputdir + f"gasdens{output_number}.dat").reshape(NY, NX)
vphi_out = np.fromfile(outputdir + f"gasvx{output_number}.dat").reshape(NY, NX)
cs = np.fromfile(outputdir + f"gasenergy{output_number}.dat").reshape(NY, NX)
#plt.plot(r, abs(v_out), 'k--', label='Condición inicial') # Graficamos abs(vr) para usar lo
for i in range(9,10,1):
    vr = np.fromfile(f"{outputdir}/gasvy{i}.dat")
    vphi = np.fromfile(f"{outputdir}/gasvx{i}.dat")
    dens = np.fromfile(f"{outputdir}/gasdens{i}.dat")
    M = -2*np.pi*Rmin*unit_length*dens*unit_surf_density*vr*unit_velocity/(2*10**(30))
    D = 9/4*(alpha*(cs*unit_velocity)**2/vphi)*dens*unit_surf_density*vphi*unit_velocity*((G*M*unit_mass)/r*unit_length)**2
    L = G*M*unit_mass/r*unit_length
    plt.plot(r, L*unit_energy, label=i)
    #plt.xlim(0.4,2.5)
    plt.xlabel('$r$ [UA]')
    plt.ylabel('$L$ [erg]')
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(ncol=3, title='Output')
plt.tight_layout()
plt.savefig('L vs r.png', dpi=300)
plt.show()
