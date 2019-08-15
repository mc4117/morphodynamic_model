#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:11:48 2019

@author: mc4117
"""

import thetis as th
import morphological_hydro_fns as morph

import pandas as pd
import numpy as np
import pylab as plt

import time

def boundary_conditions_fn(bathymetry_2d, flag, morfac = 1, t_new = 0, state = 'initial'):
    """
    Define boundary conditions for problem to be used in morphological section.
    
    Inputs:
    morfac - morphological scale factor used when calculating time dependent boundary conditions
    t_new - timestep model currently at used when calculating time dependent boundary conditions
    state - when 'initial' this is the initial boundary condition set; when 'update' these are the boundary
            conditions set during update forcings (ie. if fluc_bcs = True, this will be called)
    """
    left_bnd_id = 1
    right_bnd_id = 2

    
    
    # set boundary conditions

    gradient_flux = (-0.053 + 0.02)/6000
    gradient_flux2 = (-0.02+0.053)/(18000-6000)
    gradient_elev = (10.04414- 9.9955)/6000
    gradient_elev2 = (9.9955-10.04414)/(18000-6000)
    elev_init_const = (-max(bathymetry_2d.dat.data[:]) + 0.05436)
    swe_bnd = {}
    swe_bnd[3] = {'un': th.Constant(0.0)}
    
    if state == 'initial':
        # initial boundary conditions
        if flag == 'hydro':
            left_string = ['flux']
            right_string = ['elev', 'flux']
            flux_constant = -0.02
            elev_constant2 = elev_init_const
        
            inflow_constant = [flux_constant]
            outflow_constant = [elev_constant2, -flux_constant]
            return swe_bnd, left_bnd_id, right_bnd_id, inflow_constant, outflow_constant, left_string, right_string
        elif flag == 'morpho':
            left_string = ['flux']
            right_string = ['elev']
            flux_constant = -0.02
            elev_constant2 = elev_init_const
        
            inflow_constant = [flux_constant]
            outflow_constant = [elev_constant2]
            return swe_bnd, left_bnd_id, right_bnd_id, inflow_constant, outflow_constant, left_string, right_string

    elif state == 'update':
        # update boundary condtions
        if t_new*morfac <= 6000:
            elev_constant2 = gradient_elev*t_new*morfac + elev_init_const
            flux_constant = (gradient_flux*t_new*morfac) - 0.02
        else:
            flux_constant = (gradient_flux2*(t_new*morfac-6000)) - 0.053
            elev_constant2 = gradient_elev2*(t_new*morfac-18000) + elev_init_const      
        
        inflow_constant = [flux_constant]
        outflow_constant = [elev_constant2]

        return inflow_constant, outflow_constant

# define mesh
mesh2d = th.Mesh("meander_test_3112_test.msh")
x,y = th.SpatialCoordinate(mesh2d)

# define function spaces
V = th.FunctionSpace(mesh2d, 'CG', 1)
P1_2d = th.FunctionSpace(mesh2d, 'DG', 1)
vectorP1_2d = th.VectorFunctionSpace(mesh2d, 'DG', 1)

# define underlying bathymetry

bathymetry_2d = th.Function(V, name='Bathymetry')

gradient = th.Constant(0.0035)

L_function=th.Function(V).interpolate(th.conditional(x > 5, th.pi*4*((th.pi/2)-th.acos((x-5)/(th.sqrt((x-5)**2+(y-2.5)**2))))/th.pi, th.pi*4*((th.pi/2)-th.acos((-x+5)/(th.sqrt((x-5)**2+(y-2.5)**2))))/th.pi))
bathymetry_2d1 = th.Function(V).interpolate(th.conditional(y > 2.5, th.conditional(x < 5, (L_function*gradient) + 9.97072, -(L_function*gradient) + 9.97072), 9.97072))

init = max(bathymetry_2d1.dat.data[:])
final = min(bathymetry_2d1.dat.data[:])

bathymetry_2d2 = th.Function(V).interpolate(th.conditional(x <= 5, th.conditional(y<=2.5, -9.97072 + gradient*abs(y - 2.5) + init, 0),th.conditional(y<=2.5, -9.97072 -gradient*abs(y - 2.5) + final, 0)))
bathymetry_2d = th.Function(V).interpolate(-bathymetry_2d1 - bathymetry_2d2)

input_bathymetry_2d = th.Function(V).interpolate(bathymetry_2d)
initial_bathymetry_2d = th.Function(V).interpolate(bathymetry_2d)

# define initial elevation
elev_init = th.Function(P1_2d).interpolate(0.0544 - bathymetry_2d)
# define initial velocity
uv_init = th.Function(vectorP1_2d).interpolate(th.as_vector((0.001,0.001)))

# set up solver object
solver_obj, update_forcings_hydrodynamics = morph.hydrodynamics_only(boundary_conditions_fn, mesh2d, bathymetry_2d, uv_init, elev_init, ks = 0.003, average_size = 10**(-3), dt = 1, t_end=200, viscosity = 5*10**(-2))

# run model
solver_obj.iterate(update_forcings = update_forcings_hydrodynamics)


uv, elev = solver_obj.fields.solution_2d.split()
morph.export_final_state("hydrodynamics_meander_initial_ext_grd_linux", uv, elev)

# set up solver object

solver_obj, update_forcings_tracer, diff_bathy, diff_bathy_file = morph.morphological(boundary_conditions_fn = boundary_conditions_fn, morfac = 10, morfac_transport = True, suspendedload = False, convectivevel = False,\
                    bedload = True, angle_correction = True, slope_eff = True, seccurrent = True, fluc_bcs = True,\
                 mesh2d = mesh2d, bathymetry_2d = input_bathymetry_2d, input_dir = 'hydrodynamics_meander_initial_ext_grd_linux', viscosity_hydro = 5*10**(-2), ks = 0.003, average_size = 10**(-3), dt = 2, final_time = 18000,\
                 beta_fn = 1.3, surbeta2_fn = 1/1.5, alpha_secc_fn = 0.75)

# run model
t1 = time.time()
solver_obj.iterate(update_forcings = update_forcings_tracer)
t2 = time.time()

print(t2-t1)

# find difference between initial and final bedlevel
diff_bathy.interpolate(-solver_obj.fields.bathymetry_2d + initial_bathymetry_2d)
diff_bathy_file.write(diff_bathy)


scaled_evolution = th.Function(V).interpolate(diff_bathy/0.0544)

from matplotlib import colors as mcolors

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# plot 90 degree cross-section comparing thetis with sisyphe and experimental data

y_array = np.linspace(6, 7, 100)

scaled_evolution_list = []

for i in y_array[0:len(y_array)-1]:
    scaled_evolution_list.append(scaled_evolution.at([4.5, i]))
    
plt.plot(y_array[0:len(y_array)-1]-6, scaled_evolution_list, color = colors['mediumblue'], label = 'Thetis')


paper_90 = pd.read_csv('data/paper_data_90.csv', header = 1)

plt.plot(paper_90['Sim Distance from inner bank'], paper_90['Sim Evolution'], color = colors['orange'], label = 'Sisyphe')

plt.plot(paper_90['Exp Distance from inner bank'], paper_90['Exp Evolution'], '--', linewidth =2, label = 'Experimental Data')

plt.xlabel('Distance from inner bank (m)')
plt.ylabel('Scaled Bedlevel Evolution')
plt.legend()

plt.show()

# plot 180 degree cross-section comparing thetis with sisyphe and experimental data

scaled_evolution_list_180 = []

x_array = np.linspace(8, 9, 100)

for i in x_array[0:len(x_array)-1]:
    scaled_evolution_list_180.append(scaled_evolution.at([i, 2.5]))
    
plt.plot(x_array[0:len(x_array)-1]-8, scaled_evolution_list_180, color = colors['mediumblue'], label = 'Thetis')

paper_180 = pd.read_csv('data/paper_data_180.csv', header = 0)
plt.plot(paper_180['Sim Distance from inner bank'], paper_180['Sim Evolution'], color = colors['orange'], label = 'Sisyphe')
plt.plot(paper_180['Exp Distance from inner bank'], paper_180['Exp Evolution'], '--', linewidth =2, label = 'Experimental Data')

plt.xlabel('Distance from inner bank (m)')
plt.ylabel('Scaled Bedlevel Evolution')
plt.legend()
plt.show()
