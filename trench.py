#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:13:47 2019

@author: mc4117
"""

import thetis as th
import morphological_hydro_fns as morph
import numpy as np
import pandas as pd
import pylab as plt


def boundary_conditions_fn_trench(bathymetry_2d, flag, morfac=1, t_new=0, state='initial'):
    """
    Define boundary conditions for problem to be used in morphological section.

    Inputs:
    bathymetry_2d - underlying bathymetry
    flag - switch between morphological and hydrodynamic conditions (not used in this example)
    morfac - morphological scale factor used when calculating time dependent boundary conditions
    t_new - timestep model currently at used when calculating time dependent boundary conditions
    state - when 'initial' this is the initial boundary condition set; when 'update' these are the boundary
            conditions set during update forcings (ie. if fluc_bcs = True, this will be called)
    """
    left_bnd_id = 1
    right_bnd_id = 2
    left_string = ['flux']
    right_string = ['elev']

    # set boundary conditions
    swe_bnd = {}

    flux_constant = -0.22
    elev_constant2 = 0.397

    inflow_constant = [flux_constant]
    outflow_constant = [elev_constant2]
    return swe_bnd, left_bnd_id, right_bnd_id, inflow_constant, outflow_constant, left_string, right_string


# define mesh
lx = 16
ly = 1.1
nx = 32
ny = 5
mesh2d = th.RectangleMesh(nx, ny, lx, ly)

x, y = th.SpatialCoordinate(mesh2d)

# define function spaces
V = th.FunctionSpace(mesh2d, 'CG', 1)
P1_2d = th.FunctionSpace(mesh2d, 'DG', 1)
vectorP1_2d = th.VectorFunctionSpace(mesh2d, 'DG', 1)

# define underlying bathymetry
bathymetry_2d = th.Function(V, name='Bathymetry')
initialdepth = th.Constant(0.397)
depth_riv = th.Constant(initialdepth - 0.397)
depth_trench = th.Constant(depth_riv - 0.15)
depth_diff = depth_trench - depth_riv

trench = th.conditional(th.le(x, 5), depth_riv, th.conditional(th.le(x, 6.5), (1/1.5)*depth_diff*(x-6.5) + depth_trench,
                                                               th.conditional(th.le(x, 9.5), depth_trench, th.conditional(th.le(x, 11), -(1/1.5)*depth_diff*(x-11) + depth_riv,
                                                                                                                          depth_riv))))
bathymetry_2d.interpolate(-trench)

# if the hydrodynamics simulation has already been run with the required parameters
# then hydro can be set to False. This is because we do not need to run the initial
# hydrodyamics simulation again as the previous run can be used to initialise
# the full simulation
hydro = True

if hydro:
    # define initial elevation
    elev_init = th.Function(P1_2d).interpolate(th.Constant(0.4))
    uv_init = th.as_vector((0.51, 0.0))

    # set up solver object
    solver_obj, update_forcings_hydrodynamics = morph.hydrodynamics_only(boundary_conditions_fn_trench, mesh2d, bathymetry_2d, uv_init, elev_init, average_size=160 * (10**(-6)), dt=0.25, t_end=500)

    # run model
    solver_obj.iterate(update_forcings=update_forcings_hydrodynamics)

    uv, elev = solver_obj.fields.solution_2d.split()
    morph.export_final_state("hydrodynamics_trench", uv, elev)

# set up solver object for full simulation
solver_obj, update_forcings_tracer, diff_bathy, diff_bathy_file = morph.morphological(boundary_conditions_fn=boundary_conditions_fn_trench, morfac=10,
                                                                  morfac_transport=True, suspendedload=True, convectivevel=True, bedload=True,
                                                                  angle_correction=False, slope_eff=True, seccurrent=False, fluc_bcs=False, mesh2d=mesh2d,
                                                                  bathymetry_2d=bathymetry_2d, input_dir='hydrodynamics_trench', viscosity_hydro=10**(-6),
                                                                  ks=0.025, average_size=160 * (10**(-6)), dt=0.6, final_time=15*3600, beta_fn=1.3, surbeta2_fn=1/1.5, alpha_secc_fn=0.75)

# run model
solver_obj.iterate(update_forcings=update_forcings_tracer)

# extract final bedlevel results

xaxisthetis1 = []
bathymetrythetis1 = []

for i in np.linspace(0, 15.8, 80):
    xaxisthetis1.append(i)
    bathymetrythetis1.append(-solver_obj.fields.bathymetry_2d.at([i, 0.55]))

# plot thetis results against sisyphe and experiment
from matplotlib import colors as mcolors
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

plt.rc('font', family='Helvetica')

data = pd.read_excel('data/experimental_data.xlsx', sheet_name='recreatepaperrun', header=None)
diff_sisyphe = pd.read_excel('data/sisyphe_results.xlsx')
thetisdf = pd.read_csv('model_outputs/trench_bed_output.csv')

plt.scatter(data[0], data[1], label='Experimental Data')

plt.plot(thetisdf['x'], thetisdf['bathymetry'], '--', linewidth = 3, c = colors['darkgreen'], label = r'Thetis ($\Delta x = 0.2m$)')
plt.plot(xaxisthetis1, bathymetrythetis1, colors['mediumblue'], label = r'Thetis ($\Delta x = 0.5m$)')
plt.plot(diff_sisyphe['x'][diff_sisyphe['y'] == 0.55], -diff_sisyphe['Sisyphe'][diff_sisyphe['y'] == 0.55], color=colors['orange'], label='Sisyphe')

plt.xlabel('x (m)', fontsize=12)
plt.ylabel(r'$z_{b}$ (m)', fontsize=12)
plt.xlim([0, 16])
plt.ylim([-0.155, 0.01])
plt.legend(loc=3)
plt.show()

# record model output in csv file
df = pd.concat([pd.DataFrame(xaxisthetis1, columns = ['x']), pd.DataFrame(bathymetrythetis1, columns = ['bathymetry'])], axis=1)

df.to_csv('model_outputs/trench_bed_output_coarse.csv', index = False)
