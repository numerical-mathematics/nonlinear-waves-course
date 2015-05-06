#!/usr/bin/env python
# encoding: utf-8
r"""
Force-field challenge
=====================

Solve the one-dimensional Euler equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u)_x & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x & = 0 \\
    E_t + (u (E + p) )_x & = 0.

The fluid is an ideal gas, with pressure given by :math:`p=\rho (\gamma-1)e` where
e is internal energy.

A blast wave is approaching!  Can you stop it by heating the air between you (at x=9)
and it (at 1<x<2) ?

To run the code in IPython:

    >>> run forcefield.py

This will print out the maximum pressure reached at x=9.
You can plot pressure versus time via

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t,p)

You can also plot the solution itself via

    >>> claw.plot()
"""
import numpy as np
from clawpack import riemann
from clawpack.riemann.euler_with_efix_1D_constants import *
gamma = 1.4 # Ratio of specific heats

def gauge_pressure(q,aux):
    internal_energy = q[energy] - 0.5*q[momentum]**2/q[density]
    pressure = q[density]*(gamma-1)*internal_energy
    return [pressure]

def setplot(plotdata):
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = 'Density'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = density
    plotitem.kwargs = {'linewidth':3}
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'Energy'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = energy
    plotitem.kwargs = {'linewidth':3}
    
    return plotdata


from clawpack import pyclaw

rs = riemann.euler_with_efix_1D

solver = pyclaw.ClawSolver1D(rs)
solver.kernel_language = 'Fortran'

solver.bc_lower[0]=pyclaw.BC.extrap
solver.bc_upper[0]=pyclaw.BC.extrap

mx = 800;
x = pyclaw.Dimension('x',0.0,10.0,mx)
domain = pyclaw.Domain([x])
state = pyclaw.State(domain,num_eqn)

state.problem_data['gamma'] = gamma
xc = state.grid.x.centers

pressure = 1. + (xc>1)*(xc<2)* 1000.  # Modify this line
# You can change anything except the velocity in the interval (3,8)

state.q[density ,:] = 1.
state.q[momentum,:] = 0.
state.q[energy  ,:] = pressure / (state.q[density,:]*(gamma - 1.))

grid = state.grid
grid.add_gauges([ [9.0] ])
#grid.setup_gauge_files()
solver.compute_gauge_values = gauge_pressure
state.keep_gauges = True

claw = pyclaw.Controller()
claw.tfinal = 2.
claw.solution = pyclaw.Solution(state,domain)
claw.solver = solver
claw.num_output_times = 10
claw.setplot = setplot
claw.keep_copy = True

claw.run()
gauge_pressure
t = claw.solution.state.gauge_data[0][:,0]
p = claw.solution.state.gauge_data[0][:,1]
print np.max(p)
