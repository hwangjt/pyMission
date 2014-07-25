"""
MISSION ANALYSIS/TRAJECTORY OPTIMIZATION
This is the runscript used for plotting the history of the trajectory 
optimization problem. The history plotting can be done simultaneously as
the trajectory optimization runscript is being ran. At the end of the
figure generation, a video of the history is also produced.
The mission analysis and trajectory optimization tool was developed by:
    Jason Kao*
    John Hwang*

* University of Michigan Department of Aerospace Engineering,
  Multidisciplinary Design Optimization lab
  mdolab.engin.umich.edu

copyright July 2014
"""

import numpy
import os
import time
from subprocess import call
from history import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab

# USER SPECIFIED INPUTS:

num_elem = 2000
num_cp = 200
x_range = 150.0
step = 1
initial_ind = 0
file_index = 0
video = True
folder_name = '/home/jason/Documents/Results/test-'
fuel_guess = 10000.0

# END USER SPECIFIED INPUTS

fig = matplotlib.pylab.figure(figsize=(18.0,8.0))
index = initial_ind
folder_name = folder_name + 'dist'+str(int(x_range))+'km-'\
    +str(num_cp)+'-'+str(num_elem)+'-'+str(file_index)+'/'
if not os.path.exists(folder_name):
    raise('ERROR: SPECIFIED CASE DOES NOT EXIST')

max_name = str(int(x_range))+'km-'+str(num_cp)+\
    '-'+str(num_elem)+'-maxmin.dat'
rnd = numpy.around
fplot = fig.add_subplot
file_name = '%ikm-%i-%i-%04i' % (int(x_range),
                                 num_cp,
                                 num_elem,
                                 index)
sleep = False

while ((not os.path.isfile(folder_name+max_name))
       or (os.path.isfile(folder_name+file_name+'.dat'))):
    if os.path.isfile(folder_name+'fig-'+file_name+'.png'):
        index += step
        file_name = '%ikm-%i-%i-%04i' % (int(x_range),
                                         num_cp,
                                         num_elem,
                                         index)
    else:
        if os.path.isfile(folder_name+file_name+'.dat'):
            if sleep == True:
                time.sleep(0.1)
                sleep = False

            [dist, altitude, speed, alpha, throttle, eta, fuel,
             rho, lift_c, drag_c, thrust, gamma, weight, temp,
             SFC] = numpy.loadtxt(folder_name+file_name+'.dat')
            dist = dist/1e3
            mach = speed / numpy.sqrt(1.4*288.0*temp)
            altitude *= 3.28
            speed *= 1.94
            fuel *= 0.225
            thrust *= 0.225
            weight *= 0.225

            print 'Printing fig: ', folder_name+file_name+'...'
            fig.clf()
            nr, nc = 4, 3

            values = [altitude/1e3, speed, eta, 
                      gamma, mach, alpha,
                      rho, throttle, lift_c,
                      fuel/1e3, thrust/1e3, drag_c]
            labels = ['Altitude (*10^3 ft)', 'TAS (knots)', 'Trim (deg)',
                      'Path Angle (deg)', 'Mach Number', 'AoA (deg)',
                      'Density (kg/m^3)', 'Throttle', 'C_L',
                      'Fuel wt. (10^3 lb)', 'Thrust (10^3 lb)', 'C_D']
            limits = [[-1, 51], [100, 600], [-10, 10],
                      [-32.0, 32.0], [0.05, 1.2], [-5, 10],
                      [0.0, 1.3], [-0.1, 1.1], [0.0, 0.8],
                      [-100.0/1e3, fuel_guess/1e3], [0.0, 250.0], [0.01*3, 0.05*3]]

            fplot = fig.add_subplot
            for i in xrange(12):
                fplot(nr, nc, i+1).plot(dist, values[i])
                fplot(nr, nc, i+1).set_ylabel(labels[i])
                fplot(nr, nc, i+1).set_xlim([-100.0, rnd(x_range, -2)+100.0])
                fplot(nr, nc, i+1).set_ylim(limits[i])

            fplot(nr, nc, 10).set_xlabel('Distance (km)')
            fplot(nr, nc, 11).set_xlabel('Distance (km)')
            fplot(nr, nc, 12).set_xlabel('Distance (km)')
            fig.savefig(folder_name+'fig-'+file_name+'.png')

            index += 1
            file_name = '%ikm-%i-%i-%04i' % (int(x_range),
                                             num_cp,
                                             num_elem,
                                             index)

        else:
            sleep = True
            time.sleep(0.1)

if video == True:
    file_name = '%ikm-%i-%i' % (int(x_range),
                                     num_cp,
                                     num_elem)
    call(["mencoder", "mf://"+folder_name+'fig-*.png', "-mf", 
          "fps=10:type=png", "-ovc", "x264", "-x264encopts", 
          "bitrate=15000", "-o", folder_name+file_name+".avi"])


