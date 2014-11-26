from history import *
from framework import *
from optimization import *
from bsplines import *
from atmospherics import *
from coupled_analysis import *
from functionals import *
from aerodynamics import *
from propulsion import *
from aeroTripan import *
from alloc_func import *
import time
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab

##########################
# COMMON DECLARATIONS
##########################

num_elem = 500
num_cp = 50
fileloc = open('./path.txt', 'r')
folder_path = fileloc.readlines()[0][:-1]
fileloc.close()

# define bounds for the flight path angle
gamma_lb = numpy.tan(-35.0 * (numpy.pi/180.0))/1e-1
gamma_ub = numpy.tan(35.0 * (numpy.pi/180.0))/1e-1
                
#execfile('problem_11rt_2ac_v2.py')
execfile('problem_3rt_2ac.py')

############################
# INITIALIZE MISSION OBJECTS
############################

new_ac = {}

execfile('crm_params.py')
new_ac['CRM'] = params

execfile('bwb_params.py')
new_ac['BWB'] = params

num_routes = rt_data['number']
num_existing_ac = len(ac_data['existing_ac'])
num_new_ac = len(ac_data['new_ac'])
num_ac = num_existing_ac + num_new_ac

wt_pax = misc_data['weight/pax']

aero_surrogates = {}
for iac in xrange(num_new_ac):
    ac_name = ac_data['new_ac'][iac]
    params = new_ac[ac_name]

    surr = {}
    surr['CL'], surr['CD'], surr['CM'], surr['nums'] = setup_surrogate(params['surr'])

    aero_surrogates[ac_name] = surr

subsystems = []
subsystems.append(
    SerialSystem('allocation_param',
                 NL='NLN_GS',
                 LN='LIN_GS',
                 LN_ilimit=1,
                 NL_ilimit=1,
                 output=True,
                 subsystems=[
                     IndVar('pax/flight', val=10, size=num_routes*num_ac),
                     IndVar('flights/day', val=0.1, size=num_routes*num_ac),
                 ])
)

for irt in xrange(num_routes):
    num_elem = rt_data['num_elem'][irt]
    num_cp = rt_data['num_cp'][irt]
    mrange = rt_data['range'][irt] * 1.852

    lins = numpy.linspace(0, 1, num_cp)
    x_pt = mrange * 1e3 * \
           (1-numpy.cos(lins*numpy.pi))/2/1e6
    M_pt = numpy.ones(num_cp) * 0.82
    h_pt = 10 * numpy.sin(numpy.pi * x_pt / (mrange/1e3))
    v_pt = numpy.zeros(num_cp)

    # setup_MBI
    num_pts = num_elem + 1

    alt = numpy.linspace(0, 16, num_pts)
    x_dist = numpy.linspace(0, x_pt[-1], num_pts)/1e6

    arr = MBI.MBI(alt, [x_dist], [num_cp], [4])
    jac = arr.getJacobian(0, 0)
    jacd = arr.getJacobian(1, 0)

    c_arryx = x_pt
    d_arryx = jacd.dot(c_arryx)*1e6

    lins = numpy.linspace(0, num_pts-1, num_pts).astype(int)
    diag = scipy.sparse.csc_matrix((1.0/d_arryx,
                                    (lins,lins)))
    jace = diag.dot(jacd)

    jac_h = jac
    jac_gamma = jace

    for iac in xrange(num_new_ac):
        ac_name = ac_data['new_ac'][iac]
        params = new_ac[ac_name]
        copy = irt + iac * num_routes
        aero_surr = aero_surrogates[ac_name]

        subsystems.append(
            SerialSystem('mission',
                         copy=copy,
                         NL='NLN_GS',
                         LN='LIN_GS',
                         LN_ilimit=1,
                         NL_ilimit=1,
                         NL_atol=1e-9,
                         NL_rtol=1e-9,
                         output=True,
                         subsystems=[
                             SerialSystem('bsplines',
                                          copy=copy,
                                          NL='NLN_GS',
                                          LN='LIN_GS',
                                          LN_ilimit=1,
                                          NL_ilimit=1,
                                          output=True,
                                          subsystems=[
                                              IndVar('x_pt', copy=copy, val=x_pt, lower=0),
                                              IndVar('h_pt', copy=copy, val=h_pt, lower=0),
                                              IndVar('M_pt', copy=copy, val=M_pt, lower=0),
                                              IndVar('v_pt', copy=copy, val=v_pt, lower=0),
                                              SysXBspline('x', copy=copy,
                                                          num_elem=num_elem,
                                                          num_cp=num_cp,
                                                          x_init=x_pt,
                                                          x_0=numpy.linspace(0.0, x_pt[-1], num_elem+1),
                                                          jac_h=jac_h),
                                              SysHBspline('h', copy=copy,
                                                          num_elem=num_elem,
                                                          num_cp=num_cp,
                                                          x_init=None,
                                                          jac_h=jac_h),
                                              SysMVBspline('MV', copy=copy,
                                                           num_elem=num_elem,
                                                           num_cp=num_cp,
                                                           x_init=None,
                                                           jac_h=jac_h),
                                              SysGammaBspline('gamma', copy=copy,
                                                              num_elem=num_elem,
                                                              num_cp=num_cp,
                                                              x_init=None,
                                                              jac_gamma=jac_gamma),
                                          ]),
                             SerialSystem('atmospherics', copy=copy,
                                          NL='NLN_GS',
                                          LN='LIN_GS',
                                          LN_ilimit=1,
                                          NL_ilimit=1,
                                          output=True,
                                          subsystems=[
                                              SysSFC('SFC', copy=copy, num_elem=num_elem,
                                                     SFCSL=params['SFCSL']),
                                              SysTemp('Temp', copy=copy, num_elem=num_elem),
                                              SysRho('rho', copy=copy, num_elem=num_elem),
                                              SysSpeed('v', copy=copy, num_elem=num_elem,
                                                       v_specified=0),
                                              SysMach('M', copy=copy, num_elem=num_elem,
                                                      v_specified=0),
                                          ]),
                             SerialSystem('coupled_analysis', copy=copy,
                                          LN='KSP_PC',
                                          #PC='LIN_GS',
                                          LN_ilimit=50,
                                          GL_GS_ilimit=1,
                                          GL_NT_ilimit=8,
                                          PC_ilimit=2,
                                          GL_GS_rtol=1e-6,
                                          GL_GS_atol=1e-9,#10,
                                          GL_NT_rtol=1e-9,#14,
                                          GL_NT_atol=1e-9,#14,
                                          NL_rtol=1e-12,
                                          NL_atol=1e-12,
                                          LN_rtol=1e-20,#14,
                                          LN_atol=1e-14,#14,
                                          PC_rtol=1e-6,
                                          PC_atol=1e-10,
                                          output=True,
                                          subsystems=[
                                              SerialSystem('vert_eqlm', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='KSP_PC',
                                                           LN_ilimit=1,
                                                           NL_ilimit=1,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-10,
                                                           LN_atol=1e-10,
                                                           subsystems=[
                                                               SysCLTar('CL_tar', copy=copy,
                                                                        num_elem=num_elem,
                                                                        S=params['S'],
                                                                        wt_pax=wt_pax,
                                                                        ac_w=params['ac_w']),                                                             
                                                           ]),
                                              SerialSystem('tripan_alpha', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='LIN_GS',
                                                           LN_ilimit=18,
                                                           NL_ilimit=1,
                                                           PC_ilimit=2,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-6,
                                                           LN_atol=1e-6,
                                                           subsystems=[
                                                               SysTripanCLSurrogate('alpha', copy=copy, num_elem=num_elem, num=aero_surr['nums'], CL=aero_surr['CL']),
                                                           ]),
                                              SerialSystem('tripan_eta', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='LIN_GS',
                                                           LN_ilimit=18,
                                                           NL_ilimit=1,
                                                           PC_ilimit=2,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-6,
                                                           LN_atol=1e-6,
                                                           subsystems=[
                                                               SysTripanCMSurrogate('eta', copy=copy, num_elem=num_elem, num=aero_surr['nums'], CM=aero_surr['CM']),
                                                           ]),
                                              SerialSystem('tripan_drag', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='KSP_PC',
                                                           LN_ilimit=1,
                                                           NL_ilimit=1,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-10,
                                                           LN_atol=1e-10,
                                                           subsystems=[
                                                               SysTripanCDSurrogate('drag', copy=copy, num_elem=num_elem, num=aero_surr['nums'], CD=aero_surr['CD']),
                                                           ]),
                                              SerialSystem('hor_eqlm', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='KSP_PC',
                                                           LN_ilimit=1,
                                                           NL_ilimit=1,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-10,
                                                           LN_atol=1e-10,
                                                           subsystems=[
                                                               SysCTTar('CT_tar', copy=copy,
                                                                        num_elem=num_elem,
                                                                        S=params['S'],
                                                                        wt_pax=wt_pax,
                                                                        ac_w=params['ac_w']),
                                                           ]),
                                              SerialSystem('weight', copy=copy,
                                                           NL='NLN_GS',
                                                           LN='KSP_PC',
                                                           LN_ilimit=1,
                                                           NL_ilimit=1,
                                                           NL_rtol=1e-10,
                                                           NL_atol=1e-10,
                                                           LN_rtol=1e-10,
                                                           LN_atol=1e-10,
                                                           subsystems=[
                                                               SysFuelWeight('fuel_w', copy=copy,
                                                                             num_elem=num_elem,
                                                                             fuel_w_0=numpy.linspace(1.0, 0.0, num_elem+1),
                                                                             S=params['S']),
                                                           ]),
                                          ]),
                             SerialSystem('functionals', copy=copy,
                                          NL='NLN_GS',
                                          LN='LIN_GS',
                                          NL_ilimit=1,
                                          LN_ilimit=1,
                                          NL_rtol=1e-10,
                                          NL_atol=1e-10,
                                          LN_rtol=1e-10,
                                          LN_ato1=1e-10,
                                          output=True,
                                          subsystems=[
                                              SysTau('tau', copy=copy, num_elem=num_elem,
                                                     thrust_sl=params['thrust_sl'],
                                                     S=params['S']),
                                              SysFuelObj('fuelburn', copy=copy, num_elem=num_elem),
                                              SysHi('h_i', copy=copy,
                                                    num_elem=num_elem,
                                                    jac_h=jac_h),
                                              SysHf('h_f', copy=copy,
                                                    num_elem=num_elem,
                                                    jac_h=jac_h),
                                              SysTmin('Tmin', copy=copy, num_elem=num_elem),
                                              SysTmax('Tmax', copy=copy, num_elem=num_elem),
                                              SysMi('M_i', copy=copy),
                                              SysMf('M_f', copy=copy, num_elem=num_elem),
                                              SysBlockTime('time', copy=copy, num_elem=num_elem),
                                          ]),
                         ]),
        )

subsystems.append(
    SerialSystem('allocation_functionals',
                 NL='NLN_GS',
                 LN='LIN_GS',
                 LN_ilimit=1,
                 NL_ilimit=1,
                 output=True,
                 subsystems=[
                     Profit('profit',
                            misc_data=misc_data,
                            ac_data=ac_data,
                            rt_data=rt_data),
                     PaxCon('pax_con',
                            misc_data=misc_data,
                            ac_data=ac_data,
                            rt_data=rt_data),
                     AircraftCon('ac_con',
                            misc_data=misc_data,
                            ac_data=ac_data,
                            rt_data=rt_data),
                     IntegerCon('int_con',
                                size=num_routes*num_ac),
                 ])
)

main = SerialSystem('allocation',
                    NL='NLN_GS',
                    LN='LIN_GS',
                    LN_ilimit=1,
                    NL_ilimit=1,
                    NL_rtol=1e-6,
                    NL_atol=1e-9,
                    LN_rtol=1e-6,
                    LN_atol=1e-10,
                    output=True,
                    subsystems=subsystems).setup()

if False:
    main.compute(output=True)
    print
    print
    print
    for ind in xrange(11):
        main.compute_derivatives('rev', ('pax_con',0), ind, output=True)
    #print main.vec['u']
    main.check_derivatives_all()#elemsys_ids=['CL_tar', 'CT_tar', 'time', 'profit', 'pax_con', 'ac_con', 'int_con', 'Tmin', 'Tmax'])
    print main.vec['u'].array.shape[0]
else:
    pax_upper = numpy.zeros(num_routes * num_ac)
    upper = pax_upper.reshape((num_routes, num_ac), order='F')
    for iac in xrange(num_ac):
        if iac < num_existing_ac:
            ac_name = ac_data['existing_ac'][iac]
        else:
            inac = iac - num_existing_ac
            ac_name = ac_data['new_ac'][inac]
        for irt in xrange(num_routes):
            upper[irt, iac] = ac_data['capacity', ac_name]

    demand = numpy.zeros(num_routes)
    for irt in xrange(num_routes):
        demand[irt] = rt_data['demand'][irt]

    avail = numpy.zeros(num_ac)
    for iac in xrange(num_ac):
        if iac < num_existing_ac:
            ac_name = ac_data['existing_ac'][iac]
        else:
            inac = iac - num_existing_ac
            ac_name = ac_data['new_ac'][inac]
        avail[iac] = 24 * ac_data['number', ac_name]



    main.set_initial_var_values()


    if False:
        main('pax/flight').value = numpy.array(
            [[0,0,0],
             [0,0,0],
             [300,0,300],
             [400,400,400]]).flatten(order='F')
        main('flights/day').value = numpy.array(
            [[0,0,0],
             [0,0,0],
             [2,0,1],
             [1,1,1]]).flatten(order='F')


    for copy in xrange(num_routes * num_new_ac):
        opt = Optimization(main(('mission', copy)))
        opt.add_objective(('fuelburn', copy))
        opt.add_design_variable(('h_pt', copy), #scale=5e-2,
                                lower=0.0, upper=15.0)
        opt.add_constraint(('h_i', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_i', copy)).get_jacs, linear=True)
        opt.add_constraint(('h_f', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_f', copy)).get_jacs, linear=True)
        opt.add_constraint(('Tmin', copy), upper=0.0)
        opt.add_constraint(('Tmax', copy), upper=0.0)
        opt.add_constraint(('gamma', copy), lower=gamma_lb, upper=gamma_ub,
                           get_jacs=main(('gamma', copy)).get_jacs, linear=True)
        opt('SNOPT')
        call(['mv', 'SNOPT_print.out', 'results/traj_'+str(copy)+'_SNOPT.out'])
        call(['mv', 'hist.hst', folder_path+'/traj_'+str(copy)+'_hist.hst'])
        call(['rm', 'SNOPT_summary.out'])


    if False:
        main.compute(output=True)
        print 'Profit', main.vec['u']('profit')
        for i in xrange(6):
            print 'Fuel ' + str(i), main.vec['u'](('fuelburn',i)) * 1e-1 / 9.81 * 2.2
        exit()


    opt = Optimization(main)
    opt.add_objective('profit')
    opt.add_design_variable('pax/flight', lower=0, upper=pax_upper)
    opt.add_design_variable('flights/day', lower=0, upper=10)
    opt.add_constraint('pax_con', lower=0.1*demand, upper=demand,
                       get_jacs=main('pax_con').get_jacs)
    opt.add_constraint('ac_con', upper=avail)
    for copy in xrange(num_routes * num_new_ac):
        opt.add_design_variable(('h_pt', copy), #scale=5e-2,
                                lower=0.0, upper=15.0)
        opt.add_constraint(('h_i', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_i', copy)).get_jacs, linear=True)
        opt.add_constraint(('h_f', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_f', copy)).get_jacs, linear=True)
        opt.add_constraint(('Tmin', copy), upper=0.0, sys=main(('mission', copy)))
        opt.add_constraint(('Tmax', copy), upper=0.0, sys=main(('mission', copy)))
        opt.add_constraint(('gamma', copy), lower=gamma_lb, upper=gamma_ub,
                           get_jacs=main(('gamma', copy)).get_jacs, linear=True)
    opt('SNOPT')
    call(['mv', 'SNOPT_print.out', 'results/relaxed_SNOPT.out'])
    call(['mv', 'hist.hst', folder_path+'/relaxed_hist.hst'])
    call(['rm', 'SNOPT_summary.out'])



    opt = Optimization(main)
    opt.add_objective('profit')
    opt.add_design_variable('pax/flight', lower=0, upper=pax_upper)
    opt.add_design_variable('flights/day', lower=0, upper=10)
    opt.add_constraint('pax_con', lower=0.1*demand, upper=demand,
                       get_jacs=main('pax_con').get_jacs)
    opt.add_constraint('ac_con', upper=avail)
    opt.add_constraint('int_con', lower=0, upper=0,
                       get_jacs=main('int_con').get_jacs)
    for copy in xrange(num_routes * num_new_ac):
        opt.add_design_variable(('h_pt', copy), #scale=5e-2,
                                lower=0.0, upper=15.0)
        opt.add_constraint(('h_i', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_i', copy)).get_jacs, linear=True)
        opt.add_constraint(('h_f', copy), lower=0.0, upper=0.0,
                           get_jacs=main(('h_f', copy)).get_jacs, linear=True)
        opt.add_constraint(('Tmin', copy), upper=0.0, sys=main(('mission', copy)))
        opt.add_constraint(('Tmax', copy), upper=0.0, sys=main(('mission', copy)))
        opt.add_constraint(('gamma', copy), lower=gamma_lb, upper=gamma_ub,
                           get_jacs=main(('gamma', copy)).get_jacs, linear=True)

    opt('SNOPT')
    call(['mv', 'SNOPT_print.out', 'results/discrete_SNOPT.out'])
    call(['mv', 'hist.hst', folder_path+'/discrete_hist.hst'])
    call(['rm', 'SNOPT_summary.out'])
