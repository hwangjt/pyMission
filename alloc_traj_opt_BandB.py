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
import timeit
import copy as cpy
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab

##########################
# COMMON DECLARATIONS
##########################

start_time = timeit.default_timer()
num_elem = 500
num_cp = 50
fileloc = open('./path.txt', 'r')
folder_path = fileloc.readlines()[0][:-1]
fileloc.close()

# define bounds for the flight path angle
gamma_lb = numpy.tan(-35.0 * (numpy.pi/180.0))/1e-1
gamma_ub = numpy.tan(35.0 * (numpy.pi/180.0))/1e-1

#execfile('problem_3rt_1ac_BWB.py')
execfile('problem_3rt_3ex1newBWB.py')

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
            IndVar('pax/flight', val=10, size=num_routes*num_ac),#Aset[Fsub_i].pax_flt),
            IndVar('flights/day', val=0.1, size=num_routes*num_ac),#Aset[Fsub_i].flt_day),
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
                                             num_existing_ac=num_existing_ac,
                                             num_routes=num_routes,
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
                                             num_existing_ac=num_existing_ac,
                                             num_routes=num_routes,
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
    main.compute(output=True)
    print
    print
    print
    #for ind in xrange(11):                                                                                                                                                                  
    #    main.compute_derivatives('rev', ('pax_con',0), ind, output=True)                                                                                                                    
    print main.vec['u']
    main.check_derivatives_all()
    print main.vec['u'].array.shape[0]


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

h_pts = []
for copy in xrange(num_routes*num_new_ac):
    h_pts.append(main.vec['u'](('h_pt',copy)))

##############################
#   Branch and Bound Algorithm
##############################
class Problem(object):
    'Empty container to hold the active sets'
    pass

#Parameters for branch and bound algorithm
num_int = num_ac*num_routes
_iter = 0
iter_max = 10000000
funCall = 0
eflag = 0
#U_best = numpy.inf  # If you are testing the code use this line, instead of U_best = 0
U_best = 0 #Do nothing: zero profit (worst case) [For airline allocation]
xopt = []
fopt = 0
can_x = []
can_F = []
ter_crit = 0
opt_cr = 0.03
strategy = 1
node_num = 1
tree = [1]

prob = Problem()
prob.flt_day_lb   = numpy.zeros(num_int)
prob.flt_day_ub   = 15*numpy.ones(num_int)
prob.b_F  = 0
prob.flt_day  = 0.1*numpy.ones(num_int)
prob.pax_flt  = 10*numpy.ones(num_int)
prob.node = node_num
prob.tree = tree

Aset = []
Aset.append(prob)

print '\nStarting Branch and Bound Algorithm....'
while len(Aset) > 0 and _iter <= iter_max:
    _iter = _iter + 1

    # Pick A subproblem
    if strategy == 1: # Strategy 1: Depth first search
        # Preference given to nodes with highest tree length
        max_tree = 0
        for ii in range(len(Aset)):
            if len(Aset[ii].tree) > max_tree:
                Fsub_i = ii
                max_tree = len(Aset[ii].tree)
    elif strategy == 2: # Strategy 2: Best first search
        # Preference given to nodes with the best objective value
        Fsub = numpy.inf
        for ii in range(len(Aset)):
            if Aset(ii).b_F < Fsub:
                Fsub_i = ii
                Fsub = Aset[ii].b_F

    if False:
        main.compute(output=True)
        print
        print
        print
        #for ind in xrange(11):
        #    main.compute_derivatives('rev', ('pax_con',0), ind, output=True)
        print main.vec['u']
        main.check_derivatives_all(elemsys_ids=['CL_tar', 'CT_tar', 'time', 'profit', 'pax_con', 'ac_con', 'int_con'])
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
            avail[iac] = 12 * ac_data['number', ac_name]

        main.set_initial_var_values()

        print '*********************'
        print 'Starting SNOPT....'
        print '*********************'
        funCall = funCall +1
        opt = Optimization(main)
        opt.add_objective('profit')
        opt.add_design_variable('pax/flight', value=Aset[Fsub_i].pax_flt,
                                lower=0, upper=pax_upper)
        opt.add_design_variable('flights/day', value=Aset[Fsub_i].flt_day,
                                lower=Aset[Fsub_i].flt_day_lb,
                                upper=Aset[Fsub_i].flt_day_ub)
        opt.add_constraint('pax_con', lower=0.1*demand, upper=demand,
                           get_jacs=main('pax_con').get_jacs)
        opt.add_constraint('ac_con', upper=avail)
        for copy in xrange(num_routes * num_new_ac):
            opt.add_design_variable(('h_pt', copy), value=h_pts[copy], #scale=5e-2,
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
        call(['mv', 'SNOPT_print.out', 'results/relaxed_'+str(funCall-1)+'_SNOPT.out'])
        call(['mv', 'hist.hst', folder_path+'/relaxed_'+str(funCall-1)+'_hist.hst'])
        call(['rm', 'SNOPT_summary.out'])

    Aset[Fsub_i].pax_flt = numpy.array(main.vec['u']('pax/flight'))
    Aset[Fsub_i].flt_day = numpy.array(main.vec['u']('flights/day'))
    Aset[Fsub_i].b_F = main.vec['u']('profit')[0]
    Aset[Fsub_i].eflag = opt.exit_flag

    # write status file
    status = numpy.ones((5*len(Aset), num_routes*num_ac)) * float('nan')
    for ind in xrange(len(Aset)):
        status[5*ind+0, :] = Aset[ind].flt_day_ub
        status[5*ind+1, :] = Aset[ind].flt_day
        status[5*ind+2, :] = Aset[ind].flt_day_lb
    numpy.savetxt('status.dat', status)

    if _iter == 1:
        flt_best_relax = numpy.array(Aset[Fsub_i].flt_day)
        pax_best_relax = numpy.array(Aset[Fsub_i].pax_flt)
        b_F_best_relax = Aset[Fsub_i].b_F

    if ((Aset[Fsub_i].eflag >= 1) and (Aset[Fsub_i].b_F < U_best)):
        aa = numpy.where(numpy.abs(numpy.round(Aset[Fsub_i].flt_day) - Aset[Fsub_i].flt_day) <= 1e-06)
        Aset[Fsub_i].flt_day[aa] = numpy.round(Aset[Fsub_i].flt_day[aa])
        if numpy.linalg.norm(Aset[Fsub_i].flt_day - numpy.round(Aset[Fsub_i].flt_day)) <= 1e-06:
            print '======================='
            print 'New solution found!'
            print '======================='
            # can_x = [can_x, Aset[Fsub_i].x_F]
            can_F = [can_F, Aset[Fsub_i].b_F]
            pax_flt_best = numpy.array(Aset[Fsub_i].pax_flt)
            flt_day_best = numpy.array(Aset[Fsub_i].flt_day)
            f_best = Aset[Fsub_i].b_F

            # Discard nodes within percentage of the tolerance gap of the best feasible solution (integer)
            # Optimal solution will be opt_cr% of the best feasible solution
            U_best = f_best/(1+numpy.sign(f_best)*opt_cr)
            del Aset[Fsub_i]  # Fathom by integrality
            ter_crit = 1

        else:
            # Branching
            x_ind_maxfrac = numpy.argmax(numpy.abs(Aset[Fsub_i].flt_day - numpy.round(Aset[Fsub_i].flt_day)))
            x_split = Aset[Fsub_i].flt_day[x_ind_maxfrac]
            print 'Branching at node: %d at x%d = %f' % (Aset[Fsub_i].node, x_ind_maxfrac+1, x_split)
            F_sub = [None, None]
            for jj in 0, 1:
                F_sub[jj] = cpy.deepcopy(Aset[Fsub_i])
                if jj == 0:
                    ub_new = numpy.floor(x_split)
                    if ub_new < F_sub[jj].flt_day_ub[x_ind_maxfrac]:
                        F_sub[jj].flt_day_ub[x_ind_maxfrac] = ub_new
                        F_sub[jj].flt_day[x_ind_maxfrac] = ub_new
                elif jj == 1:
                    lb_new = numpy.ceil(x_split)
                    if lb_new > F_sub[jj].flt_day_lb[x_ind_maxfrac]:
                        F_sub[jj].flt_day_lb[x_ind_maxfrac] = lb_new
                        F_sub[jj].flt_day[x_ind_maxfrac] = lb_new

                F_sub[jj].tree.append( jj+1)
                node_num = node_num + 1
                F_sub[jj].node = node_num
            del Aset[Fsub_i]  # Fathomed by branching
            Aset.extend(F_sub)
    else:
        del Aset[Fsub_i]  # Fathomed by infeasibility or bounds

print '\nTerminating algorithm...'
if ter_crit ==1:
    eflag = 1
    flt_day_opt = flt_day_best
    pax_flt_opt = pax_flt_best
    fopt = f_best
    print 'Optimal solution found!!'
    print flt_day_opt
    print pax_flt_opt
    print 'Profit (-ve): $',fopt
    print 'Total no.of integer solutions found:', len(can_F)
    print 'Total nodes evaluated: ', funCall
    print 'Relaxed solution'
    print flt_best_relax
    print pax_best_relax
    print b_F_best_relax
else:
    print 'No solution found!!'
    print 'Total nodes evaluated: ', funCall
    if _iter > iter_max:
        print 'Maximum number of iterations reached!'

end_time = timeit.default_timer()
tot_time = end_time - start_time
print 'Total execution time [min] ',tot_time/60
