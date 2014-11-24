import numpy
import sys
import shelve


for file_name in [sys.argv[2]]:
    f = shelve.open(file_name, 'r')

    sol = f[f['last']]['xuser']
    num_ac = int(sys.argv[1])
    num_rt = sol['pax/flight_0'].shape[0] / num_ac

    print 'Design variable values'
    print '----------------------'
    print
    pax_flt = numpy.around(sol['pax/flight_0'], 2).reshape((num_rt, num_ac), order='F')
    flt_day = numpy.around(sol['flights/day_0'], 2).reshape((num_rt, num_ac), order='F')
    for iac in xrange(num_ac):
        print '%03i p/f:'%(iac), pax_flt[:, iac]
        print '%03i f/d:'%(iac), flt_day[:, iac]
        print

    if len(sys.argv) > 3 and sys.argv[3]:
        print
        copy = 0
        while 'h_pt_' + str(copy) in sol:
            print 'h_pt' + str(copy) + ':', numpy.around(sol['h_pt_' + str(copy)], 2)
            copy += 1

    index = int(f['last'])
    while 'funcs' not in f[str(index)]:
        index -= 1
    sol = f[str(index)]['funcs']
    print
    funcs = ['pax_con_0', 'ac_con_0', 'int_con_0']
    copy = 0
    while 'gamma_' + str(copy) in sol:
        funcs.append('Tmin_' + str(copy))
        funcs.append('Tmax_' + str(copy))
        funcs.append('h_i_' + str(copy))
        funcs.append('h_f_' + str(copy))
        funcs.append('gamma_' + str(copy))
        copy += 1
    for func in funcs:
	if func in sol:
            print func, numpy.min(sol[func]), numpy.max(sol[func])

    f.close()
