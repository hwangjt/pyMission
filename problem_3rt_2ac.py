##########################
# misc parameters
##########################

misc_data = {}
misc_data['cost/fuel'] = 0.2431 # $/lb
misc_data['turnaround'] = 1.0 # hr
misc_data['weight/pax'] = 84 # kg

##########################
# a/c parameters
##########################

ac_data = {}
ac_data['existing_ac'] = ['B737', 'B738']
ac_data['new_ac'] = ['CRM', 'BWB']

#
ac_data['capacity', 'B737'] = 107
ac_data['capacity', 'B738'] = 122
ac_data['capacity', 'CRM'] = 300
ac_data['capacity', 'BWB'] = 400

#
ac_data['number', 'B737'] = 6
ac_data['number', 'B738'] = 4
ac_data['number', 'CRM'] = 4
ac_data['number', 'BWB'] = 4

# maintenance hours / block hours
ac_data['MH', 'B737'] = 0.936
ac_data['MH', 'B738'] = 0.948
ac_data['MH', 'CRM'] = 0.94
ac_data['MH', 'BWB'] = 0.94

# lb
ac_data['fuel', 'B737'] = numpy.array([24794.94, 18617.62, 12722.44])
ac_data['fuel', 'B738'] = numpy.array([25774.5, 19431.4, 13370.31])

# hr
ac_data['block time', 'B737'] = numpy.array([4.8891, 3.7773, 2.6617])
ac_data['block time', 'B738'] = numpy.array([4.8491, 3.7473, 2.6417])

# $
ac_data['flight cost no fuel', 'B737'] = numpy.array([24050.53, 18864.80, 13686.23])
ac_data['flight cost no fuel', 'B738'] = numpy.array([29528.54, 23558.28, 17544.29])
ac_data['flight cost no fuel', 'CRM'] = numpy.array([27000., 21000., 15500.])
ac_data['flight cost no fuel', 'BWB'] = numpy.array([27000., 21000., 15500.])

# $
ac_data['ticket price', 'B737'] = numpy.array([295.8682, 235.3350, 176.7478])
ac_data['ticket price', 'B738'] = numpy.array([308.7701, 248.2936, 188.5960])
ac_data['ticket price', 'CRM'] = numpy.array([300., 240., 180.])
ac_data['ticket price', 'BWB'] = numpy.array([300., 240., 180.])

##########################
# route parameters
##########################

rt_data = {}
rt_data['number'] = 3

# nm
rt_data['range'] = numpy.array([1999.6, 1498.8, 1000.8])

#
rt_data['demand'] = numpy.array([300.0, 700.0, 220.0])

#
rt_data['num_elem'] = numpy.array([251,252,253]) #500 * numpy.ones(3, int)
rt_data['num_cp'] = numpy.array([51,52,53]) #50 * numpy.ones(3, int)
#rt_data['num_elem'] = 25 * numpy.ones(3, int)
#rt_data['num_cp'] = 5 * numpy.ones(3, int)


