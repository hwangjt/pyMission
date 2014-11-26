##########################
# misc parameters
##########################

misc_data = {}
misc_data['cost/fuel'] = 0.2431 # $/lb
misc_data['turnaround'] = 1.0 # hr
misc_data['weight/pax'] = 84 * 9.81 # N

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
ac_data['number', 'B737'] = 6#12
ac_data['number', 'B738'] = 4#8
ac_data['number', 'CRM'] = 4#8
ac_data['number', 'BWB'] = 4#8

# maintenance hours / block hours
ac_data['MH', 'B737'] = 0.936
ac_data['MH', 'B738'] = 0.948
ac_data['MH', 'CRM'] = 0.94
ac_data['MH', 'BWB'] = 0.94

# lb
ac_data['fuel', 'B737'] = numpy.array([3313.99, 9880.66,  12414.04, 13815.11,   16908.70,   18092.59,    26860.88,  27951.92,   28209.94,   29067.88,   29269.39])
ac_data['fuel', 'B738'] = numpy.array([5213.30, 10438.25, 13051.58, 14496.79,   17673.67,   18891.54,    27891.90,  29008.68,   29272.66,   30150.11,   30356.14])

# hr
ac_data['block time', 'B737'] = numpy.array([0.7926,    2.1123,    2.6028,    2.8778,    3.4597,    3.6813,    5.2480,    5.4356,    5.4798,    5.6296,    5.6653])
ac_data['block time', 'B738'] = numpy.array([1.0800,    2.0870,    2.5828,    2.8483,    3.4297,    3.6513,    5.2080,    5.4004,    5.4467,    5.5961,    5.6303])

# $
ac_data['flight cost no fuel', 'B737'] = numpy.array([5277.17,  11192.01,   13416.51,   14724.83,   17402.33,   18418.01,   25751.45,   26643.69,   26854.15,   27556.33,   27721.77])
ac_data['flight cost no fuel', 'B738'] = numpy.array([9334.81,  14633.87,   17230.35,   18760.10,   21864.49,   23041.04,   31487.07,   32520.73,   32765.80,   33572.59,   33760.35])
ac_data['flight cost no fuel', 'CRM'] = numpy.array([7200.,     12500.,     15100.,     16700.,     19600.,     20700.,     28500.,     29000.,     29000.,     30000.,     30000.])
ac_data['flight cost no fuel', 'BWB'] = numpy.array([7200.,     12500.,     15100.,     16700.,     19600.,     20700.,     28500.,     29000.,     29000.,     30000.,     30000.])

# $
ac_data['ticket price', 'B737'] = numpy.array([100.7932,    149.5979,   173.7622,   188.1758,   218.5249,   230.1823,   315.9648,   326.5451,   329.0442,   337.3819,   339.3462])
ac_data['ticket price', 'B738'] = numpy.array([114.9947,    160.1779,   185.5060,   200.4201,   231.3065,   243.0953,   328.7568,   339.3156,   341.8186,   350.0772,   352.0032])
ac_data['ticket price', 'CRM'] = numpy.array([107.,         155.,       179.,       194.,       225.,       237.,       322.,       332.,       335.,       342.,       344.,   ])
ac_data['ticket price', 'BWB'] = numpy.array([107.,         155.,       179.,       194.,       225.,       237.,       322.,       332.,       335.,       342.,       344.,   ])

##########################
# route parameters
##########################

rt_data = {}
rt_data['number'] = 11

# nm
rt_data['range'] = numpy.array([162,753,974,1094,1357,1455,2169,2249,2269,2337,2350])

#
rt_data['demand'] = numpy.array([41,1009,89,661,1041,358,146,97,447,194,263])

#
rt_data['num_elem'] = 500 * numpy.ones(11, int)
rt_data['num_cp'] = 50 * numpy.ones(11, int)
