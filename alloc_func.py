from __future__ import division
from framework import *
import numpy
import scipy.sparse


class IntegerCon(ExplicitSystem):

    def _declare(self):
        self.size = self.kwargs['size']

        size = self.size

        ind = range(size)
        self._declare_variable('int_con', size=size)
        self._declare_argument('flights/day', indices=ind)

    def apply_G(self):
        size = self.size

        con = self.vec['u']('int_con')
        flights = self.vec['p']('flights/day')

        con[:] = numpy.sin(numpy.pi * flights)

    def apply_dGdp(self, args):
        size = self.size

        con = self.vec['u']('int_con')
        flights = self.vec['p']('flights/day')

        dcon = self.vec['dg']('int_con')
        dflights = self.vec['dp']('flights/day')

        if self.mode == 'fwd':
            dcon[:] = 0.0
            if self.get_id('flights/day') in args:
                dcon[:] += numpy.pi * numpy.cos(numpy.pi * flights) * dflights
        elif self.mode == 'rev':
            dflights[:] = 0.0
            if self.get_id('flights/day') in args:
                dflights[:] += numpy.pi * numpy.cos(numpy.pi * flights) * dcon

    def get_jacs(self):
        lins = range(self.size)
        data = numpy.pi * numpy.cos(numpy.pi * self.vec['p']('flights/day'))
        jac = scipy.sparse.csr_matrix((data, (lins, lins)),
                                      shape=(self.size, self.size))
        return {'flights/day': jac}


class Profit(ExplicitSystem):
    """ computes the profit of the entire network """

    def _declare(self):
        self.misc_data = self.kwargs['misc_data']
        self.ac_data = self.kwargs['ac_data']
        self.rt_data = self.kwargs['rt_data']

        self.num_routes = self.rt_data['number']
        self.num_existing_ac = len(self.ac_data['existing_ac'])
        self.num_new_ac = len(self.ac_data['new_ac'])
        self.num_ac = self.num_existing_ac + self.num_new_ac

        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        ind = range(num_routes * num_ac)
        self._declare_variable('profit')
        self._declare_argument('pax/flight', indices=ind)
        self._declare_argument('flights/day', indices=ind)
        for inac in xrange(num_new_ac):
            for irt in xrange(num_routes):
                copy = irt + inac * num_routes
                self._declare_argument(('fuelburn', copy), indices=[0])

    def apply_G(self):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        cost_fuel = misc_data['cost/fuel'] * 2.2 / 9.81

        profit = self.vec['u']('profit')
        pax_flt = self.vec['p']('pax/flight').reshape((num_routes, num_ac), order='F')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')
        pvec = self.vec['p']

        profit[0] = 0
        for iac in xrange(num_ac):
            if iac < num_existing_ac:
                ac_name = ac_data['existing_ac'][iac]
            else:
                inac = iac - num_existing_ac
                ac_name = ac_data['new_ac'][inac]
            for irt in xrange(num_routes):
                if iac < num_existing_ac:
                    fuel = ac_data['fuel', ac_name][irt]
                else:
                    copy = irt + inac * num_routes
                    fuel = pvec(('fuelburn', copy)) * 1e5

                cost_nf = ac_data['flight cost no fuel', ac_name][irt]
                prc_pax = ac_data['ticket price', ac_name][irt]

                profit[0] -= prc_pax * pax_flt[irt, iac] * flt_day[irt, iac] / 1e6
                profit[0] += (cost_fuel * fuel + cost_nf) * flt_day[irt, iac] / 1e6

    def apply_dGdp(self, args):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        cost_fuel = misc_data['cost/fuel'] * 2.2 / 9.81

        profit = self.vec['u']('profit')
        pax_flt = self.vec['p']('pax/flight').reshape((num_routes, num_ac), order='F')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')
        pvec = self.vec['p']

        dprofit = self.vec['dg']('profit')
        dpax_flt = self.vec['dp']('pax/flight').reshape((num_routes, num_ac), order='F')
        dflt_day = self.vec['dp']('flights/day').reshape((num_routes, num_ac), order='F')
        dpvec = self.vec['dp']

        if self.mode == 'fwd':
            dprofit[0] = 0
            for iac in xrange(num_ac):
                if iac < num_existing_ac:
                    ac_name = ac_data['existing_ac'][iac]
                else:
                    inac = iac - num_existing_ac
                    ac_name = ac_data['new_ac'][inac]
                for irt in xrange(num_routes):
                    if iac < num_existing_ac:
                        fuel = ac_data['fuel', ac_name][irt]
                    else:
                        copy = irt + inac * num_routes
                        fuel = pvec(('fuelburn', copy)) * 1e5
                        dfuel = dpvec(('fuelburn', copy))
                        if self.get_id(('fuelburn', copy)) in args:
                            dprofit[0] += cost_fuel * flt_day[irt, iac] * dfuel * 1e5 / 1e6

                #profit[0] -= prc_pax * pax_flt[irt, iac] * flt_day[irt, iac] / 1e6
                #profit[0] += (cost_fuel * fuel + cost_nf) * flt_day[irt, iac] / 1e6

                    cost_nf = ac_data['flight cost no fuel', ac_name][irt]
                    prc_pax = ac_data['ticket price', ac_name][irt]

                    if self.get_id('pax/flight') in args:
                        dprofit[0] -= prc_pax * flt_day[irt, iac] * dpax_flt[irt, iac] / 1e6
                    if self.get_id('flights/day') in args:
                        dprofit[0] -= prc_pax * pax_flt[irt, iac] * dflt_day[irt, iac] / 1e6
                        dprofit[0] += (cost_fuel * fuel + cost_nf) * dflt_day[irt, iac] / 1e6
        elif self.mode == 'rev':
            dpax_flt[:, :] = 0.0
            dflt_day[:, :] = 0.0
            for iac in xrange(num_ac):
                if iac < num_existing_ac:
                    ac_name = ac_data['existing_ac'][iac]
                else:
                    inac = iac - num_existing_ac
                    ac_name = ac_data['new_ac'][inac]
                for irt in xrange(num_routes):
                    if iac < num_existing_ac:
                        fuel = ac_data['fuel', ac_name][irt]
                    else:
                        copy = irt + inac * num_routes
                        fuel = pvec(('fuelburn', copy)) * 1e5
                        dfuel = dpvec(('fuelburn', copy))
                        dfuel[0] = 0.0
                        if self.get_id(('fuelburn', copy)) in args:
                            dfuel[0] += cost_fuel * flt_day[irt, iac] * dprofit[0] * 1e5 / 1e6

                    cost_nf = ac_data['flight cost no fuel', ac_name][irt]
                    prc_pax = ac_data['ticket price', ac_name][irt]

                    if self.get_id('pax/flight') in args:
                        dpax_flt[irt, iac] -= prc_pax * flt_day[irt, iac] * dprofit[0] / 1e6
                    if self.get_id('flights/day') in args:
                        dflt_day[irt, iac] -= prc_pax * pax_flt[irt, iac] * dprofit[0] / 1e6
                        dflt_day[irt, iac] += (cost_fuel * fuel + cost_nf) * dprofit[0] / 1e6


class PaxCon(ExplicitSystem):

    def _declare(self):
        self.misc_data = self.kwargs['misc_data']
        self.ac_data = self.kwargs['ac_data']
        self.rt_data = self.kwargs['rt_data']

        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data

        self.num_routes = rt_data['number']
        self.num_existing_ac = len(ac_data['existing_ac'])
        self.num_new_ac = len(ac_data['new_ac'])
        self.num_ac = self.num_existing_ac + self.num_new_ac

        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        ind = range(num_routes * num_ac)
        self._declare_variable('pax_con', size=num_routes)
        self._declare_argument('pax/flight', indices=ind)
        self._declare_argument('flights/day', indices=ind)

        self.rows = numpy.zeros(num_routes * num_ac, int)
        self.cols = numpy.zeros(num_routes * num_ac, int)
        for irt in xrange(num_routes):
            for iac in xrange(num_ac):
                copy = irt + iac * num_routes
                self.rows[copy] = irt
                self.cols[copy] = copy

    def apply_G(self):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        pax_con = self.vec['u']('pax_con')
        pax_flt = self.vec['p']('pax/flight').reshape((num_routes, num_ac), order='F')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')

        pax_con[:] = 0.0
        for irt in xrange(self.num_routes):
            for iac in xrange(self.num_ac):
                pax_con[irt] += pax_flt[irt, iac] * flt_day[irt, iac]

    def apply_dGdp(self, args):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        pax_con = self.vec['u']('pax_con')
        pax_flt = self.vec['p']('pax/flight').reshape((num_routes, num_ac), order='F')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')

        dpax_con = self.vec['dg']('pax_con')
        dpax_flt = self.vec['dp']('pax/flight').reshape((num_routes, num_ac), order='F')
        dflt_day = self.vec['dp']('flights/day').reshape((num_routes, num_ac), order='F')

        if self.mode == 'fwd':
            dpax_con[:] = 0.0
            for irt in xrange(num_routes):
                for iac in xrange(num_ac):
                    if self.get_id('pax/flight') in args:
                        dpax_con[irt] += flt_day[irt, iac] * dpax_flt[irt, iac]
                    if self.get_id('flights/day') in args:
                        dpax_con[irt] += pax_flt[irt, iac] * dflt_day[irt, iac]
        elif self.mode == 'rev':
            dpax_flt[:, :] = 0.0
            dflt_day[:, :] = 0.0
            for irt in xrange(num_routes):
                for iac in xrange(num_ac):
                    if self.get_id('pax/flight') in args:
                        dpax_flt[irt, iac] += flt_day[irt, iac] * dpax_con[irt]
                    if self.get_id('flights/day') in args:
                        dflt_day[irt, iac] += pax_flt[irt, iac] * dpax_con[irt]

    def get_jacs(self):
        num_routes, num_ac = self.num_routes, self.num_ac
        rows, cols = self.rows, self.cols

        data = self.vec['p']('pax/flight')
        jac_f = scipy.sparse.csr_matrix((data, (rows, cols)),
                                        shape=(num_routes, num_routes*num_ac))

        data = self.vec['p']('flights/day')
        jac_p = scipy.sparse.csr_matrix((data, (rows, cols)),
                                        shape=(num_routes, num_routes*num_ac))
        return {'flights/day': jac_f, 'pax/flight': jac_p}


class AircraftCon(ExplicitSystem):

    def _declare(self):
        self.misc_data = self.kwargs['misc_data']
        self.ac_data = self.kwargs['ac_data']
        self.rt_data = self.kwargs['rt_data']

        self.num_routes = self.rt_data['number']
        self.num_existing_ac = len(self.ac_data['existing_ac'])
        self.num_new_ac = len(self.ac_data['new_ac'])
        self.num_ac = self.num_existing_ac + self.num_new_ac

        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        ind = range(self.num_routes * num_ac)
        self._declare_variable('ac_con', size=num_ac)
        self._declare_argument('flights/day', indices=ind)
        for inac in xrange(num_new_ac):
            for irt in xrange(num_routes):
                copy = irt + inac * num_routes
                self._declare_argument(('time', copy), indices=[0])

    def apply_G(self):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        turnaround = misc_data['turnaround']

        ac_con = self.vec['u']('ac_con')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')
        pvec = self.vec['p']

        ac_con[:] = 0.0
        for iac in xrange(num_ac):
            if iac < num_existing_ac:
                ac_name = ac_data['existing_ac'][iac]
            else:
                inac = iac - num_existing_ac
                ac_name = ac_data['new_ac'][inac]
            maintenance = ac_data['MH', ac_name]
            for irt in xrange(num_routes):
                if iac < num_existing_ac:
                    time = ac_data['block time', ac_name][irt]
                else:
                    copy = irt + inac * num_routes
                    time = pvec(('time', copy)) / 3.6e3 * 1e4

                ac_con[iac] += flt_day[irt, iac] * \
                               (time * (1 + maintenance) + turnaround)

    def apply_dGdp(self, args):
        misc_data, ac_data, rt_data = self.misc_data, self.ac_data, self.rt_data
        num_routes, num_ac = self.num_routes, self.num_ac
        num_existing_ac, num_new_ac = self.num_existing_ac, self.num_new_ac

        turnaround = misc_data['turnaround']

        ac_con = self.vec['u']('ac_con')
        flt_day = self.vec['p']('flights/day').reshape((num_routes, num_ac), order='F')
        pvec = self.vec['p']

        dac_con = self.vec['dg']('ac_con')
        dflt_day = self.vec['dp']('flights/day').reshape((num_routes, num_ac), order='F')
        dpvec = self.vec['dp']

        if self.mode == 'fwd':
            dac_con[:] = 0.0
            for iac in xrange(num_ac):
                if iac < num_existing_ac:
                    ac_name = ac_data['existing_ac'][iac]
                else:
                    inac = iac - num_existing_ac
                    ac_name = ac_data['new_ac'][inac]
                maintenance = ac_data['MH', ac_name]
                num_avail = ac_data['number', ac_name]
                for irt in xrange(num_routes):
                    if iac < num_existing_ac:
                        time = ac_data['block time', ac_name][irt]
                    else:
                        copy = irt + inac * num_routes
                        time = pvec(('time', copy)) / 3.6e3 * 1e4
                        dtime = dpvec(('time', copy))
                        if self.get_id(('time', copy)) in args:
                            dac_con[iac] += flt_day[irt, iac] * (1 + maintenance) * dtime / 3.6e3 * 1e4

                    if self.get_id('flights/day') in args:
                        dac_con[iac] += (time * (1 + maintenance) + turnaround) * \
                            dflt_day[irt, iac]
        elif self.mode == 'rev':
            dflt_day[:, :] = 0.0
            for iac in xrange(num_ac):
                if iac < num_existing_ac:
                    ac_name = ac_data['existing_ac'][iac]
                else:
                    inac = iac - num_existing_ac
                    ac_name = ac_data['new_ac'][inac]
                maintenance = ac_data['MH', ac_name]
                num_avail = ac_data['number', ac_name]
                for irt in xrange(num_routes):
                    if iac < num_existing_ac:
                        time = ac_data['block time', ac_name][irt]
                    else:
                        copy = irt + inac * num_routes
                        time = pvec(('time', copy)) / 3.6e3 * 1e4
                        dtime = dpvec(('time', copy))
                        dtime[0] = 0.0
                        if self.get_id(('time', copy)) in args:
                           dtime[0] += flt_day[irt, iac] * (1 + maintenance) * dac_con[iac] / 3.6e3 * 1e4

                    if self.get_id('flights/day') in args:
                        dflt_day[irt, iac] += (time * (1 + maintenance) + turnaround) * \
                            dac_con[iac]
