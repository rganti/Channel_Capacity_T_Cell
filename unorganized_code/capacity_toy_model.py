import numpy as np

# from toy_model import parameters

'''Defines a class that numerically calculates the capacity as a function of LsT/LfT and ks/kf.'''

parameters = {'kp': 0.1, 'koff': 0.05, 'koffs': 0.05, 'kon': 0.0022, 'kons': 0.1, 'kf': 0.2,
              'R': 30000.0, 'lfT': 10.0, 'M': 15, 'St': 10000.0}


class ComputeCapacity(object):

    def __init__(self, ls_multiplier=100.0, ks_multiplier=10.0):
        self.St = parameters['St']

        self.kp = parameters['kp']
        self.koff = self.koffs = parameters['koff']
        self.kon = parameters['kon']

        self.kons = parameters['kons']
        self.kf = parameters['kf']
        self.ks = ks_multiplier * self.kf
        self.R = parameters['R']
        self.lft = parameters['lfT']
        self.lst = ls_multiplier * self.lft

        self.gamma = self.kp / (self.kp + self.koff)
        self.alpha = self.kp / (self.kp + self.ks)
        self.beta = self.kp / (self.kp + self.kf)

        self.M = parameters['M']
        self.Nl_list = np.arange(1, self.M - 3 + 1)

        self.mu = np.log(self.lst)
        self.sigma = 1.0

        # bin_size = 0.01
        # sp = np.arange(0, max_bin, bin_size)
        # p_sp_list = []

    '''Function below calculates probability of sp given lf = 0.'''

    def compute_p_sp(self, sp, n):

        Nl = n
        f = self.kp / self.koff
        a = (self.kons / (self.kp + self.koffs)) * (self.kon / self.koff) * \
            self.gamma ** (self.M - 3 - Nl) * (self.alpha ** (Nl)) * self.R
        omega = (1 + self.R * (self.kon / self.kp) * self.alpha * ((1 - self.alpha ** (Nl)) / (1 - self.alpha)))

        l = omega * (sp / (f * a * self.St - sp * a * (1 + f)))
        dldsp = omega * (((f * a * self.St - sp * a * (1 + f)) + sp * a * (1 + f)) / (
                f * a * self.St - sp * a * (1 + f)) ** (2))

        p_sp = np.nan_to_num((1 / (self.sigma * l * np.sqrt(2 * np.pi))) *
                             np.exp(-(np.log(l) - self.mu) ** (2) / (2 * (self.sigma ** 2)))) * dldsp

        return p_sp

    '''Function below calculates probability of sp given lf > 0.'''

    def compute_p_sp_lf(self, sp, n):

        Nl = n
        f = self.kp / self.koff
        b = (self.kons / (self.kp + self.koffs)) * (self.kon / self.koff) * self.gamma ** (self.M - 3 - Nl) * self.R
        omega = (1 + self.R * (self.kon / self.kp) * self.alpha * ((1 - self.alpha ** (Nl)) / (1 - self.alpha)))

        lf = self.lft / (1 + self.R * (self.kon / self.kp) * self.beta * ((1 - self.beta ** (Nl)) / (1 - self.beta)))
        g = sp + sp * b * (self.beta ** (Nl)) * lf * (1 + f) - f * b * (self.beta ** (Nl)) * lf * self.St
        h = f * b * (self.alpha ** (Nl)) * self.St - sp * b * (self.alpha ** (Nl)) * (1 + f)

        l = omega * (g / h)

        gprime = 1 + b * (self.beta ** (Nl)) * lf * (1 + f)
        hprime = -b * (self.alpha ** (Nl)) * (1 + f)

        dldsp = omega * ((gprime * h - g * hprime) / (h ** 2))

        p_sp_lf = np.nan_to_num((1 / (self.sigma * l * np.sqrt(2 * np.pi))) *
                                np.exp(-(np.log(l) - self.mu) ** (2) / (2 * (self.sigma ** 2)))) * dldsp

        return p_sp_lf

    def compute_capacity(self):
        C_list = []

        for n in self.Nl_list:
            max_bin = 200000
            bin_size = 1.0
            p_O_integral = 0
            C = 0
            print("Nl = " + str(n))

            while p_O_integral < 0.99 or round(p_O_integral, 2) > 1.0:
                print("bin size = " + str(bin_size))
                sp = np.arange(0, max_bin, bin_size)
                print("num bins = " + str(len(sp)))

                p_sp_lf_0 = self.compute_p_sp(sp, n)
                p_sp_lf_greater_0 = self.compute_p_sp_lf(sp, n)

                p_O = 0.5 * (p_sp_lf_0 + p_sp_lf_greater_0)
                dsp = sp[1] - sp[0]

                p_O_integral = np.trapz(p_O, dx=dsp)
                print("p(O) integral = " + str(p_O_integral))
                term_1_c0 = 0.5 * p_sp_lf_greater_0 * np.nan_to_num(np.log2(p_sp_lf_greater_0 / p_O))
                term_2_d0 = 0.5 * p_sp_lf_0 * np.nan_to_num(np.log2(p_sp_lf_0 / p_O))

                C = np.trapz(term_1_c0 + term_2_d0, dx=dsp)
                print("C = " + str(C))

                if np.round(p_O_integral, 3) == np.round(C, 3):
                    print("C == P(O) integral: " + str(np.round(p_O_integral, 3) == np.round(C, 3)))
                    C = 1.00
                    print("New C " + str(C))
                    break

                bin_size /= 2

            C_list.append(C)

        return C_list
