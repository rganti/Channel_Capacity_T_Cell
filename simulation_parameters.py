class DefineRegion(object):
    def __init__(self):
        self.x = self.y = 20.0
        self.depth = 1.0
        self.subvolume_edge = 1.0
        self.cytosol_depth = 10.0
        self.num_chambers = (self.x / self.subvolume_edge) * (self.y / self.subvolume_edge) * (
                    self.depth / self.subvolume_edge)
        self.num_cytosol_chambers = self.num_chambers * self.cytosol_depth

    def define_region(self, f):
        f.write('region World box width 1 height 1 depth 1\n')
        f.write('subvolume edge 1\n\n')

    def define_membrane_region(self, f):
        f.write('region Plasma \n \t box width {0} height {1} depth {2}\n\n'.format(self.x, self.y, self.depth))
        f.write('region Cytosol \n ')
        f.write('\t move 0 0 -{0} \n'.format(self.cytosol_depth / 2.0))
        f.write('\t\t box width {0} height {1} depth {2}\n\n'.format(self.x, self.y, self.cytosol_depth))
        f.write('subvolume edge {0}\n\n'.format(self.subvolume_edge))


class InitialConcentrations(object):
    def __init__(self):
        self.r_0 = 305  # 00
        self.lck_0 = 305  # 00
        self.zap_0 = 80000  # 00

        self.lat_0 = 150000
        self.sos_0 = 1000

        self.ras_gdp_0 = 1000
        self.ras_gap_0 = 10


class BindingParameters(object):
    def __init__(self):
        # Initial Conditions
        self.initial = InitialConcentrations()

        # Zeroth Cycle ligand binding
        self.k_L_on = 0.0052
        self.k_foreign_off = 0.2
        self.k_self_off = 2.0  # 10.0 * self.k_foreign_off

        # First Cycle Lck binding
        self.on_rate = 2.0

        self.k_lck_on_R_pmhc = (self.on_rate / self.initial.lck_0)  # * 100.0
        self.k_lck_off_R_pmhc = self.k_foreign_off / 40.0

        self.k_lck_on_R = 0.0  # 1.0 / 100000.0
        self.k_lck_off_R = 20.0

        # Second Cycle Lck phosphorylation
        self.k_p_on_R_pmhc = self.on_rate
        self.k_p_off_R_pmhc = self.k_foreign_off / 40.0  # self.k_foreign_off / 10.0

        self.k_lck_on_RP = 1.0 / 100000.0
        self.k_lck_off_RP = 20.0

        self.k_p_on_R = 0.0  # 1.0 / 100000.0
        self.k_p_off_R = 20.0

        # Third Cycle Zap binding
        self.k_zap_on_R_pmhc = (self.on_rate / self.initial.zap_0)  # * 100.0
        self.k_zap_off_R_pmhc = self.k_lck_off_R_pmhc

        self.k_lck_on_zap_R = self.k_lck_on_RP
        self.k_lck_off_zap_R = 20.0

        self.k_zap_on_R = 1.0 / 100000.0
        self.k_zap_off_R = 20.0

        # Fourth Cycle phosphorylate zap
        self.k_p_on_zap_species = self.on_rate  # / 10.0
        self.k_p_off_zap_species = self.k_p_off_R_pmhc

        self.k_lck_on_zap_p = self.k_lck_on_RP
        self.k_lck_off_zap_p = self.k_lck_off_zap_R

        # Fifth Negative Feedback Loop
        self.k_negative_loop = 0.05  # 5.0
        self.k_lcki = 0.01

        # Sixth LAT on
        self.k_lat_on_species = self.on_rate / self.initial.lat_0
        self.k_lat_off_species = self.k_lck_off_R_pmhc

        self.k_lat_on_rp_zap = self.k_zap_on_R
        self.k_lat_off_rp_zap = 20.0

        # Seventh first LAT phosphorylation
        self.k_p_lat_1 = self.on_rate

        # Eighth second LAT phosphorylation
        self.k_p_lat_on_species = self.on_rate / 10.0  # / 10.0  # self.k_p_on_zap_species / 100
        self.k_p_lat_off_species = self.k_p_off_R_pmhc

        self.k_p_on_lat = self.k_p_on_R
        self.k_p_off_lat = 20.0

        # Eighth Sos on
        self.k_sos_on = 0.001
        self.k_sos_off = 0.005

        self.multiplier = 10.0
        # Ninth Sos RasGDP and RasGTP
        self.k_sos_on_rgdp = 0.0024 * self.multiplier
        self.k_sos_off_rgdp = 3.0 * self.multiplier

        self.k_sos_on_rgtp = 0.0022 * self.multiplier
        self.k_sos_off_rgtp = 0.4 * self.multiplier

        # Tenth
        # sos_rgtp + rgdp
        self.k_rgdp_on_sos_rgtp = 0.001 * self.multiplier
        self.k_rgdp_off_sos_rgtp = 0.1 * self.multiplier
        self.k_cat_3 = 0.038 * 0.9 * self.multiplier

        # sos_rgdp + rgdp
        self.k_rgdp_on_sos_rgdp = 0.0014 * self.multiplier
        self.k_rgdp_off_sos_rgdp = 1.0 * self.multiplier
        self.k_cat_4 = 0.003 * self.multiplier

        # rgap + rgtp -> rgdp
        self.k_rgap_on_rgtp = 0.0348 * self.multiplier
        self.k_rgap_off_rgtp = 0.2 * self.multiplier
        self.k_cat_5 = 0.1 * self.multiplier

        # Eighth LAT -> final product
        self.k_product_on = 0.008
        self.k_product_off = self.k_p_off_R_pmhc

        # Ninth: Positive Feedback Loop
        self.k_positive_loop = 0.0025

        # # Eighth Grb on
        # self.k_grb_on_species = self.k_lck_on_R_pmhc
        # self.k_grb_off_species = self.k_lck_off_R_pmhc
        #
        # # Ninth Sos on
        # self.k_sos_on_species = self.k_lck_on_R_pmhc
        # self.k_sos_off_species = self.k_lck_off_R_pmhc
        #
        # # Tenth Ras-GDP on
        # self.k_ras_on_species = self.k_lck_on_R_pmhc
        # self.k_ras_off_species = self.k_lck_off_R_pmhc


class MembraneInitialConcentrations(InitialConcentrations):
    def __init__(self):
        InitialConcentrations.__init__(self)

        # self.r_0 = 300  # 75
        # self.lck_0 = 300  # 2705
        # self.zap_0 = 12000  # 12180
        self.lat_0 = 1522


class DiffusionRates(object):
    def __init__(self):
        self.region = DefineRegion()
        self.d_r = self.d_l = self.d_rl = 0.13 * self.region.num_chambers
        self.d_lck = 0.085 * self.region.num_chambers

        self.d_zap = 10.0


class MembraneBindingParameters(object):
    def __init__(self):
        # Initial Conditions
        self.initial = MembraneInitialConcentrations()
        self.region = DefineRegion()
        self.rates = BindingParameters()

        # Zeroth Cycle ligand binding
        self.k_L_on = self.rates.k_L_on * self.region.num_chambers
        self.k_foreign_off = self.rates.k_foreign_off
        self.k_self_off = self.rates.k_self_off

        # Lck binding
        self.on_rate = self.rates.on_rate

        self.k_lck_on_R_pmhc = (self.on_rate / self.initial.lck_0) * self.region.num_chambers
        self.k_lck_off_R_pmhc = self.rates.k_lck_off_R_pmhc

        self.k_lck_on_R = self.rates.k_lck_on_R * self.region.num_chambers
        self.k_lck_off_R = self.rates.k_lck_off_R

        # ITAM phosphorylation

        self.k_p_on_R_pmhc = self.on_rate
        self.k_p_off_R_pmhc = self.rates.k_p_off_R_pmhc

        self.k_lck_on_RP = self.rates.k_lck_on_RP * self.region.num_chambers
        self.k_lck_off_RP = self.rates.k_lck_off_RP

        self.k_p_on_R = self.rates.k_p_on_R
        self.k_p_off_R = self.rates.k_p_off_R

        # Zap 70 binding to ITAMs

        self.k_zap_on_R_pmhc = (self.on_rate / self.initial.zap_0) * self.region.num_chambers  # * 100.0
        self.k_zap_off_R_pmhc = self.rates.k_zap_off_R_pmhc

        self.k_lck_on_zap_R = self.rates.k_lck_on_zap_R * self.region.num_chambers
        self.k_lck_off_zap_R = self.rates.k_lck_off_zap_R

        self.k_zap_on_R = self.rates.k_zap_on_R * self.region.num_chambers
        self.k_zap_off_R = self.rates.k_zap_off_R
