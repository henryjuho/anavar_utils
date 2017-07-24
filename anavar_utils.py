#!/usr/bin/env python


class Snp1ControlFile(object):

    def __init__(self):

        # create holders
        self.alg_opts = ''
        self.model_opts = ''
        self.dfe_opts = ''
        self.constraint_opts = ''

    def set_alg_opts(self, alg='NLOPT_LD_LBFGS', search=500):

        """
        sets algorithm options in control file
        :param alg: str
        :param search: int
        :return: sets algorithm string
        """

        algs = {'NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LD_VAR1', 'NLOPT_LD_VAR2', 'NLOPT_LN_NELDERMEAD'}

        if alg not in algs:
            raise ValueError('Unsupported algorithm')

        alg_string = ('[algorithm_commands]\n'
                      'search_algorithm: {}\n'
                      'maxeval: 100000\n'
                      'maxtime: 600\n'
                      'num_searches: {}\n'
                      'nnoimp: 1\n'
                      'maximp: 3\n'
                      'optional: true\n'
                      'size: 10000\n'
                      'key: 3\n'
                      'epsabs: 1e-20\n'
                      'epsrel: 1e-9\n'
                      'rftol: 1e-9\n'
                      '\n').format(alg, search)

        self.alg_opts = alg_string

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: list(tuple(list(float), int), ..)
        :param n: int
        :param snp_fold: bool
        :param dfe: str
        :param c: int
        :param theta_r: tuple(float, float)
        :param gamma_r: tuple(float, float)
        :param error_r: tuple(float, float)
        :param shape_r: tuple(float, float)
        :param scale_r: tuple(float, float)
        :param r_r: tuple(float, float)
        :return: NA
        """

        model_ctrls = ('[model_commands]\n'
                       'model: SNP_1\n'
                       'n: {n}\n'
                       'm: {m}\n'
                       'folded: {fold}\n'
                       'sfs: {sfs}\n'.format(n=n, m=sfs_m[0][1],
                                             fold=str(snp_fold).lower(),
                                             sfs=', '.join([str(x) for x in sfs_m[0][0]])))

        self.model_opts = model_ctrls

        if dfe == 'discrete':
            dfe_param = ('dfe: discrete\n'
                         'c: {c}\n'
                         'theta_range: {t1}, {t2}\n'
                         'gamma_range: {g1}, {g2}\n'
                         'e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(c=c,
                                         t1=theta_r[0], t2=theta_r[1],
                                         g1=gamma_r[0], g2=gamma_r[1],
                                         e1=error_r[0], e2=error_r[1])
        else:
            dfe_param = ('dfe: continuous\n'
                         'distribution: reflected_gamma\n'
                         'theta_range: {t1}, {t2}\n'
                         'shape_range: {sh1}, {sh2}\n'
                         'scale_range: {sc1}, {sc2}\n'
                         'e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(t1=theta_r[0], t2=theta_r[1],
                                         sh1=shape_r[0], sh2=shape_r[1],
                                         sc1=scale_r[0], sc2=scale_r[1],
                                         e1=error_r[0], e2=error_r[1])
        self.dfe_opts = dfe_param

    def set_constraint(self):
        constraint = 'constraint: none\n'

        self.constraint_opts = constraint

    def construct(self):

        # check data set
        if self.model_opts == '' or self.dfe_opts == '':
            raise ValueError('Data not set')

        # if data is set then check if other sections are set, if not set to defaults
        if self.alg_opts == '':
            self.set_alg_opts()
        if self.constraint_opts == '':
            self.set_constraint()

        control_contents = self.alg_opts + self.model_opts + self.dfe_opts + self.constraint_opts
        return control_contents


class Indel1ControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: list(tuple(list(float), int), ..)
        :param n: int
        :param snp_fold: bool
        :param dfe: str
        :param c: int
        :param theta_r: tuple(float, float)
        :param gamma_r: tuple(float, float)
        :param error_r: tuple(float, float)
        :param shape_r: tuple(float, float)
        :param scale_r: tuple(float, float)
        :param r_r: tuple(float, float)
        :return: NA
        """

        model_ctrls = ('[model_commands]\n'
                       'model: INDEL_1\n'
                       'n: {n}\n'
                       'm: {m1}\n'
                       'ins_sfs: {sfs1}\n'
                       'del_sfs: {sfs2}\n').format(n=n, m1=sfs_m[0][1],
                                                   sfs1=', '.join([str(x) for x in sfs_m[0][0]]),
                                                   sfs2=', '.join([str(x) for x in sfs_m[1][0]]))

        self.model_opts = model_ctrls

        if dfe == 'discrete':
            dfe_param = ('dfe: discrete\n'
                         'c: {c}\n'
                         'ins_theta_range: {t1}, {t2}\n'
                         'ins_gamma_range: {g1}, {g2}\n'
                         'ins_e_range: {e1}, {e2}\n'
                         'del_theta_range: {t1}, {t2}\n'
                         'del_gamma_range: {g1}, {g2}\n'
                         'del_e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(c=c,
                                         t1=theta_r[0], t2=theta_r[1],
                                         g1=gamma_r[0], g2=gamma_r[1],
                                         e1=error_r[0], e2=error_r[1])

        else:
            dfe_param = ('dfe: continuous\n'
                         'distribution: reflected_gamma\n'
                         'ins_theta_range: {t1}, {t2}\n'
                         'ins_shape_range: {sh1}, {sh2}\n'
                         'ins_scale_range: {sc1}, {sc2}\n'
                         'ins_e_range: {e1}, {e2}\n'
                         'del_theta_range: {t1}, {t2}\n'
                         'del_shape_range: {sh1}, {sh2}\n'
                         'del_scale_range: {sc1}, {sc2}\n'
                         'del_e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(t1=theta_r[0], t2=theta_r[1],
                                         sh1=shape_r[0], sh2=shape_r[1],
                                         sc1=scale_r[0], sc2=scale_r[1],
                                         e1=error_r[0], e2=error_r[1])

        self.dfe_opts = dfe_param


class GbgcControlFile(Snp1ControlFile):
    def __init__(self):

        Snp1ControlFile.__init__(self)

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: list(tuple(list(float), int), ..)
        :param n: int
        :param snp_fold: bool
        :param dfe: str
        :param c: int
        :param theta_r: tuple(float, float)
        :param gamma_r: tuple(float, float)
        :param error_r: tuple(float, float)
        :param shape_r: tuple(float, float)
        :param scale_r: tuple(float, float)
        :param r_r: tuple(float, float)
        :return: NA
        """

        model_ctrls = ('[model_commands]\n'
                       'model: gBGC_GLEMIN_EXTENDED_M1*\n'
                       'n: {n}\n'
                       'r_range: {r1}, {r2}\n'
                       '\n'
                       '[neutral_SNPs]\n'
                       'm: {m1}\n'
                       'sfs: {sfs1}\n'
                       'theta_range: {t1}, {t2}\n'
                       'gamma_range: {g1}, {g2}\n'
                       'e_range: {e1}, {e2}\n'
                       '\n'
                       '[ws_SNPs]\n'
                       'm: {m2}\n'
                       'sfs: {sfs2}\n'
                       'theta_range: {t1}, {t2}\n'
                       'gamma_range: {g1}, {g2}\n'
                       'e_range: {e1}, {e2}\n'
                       '\n'
                       '[sw_SNPs]\n'
                       'm: {m3}\n'
                       'sfs: {sfs3}\n'
                       'theta_range: {t1}, {t2}\n'
                       'gamma_range: {g1}, {g2}\n'
                       'e_range: {e1}, {e2}')

        model_ctrls = model_ctrls.format(n=n,
                                         r1=r_r[0], r2=r_r[1],
                                         m1=sfs_m[0][1], sfs1=', '.join([str(x) for x in sfs_m[0][0]]),
                                         m2=sfs_m[1][1], sfs2=', '.join([str(x) for x in sfs_m[1][0]]),
                                         m3=sfs_m[2][1], sfs3=', '.join([str(x) for x in sfs_m[2][0]]),
                                         t1=theta_r[0], t2=theta_r[1],
                                         g1=gamma_r[0], g2=gamma_r[1],
                                         e1=error_r[0], e2=error_r[1])

        self.model_opts = model_ctrls
        self.dfe_opts = '\n'


class IndelNeuSelControlFile(Snp1ControlFile):
    def __init__(self):

        Snp1ControlFile.__init__(self)

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: list(tuple(list(float), int), ..)
        :param n: int
        :param snp_fold: bool
        :param dfe: str
        :param c: int
        :param theta_r: tuple(float, float)
        :param gamma_r: tuple(float, float)
        :param error_r: tuple(float, float)
        :param shape_r: tuple(float, float)
        :param scale_r: tuple(float, float)
        :param r_r: tuple(float, float)
        :return: NA
        """

        model_ctrls = ('[model_commands]\n'
                       'model: neutralINDEL_vs_selectedINDEL\n'
                       'n: {n}\n'
                       'r_range: {r1}, {r2}\n'
                       'neu_indel_m: {m1}\n'
                       'neu_ins_sfs: {sfs1}\n'
                       'neu_del_sfs: {sfs2}\n'
                       'neu_ins_theta_range: {t1}, {t2}\n'
                       'neu_ins_e_range: {e1}, {e2}\n'
                       'neu_del_theta_range: {t1}, {t2}\n'
                       'neu_del_e_range: {e1}, {e2}\n'
                       'sel_indel_m: {m3}\n'
                       'sel_ins_sfs: {sfs3}\n'
                       'sel_del_sfs: {sfs4}\n')

        model_ctrls = model_ctrls.format(n=n,
                                         r1=r_r[0], r2=r_r[1],
                                         m1=sfs_m[0][1],
                                         sfs1=', '.join([str(x) for x in sfs_m[0][0]]),
                                         sfs2=', '.join([str(x) for x in sfs_m[1][0]]),
                                         m3=sfs_m[2][1],
                                         sfs3=', '.join([str(x) for x in sfs_m[2][0]]),
                                         sfs4=', '.join([str(x) for x in sfs_m[3][0]]),
                                         t1=theta_r[0], t2=theta_r[1],
                                         e1=error_r[0], e2=error_r[1])

        self.model_opts = model_ctrls

        if dfe == 'discrete':
            dfe_param = ('dfe: discrete\n'
                         'c: {c}\n'
                         'ins_theta_range: {t1}, {t2}\n'
                         'ins_gamma_range: {g1}, {g2}\n'
                         'ins_e_range: {e1}, {e2}\n'
                         'del_theta_range: {t1}, {t2}\n'
                         'del_gamma_range: {g1}, {g2}\n'
                         'del_e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(c=c,
                                         t1=theta_r[0], t2=theta_r[1],
                                         g1=gamma_r[0], g2=gamma_r[1],
                                         e1=error_r[0], e2=error_r[1])

        else:
            dfe_param = ('dfe: continuous\n'
                         'distribution: reflected_gamma\n'
                         'ins_theta_range: {t1}, {t2}\n'
                         'ins_shape_range: {sh1}, {sh2}\n'
                         'ins_scale_range: {sc1}, {sc2}\n'
                         'ins_e_range: {e1}, {e2}\n'
                         'del_theta_range: {t1}, {t2}\n'
                         'del_shape_range: {sh1}, {sh2}\n'
                         'del_scale_range: {sc1}, {sc2}\n'
                         'del_e_range: {e1}, {e2}\n')

            dfe_param = dfe_param.format(t1=theta_r[0], t2=theta_r[1],
                                         sh1=shape_r[0], sh2=shape_r[1],
                                         sc1=scale_r[0], sc2=scale_r[1],
                                         e1=error_r[0], e2=error_r[1])

        self.dfe_opts = dfe_param
