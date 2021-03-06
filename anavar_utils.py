#!/usr/bin/env python

from __future__ import print_function
import re
import os
from anavar_alpha import alpha


class Snp1ControlFile(object):

    def __init__(self):

        # create holders
        self.alg_opts = ''
        self.model_opts = ''
        self.dfe_opts = ''
        self.constraint_opts = ''
        self.sfs = ['SNP']
        self.dfe_optional_opts = ''
        self.constraints = ['none', 'no_pol_error', 'equal_pol_error', 'equal_mutation_rate']

    def _check_sfs_m_in(self, sfs_m):

        """
        checks that all required sfs are given in input
        :param sfs_m: dict
        :return: NA
        """

        sfs_keys = sfs_m.keys()

        sfs_not_given = []
        for key in self.sfs:
            if key not in sfs_keys:
                sfs_not_given.append(key)

        if len(sfs_not_given) > 0:
            raise KeyError('Missing the following SFS in input: {}'.format(','.join(sfs_not_given)))
            
    def set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3,
                     optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=()):

        """
        sets algorithm options in control file
        :param alg: str
        :param maxeval: int
        :param maxtime: int
        :param search: int
        :param nnoimp: int
        :param maximp: int
        :param optional: bool
        :param size: int
        :param key: int
        :param epsabs: float
        :param epsrel: float
        :param rftol: float
        :param init: tuple
        :return: sets algorithm string
        """

        algs = {'NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LD_VAR1', 'NLOPT_LD_VAR2',
                'NLOPT_LN_NELDERMEAD', 'NLOPT_LD_TNEWTON_PRECOND_RESTART'}

        if alg not in algs:
            raise ValueError('Unsupported algorithm')

        alg_string = ('[algorithm_commands]\n'
                      'search_algorithm: {}\n'
                      'maxeval: {}\n'
                      'maxtime: {}\n'
                      'num_searches: {}\n'
                      'nnoimp: {}\n'
                      'maximp: {}\n').format(alg, maxeval, maxtime, search, nnoimp, maximp)

        if optional:

            # set init string - new in 1.4
            if len(init) == 0:
                init = 'random'
            else:
                init = 'given[{}]'.format('\t'.join([str(x) for x in init]))

            alg_string += ('optional: true\n'
                           'size: {}\n'
                           'key: {}\n'
                           'epsabs: {}\n'
                           'epsrel: {}\n'
                           'rftol: {}\n'
                           'init: {}\n'
                           '\n').format(size, key, epsabs, epsrel, rftol, init)
        else:
            alg_string += 'optional: false\n\n'

        self.alg_opts = alg_string

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

        model_ctrls = ('[model_commands]\n'
                       'model: SNP_1\n'
                       'n: {n}\n'
                       'm: {m}\n'
                       'folded: {fold}\n'
                       'sfs: {sfs}\n'.format(n=n, m=sfs_m['SNP'][1],
                                             fold=str(snp_fold).lower(),
                                             sfs=', '.join([str(x) for x in sfs_m['SNP'][0]])))

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

            if self.dfe_optional_opts == '':
                self.set_dfe_optional_opts()

        self.dfe_opts = dfe_param

    def set_constraint(self, constraint):

        """
        sets constraint
        :param constraint: str
        :return: sets str
        """

        if constraint not in self.constraints:
            raise ValueError('Not a recognised constraint!')

        constraint_str = 'constraint: {}\n'.format(constraint)

        self.constraint_opts = constraint_str

    def set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-5):

        opt_string = 'optional: {}\n'.format(str(optional).lower())
        if optional is True:
            opt_string += ('fraction: {}\n'
                           'delta: {}\n'
                           'degree: {}\n'
                           '').format(fraction, delta, degree)

        self.dfe_optional_opts = opt_string

    def construct(self):

        # check data set
        if self.model_opts == '' or self.dfe_opts == '':
            raise ValueError('Data not set')

        # if data is set then check if other sections are set, if not set to defaults
        if self.alg_opts == '':
            self.set_alg_opts()
        if self.constraint_opts == '':
            self.set_constraint('none')

        control_contents = (self.alg_opts + self.model_opts +
                            self.dfe_opts + self.constraint_opts +
                            self.dfe_optional_opts)
        return control_contents


class Indel1ControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

        self.sfs = ['INS', 'DEL']
        self.constraints = ['none', 'neutral']

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

        model_ctrls = ('[model_commands]\n'
                       'model: INDEL_1\n'
                       'n: {n}\n'
                       'm: {m1}\n'
                       'ins_sfs: {sfs1}\n'
                       'del_sfs: {sfs2}\n').format(n=n, m1=sfs_m['INS'][1],
                                                   sfs1=', '.join([str(x) for x in sfs_m['INS'][0]]),
                                                   sfs2=', '.join([str(x) for x in sfs_m['DEL'][0]]))

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

            if self.dfe_optional_opts == '':
                self.set_dfe_optional_opts()

        self.dfe_opts = dfe_param


class GbgcControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

        self.sfs = ['neutral_SNPs', 'ws_SNPs', 'sw_SNPs']
        self.constraints = ['none', 'M0', 'M1', 'M0*', 'M1*']

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

        model_ctrls = ('[model_commands]\n'
                       'model: gBGC_EXTENDED_M1*\n'
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
                                         m1=sfs_m['neutral_SNPs'][1],
                                         sfs1=', '.join([str(x) for x in sfs_m['neutral_SNPs'][0]]),
                                         m2=sfs_m['ws_SNPs'][1],
                                         sfs2=', '.join([str(x) for x in sfs_m['ws_SNPs'][0]]),
                                         m3=sfs_m['sw_SNPs'][1],
                                         sfs3=', '.join([str(x) for x in sfs_m['sw_SNPs'][0]]),
                                         t1=theta_r[0], t2=theta_r[1],
                                         g1=gamma_r[0], g2=gamma_r[1],
                                         e1=error_r[0], e2=error_r[1])

        self.model_opts = model_ctrls
        self.dfe_opts = '\n'


class IndelNeuSelControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

        self.sfs = ['neutral_INS', 'neutral_DEL', 'selected_INS', 'selected_DEL']
        self.constraints = ['none', 'equal_mutation_rate']

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

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
                                         m1=sfs_m['neutral_INS'][1],
                                         sfs1=', '.join([str(x) for x in sfs_m['neutral_INS'][0]]),
                                         sfs2=', '.join([str(x) for x in sfs_m['neutral_DEL'][0]]),
                                         m3=sfs_m['selected_INS'][1],
                                         sfs3=', '.join([str(x) for x in sfs_m['selected_INS'][0]]),
                                         sfs4=', '.join([str(x) for x in sfs_m['selected_DEL'][0]]),
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

            if self.dfe_optional_opts == '':
                self.set_dfe_optional_opts()

        self.dfe_opts = dfe_param


class IndelSnpNeuSelControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

        self.sfs = ['neutral_SNP', 'selected_INS', 'selected_DEL']
        self.constraints = ['none']

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

        model_ctrls = ('[model_commands]\n'
                       'model: neutralSNP_vs_selectedINDEL\n'
                       'n: {n}\n'
                       'r_range: {r1}, {r2}\n'
                       'neu_snp_m: {m1}\n'
                       'neu_snp_sfs: {sfs1}\n'
                       'neu_snp_theta_range: {t1}, {t2}\n'
                       'neu_snp_e_range: {e1}, {e2}\n'
                       'sel_indel_m: {m3}\n'
                       'sel_ins_sfs: {sfs3}\n'
                       'sel_del_sfs: {sfs4}\n')

        model_ctrls = model_ctrls.format(n=n,
                                         r1=r_r[0], r2=r_r[1],
                                         m1=sfs_m['neutral_SNP'][1],
                                         sfs1=', '.join([str(x) for x in sfs_m['neutral_SNP'][0]]),
                                         m3=sfs_m['selected_INS'][1],
                                         sfs3=', '.join([str(x) for x in sfs_m['selected_INS'][0]]),
                                         sfs4=', '.join([str(x) for x in sfs_m['selected_DEL'][0]]),
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

            if self.dfe_optional_opts == '':
                self.set_dfe_optional_opts()

        self.dfe_opts = dfe_param


class SNPNeuSelControlFile(Snp1ControlFile):

    def __init__(self):
        Snp1ControlFile.__init__(self)

        self.sfs = ['neutral_SNP', 'selected_SNP']
        self.constraints = ['none', 'equal_mutation_rate', 'equal_mutation_rate_AND_no_pol_error']

    def set_data(self, sfs_m, n, snp_fold=False,
                 dfe='discrete', c=1,
                 theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                 shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        sets model and dfe commands in control file
        :param sfs_m: dict{str: tuple(list(float), int), ..}
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

        self._check_sfs_m_in(sfs_m)

        model_ctrls = ('[model_commands]\n'
                       'model: neutralSNP_vs_selectedSNP\n'
                       'n: {n}\n'
                       'folded: {fold}\n'
                       'r_range: {r1}, {r2}\n'
                       'neu_m: {m1}\n'
                       'neu_sfs: {sfs1}\n'
                       'neu_theta_range: {t1}, {t2}\n'
                       'neu_e_range: {e1}, {e2}\n'
                       'sel_m: {m3}\n'
                       'sel_sfs: {sfs3}\n')

        model_ctrls = model_ctrls.format(n=n,
                                         r1=r_r[0], r2=r_r[1],
                                         m1=sfs_m['neutral_SNP'][1],
                                         fold=str(snp_fold).lower(),
                                         sfs1=', '.join([str(x) for x in sfs_m['neutral_SNP'][0]]),
                                         m3=sfs_m['selected_SNP'][1],
                                         sfs3=', '.join([str(x) for x in sfs_m['selected_SNP'][0]]),
                                         t1=theta_r[0], t2=theta_r[1],
                                         e1=error_r[0], e2=error_r[1])

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

            if self.dfe_optional_opts == '':
                self.set_dfe_optional_opts()

        self.dfe_opts = dfe_param


class ResultsFile(object):

    def __init__(self, anavar_results_file):

        """
        Creates an anavar results file object
        :param anavar_results_file: file
        """

        contents = anavar_results_file.readlines()
        self.head = contents[0:7]
        self.control_location = self.head[0].split()[2]
        self.param = tuple(self.head[2].split()[4:])
        self.columns = tuple(self.head[6].split())

        self.data = contents[7:]

    def seed(self):

        """
        returns the seed anavar was run with
        :return: int
        """

        seed = int(self.head[4].rstrip().replace('seed = ', ''))

        return seed

    def free_parameters(self):

        """
        returns all free parameters in results file
        :return: tuple
        """

        return self.param

    def control_file(self):

        """
        returns control file path
        :return: str
        """

        return self.control_location

    def header(self):

        """
        returns the column headers line from the file
        :return: tuple
        """

        return self.columns

    def __results_to_dict(self, line):

        """
        converts anavar results line to dict
        :param line: str
        :return: dict
        """

        line = line.split()
        line_dict = {}

        for i in range(0, len(self.columns)):

            col_val = line[i]
            if '.' in col_val or 'e' in col_val:
                col_val = float(col_val)
            else:
                col_val = int(col_val)

            line_dict[self.columns[i]] = col_val

        return line_dict

    def estimates(self, as_string=False):

        """
        returns generator of all estimates in dictionary form
        :param as_string: bool
        :return: generator
        """

        for line in self.data:
            if as_string:
                line = line.rstrip()
            else:
                line = self.__results_to_dict(line)
            yield line

    def ml_estimate(self, as_string=False):

        """
        returns the maximum likelihood estimate
        :param as_string: bool
        :return: dict
        """

        ml = self.data[0]

        if as_string:
            return ml.rstrip()
        else:
            return self.__results_to_dict(ml)

    def converged(self):

        """
        assesses convergence for results file
        :return: bool
        """

        top_lnls = [abs(x['lnL']) for x in self.estimates()][0:2]

        # dif = (top_lnls[0] - top_lnls[1]) / (0.5 * (top_lnls[0] + top_lnls[1]))
        dif = top_lnls[0] - top_lnls[1]

        if abs(dif) < 0.1:
            return True

        else:
            return False

    def bounds_hit(self, theta_r=(1e-6, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5),
                   shape_r=(1e-3, 200), scale_r=(0.1, 1e3), r_r=(0.05, 5.0)):

        """
        lists parameters too close to estimate boundaries
        :param theta_r: tuple
        :param gamma_r: tuple
        :param error_r: tuple
        :param shape_r: tuple
        :param scale_r: tuple
        :param r_r: tuple
        :return: list
        """

        params_hitting_limits = []
        skip_params = ['run', 'imp', 'exit_code', 'lnL']

        ml_dict = self.ml_estimate()
        for param in ml_dict.keys():

            if re.search(r'theta', param):
                bounds = theta_r

            elif re.search(r'gamma', param):
                bounds = gamma_r

            elif re.search(r'scale', param):
                bounds = scale_r

            elif re.search(r'shape', param):
                bounds = shape_r

            elif re.search(r'r_\d{1,2}', param):
                bounds = r_r

            elif param in skip_params:
                continue

            else:
                bounds = error_r

            try:
                b_hit = self.__boundary_hit(ml_dict[param], bounds[1], bounds[0], 1e-5)

            # catches error = 0
            except ZeroDivisionError:
                params_hitting_limits.append(param)
                continue

            if b_hit:
                params_hitting_limits.append(param)

        return params_hitting_limits

    def __boundary_hit(self, mle, upr_b, lwr_b, min_dif):

        """
        decides if ml estimate is too close to boundary
        :param mle: float
        :param upr_b: float
        :param lwr_b: float
        :param min_dif: float
        :return: bool
        """

        for b in [abs(upr_b), abs(lwr_b)]:

            dif = (abs(mle) - b) / (0.5 * (abs(mle) + b))

            if abs(dif) < abs(min_dif):
                return True

        return False

    def data_type(self):

        """
        assigns data as either indel, snp or snp_indel
        :return: str
        """

        param_str = ','.join(self.free_parameters())

        if 'ins' in param_str and 'del' in param_str:
            if 'snp' in param_str:
                return 'snp_indel'
            else:
                return 'indel'
        else:
            return 'snp'

    def num_class(self):

        """
        returns number of classes in results file
        :return: int
        """
        params = self.free_parameters()

        thetas = set([x.split('_')[-1] for x in params if 'theta' in x])

        if 'theta' in thetas:
            thetas.remove('theta')

        if 1 >= len(thetas) >= 0:
            return 1
        else:
            return len(thetas)

    def dfe(self):

        """
        returns type of distribution
        :return: str
        """
        param_str = ','.join(self.free_parameters())

        if 'scale' in param_str and 'shape' in param_str:
            return 'continuous'
        else:
            return 'discrete'

    def num_runs(self):

        """
        gives the number of runs listed in the results file
        :return: int
        """

        return len(self.data)

    def get_alpha(self, dn, ds, var_type):

        """
        calculates alpha as described in equation 19 Barton and Zeng 2018
        :param dn: float
        :param ds: float
        :param var_type: str
        :return: str
        """

        if var_type not in ['snp', 'indel', 'ins', 'del']:
            raise ValueError('Unrecognised variant type')

        # if self.data_type() != 'indel':
        #     raise NotImplementedError('Only implemented for INDELs')

        mle = self.ml_estimate()

        theta_keys = [x for x in mle.keys() if 'theta' in x and 'sel' in x]

        if self.dfe() == 'continuous':

            gamma_keys = [x for x in mle.keys() if 'scale' in x and 'sel' in x] + \
                         [x for x in mle.keys() if 'shape' in x and 'sel' in x]

        else:

            gamma_keys = [x for x in mle.keys() if 'gamma' in x and 'sel' in x]

        # separate by class

        c = self.num_class()

        if var_type == 'indel':
            variants_to_process = ['ins', 'del']
        elif var_type == 'ins':
            variants_to_process = ['ins']
        elif var_type == 'del':
            variants_to_process = ['del']
        else:
            variants_to_process = ['']

        ts = []
        gs = []

        for var in variants_to_process:

            for i in range(1, c+1):

                if self.dfe() == 'continuous':

                    ti_key = [x for x in theta_keys if var in x]

                    dfei_keys = [x for x in gamma_keys if var in x]

                    gi = -(mle[dfei_keys[0]] * mle[dfei_keys[1]])

                else:

                    ti_key = [x for x in theta_keys if var in x and str(i) in x]

                    gi_key = [x for x in gamma_keys if var in x and str(i) in x]

                    gi = mle[gi_key[0]]

                ti = mle[ti_key[0]]

                ts.append(ti)
                gs.append(gi)

        var_alpha = alpha(ts, gs, dn, ds)
        return var_alpha


class MultiResultsFile(ResultsFile):

    def __init__(self, file_list):

        """
        takes a list of anavar results files and merges them into one,
        then calls ResultsFile instance on merged file before deleting
        merged file
        :param file_list: list
        """

        merged_name = file_list[0] + '.merged.results.txt'

        merge_results(file_list, merged_name)

        ResultsFile.__init__(self, open(merged_name))

        os.remove(merged_name)


def merge_results(file_list, name):

    """
    takes list of results files and merges them
    :param file_list: list
    :param name: str
    :return:
    """

    all_estimates = []
    merged_file = open(name, 'w')

    first = True
    for x in file_list:
        with open(x) as results:

            res = ResultsFile(results)

            if first:
                print(''.join(res.head).rstrip(), file=merged_file)
                first = False

            estimates = res.estimates(as_string=True)
            lnls = [z['lnL'] for z in res.estimates()]

            est_by_lnl = zip(lnls, estimates)

            all_estimates += est_by_lnl

    sorted_data = sorted(all_estimates, reverse=True)

    for mle in sorted_data:
        print(mle[1], file=merged_file)

    merged_file.close()
