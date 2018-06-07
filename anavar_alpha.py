from __future__ import print_function
import math


def mean_fix_prob(mle_thetas, mle_gammas):

    """
    equation 18 Barton and Zeng 2018
    :param mle_thetas: tuple
    :param mle_gammas: tuple
    :return: float
    """

    theta = sum(mle_thetas)

    numerator = 0.0

    for c in range(0, len(mle_thetas)):

        # Kai email 02/05/18: For eq. 18, you need to set gamma for beneficial site classes to something like
        # -10000, so that they contribute to de novo mutations, but not the expected divergence level due to
        # deleterious mutations.
        if mle_gammas[c] > 0:
            gamma = -10000
        else:
            gamma = mle_gammas[c]

        if gamma < -500:
            class_calc = 0
        else:
            class_calc = (mle_thetas[c] * gamma) / (1 - math.exp(-gamma))

        numerator += class_calc

    mu_mean = 1.0/theta * numerator

    return mu_mean


def alpha(mle_thetas, mle_gammas, dn, ds):

    """
    equation 19 Barton and Zeng 2018
    :param mle_thetas: tuple
    :param mle_gammas: tuple
    :param dn: int
    :param ds: int
    :return: float
    """

    mu_mean = mean_fix_prob(mle_thetas=mle_thetas, mle_gammas=mle_gammas)

    a = (float(dn) - (float(ds) * mu_mean)) / float(dn)

    return a


def debug():

    mle_ts = (1.8e-5, 7.2e-4, 5.3e-5, 0.0011)
    mle_gs = (1.98, -1566.4, -1.69, -642.5)  # -10000 actually 1.98

    mle_ts2 = (1.6e-5, 1.9e-4, 4.9e-5, 0.001)
    mle_gs2 = (-1.31, -284.1, -3.77, -454.8)

    dn = 0.000186995113194375
    ds = 0.00287748490295213

    ds_2 = 0.0055

    print('# table 5 Barton and Zeng 2018 results')
    print('alpha with non-coding reference:', alpha(mle_ts, mle_gs, dn, ds))
    print('alpha with 4fold reference:', alpha(mle_ts2, mle_gs2, dn, ds_2))
