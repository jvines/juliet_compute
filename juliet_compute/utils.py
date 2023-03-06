"""Other utilities."""
import os

import numpy as np
from scipy.special import erf
from scipy.stats import gaussian_kde, norm
from termcolor import colored


def estimate_pdf(distribution):
    """Estimates the PDF of a distribution using a gaussian KDE.

    Parameters
    ----------
    distribution: array_like
        The distribution.
    Returns
    -------
    xx: array_like
        The x values of the PDF.
    pdf: array_like
        The estimated PDF.
    """
    kde = gaussian_kde(distribution)
    xmin, xmax = distribution.min(), distribution.max()
    xx = np.linspace(xmin, xmax, 500)
    pdf = kde(xx)
    return xx, pdf


def estimate_cdf(distribution, hdr=False):
    """Estimate the CDF of a distribution."""
    h, hx = np.histogram(distribution, density=True, bins=499)
    cdf = np.zeros(500)  # ensure the first value of the CDF is 0
    if hdr:
        idx = np.argsort(h)[::-1]
        cdf[1:] = np.cumsum(h[idx]) * np.diff(hx)
    else:
        cdf[1:] = np.cumsum(h) * np.diff(hx)
    return cdf


def credibility_interval_hdr(xx, pdf, cdf, sigma=1.):
    """Calculate the highest density region for an empirical distribution.

    Reference: Hyndman, Rob J. 1996

    The HDR is capable of calculating more robust credible regions
    for multimodal distributions. It is identical to the usual probability
    regions of symmetric about the mean distributions. Using this then should
    lead to more realistic errorbars and 3-sigma intervals for multimodal
    outputs.

    Parameters
    ----------
    xx: array_like
        The x values of the PDF (and the y values of the CDF).
    pdf: array_like
        The PDF of the distribution.
    cdf: array_like
        The CDF of the distribution.
    sigma: float
        The confidence level in sigma notation. (e.g. 1 sigma = 68%)

    Returns
    -------
    best: float
        The value corresponding to the peak of the posterior distribution.
    low: float
        The minimum value of the HDR.
    high: float
        The maximum value of the HDR.
    """
    # Get best fit value
    best = xx[np.argmax(pdf)]
    z = erf(sigma / np.sqrt(2))
    # Sort the pdf in reverse order
    idx = np.argsort(pdf)[::-1]
    # Find where the CDF reaches 100*z%
    idx_hdr = np.where(cdf >= z)[0][0]
    # Isolate the HDR
    hdr = pdf[idx][:idx_hdr]
    # Get the minimum density
    hdr_min = hdr.min()
    # Get CI
    low = xx[pdf > hdr_min].min()
    high = xx[pdf > hdr_min].max()
    return best, low, high


def credibility_interval(post, sigma=1., axis=None):
    """Calculate bayesian credibility interval.

    Parameters:
    -----------
    post : array_like
        The posterior sample over which to calculate the bayesian credibility
        interval.
    sigma : float, optional
        Confidence level.
    Returns:
    --------
    med : float
        Median of the posterior.
    low : float
        Lower part of the credibility interval.
    up : float
        Upper part of the credibility interval.

    """
    z = erf(sigma / np.sqrt(2))

    lower_percentile = 100 * (1 - z) / 2
    upper_percentile = 100 * (1 + z) / 2
    low, med, up = np.percentile(
        post, [lower_percentile, 50, upper_percentile], axis=axis
    )
    return med, low, up


def create_dir(path, verbose=False):
    """Create a directory."""
    try:
        os.mkdir(path)
    except OSError:
        err_msg = f"Creation of the directory {path:s} failed. "
        err_msg += "It might already exist"
        if verbose:
            print(colored(err_msg, 'red'))
        pass
    else:
        if verbose:
            print(colored(f"Created the directory {path:s}", 'blue'))
        pass
    pass


def montecarlo(loc, scale, size):
    """Perform a gaussian montecarlo simulation with loc and scale."""
    if scale == 0:
        return np.ones(size) * loc
    return norm.rvs(loc=loc, scale=scale, size=size)


def param_array(loc, scale, size):
    """Generate a parameter array with a loc scale and size."""
    if loc is not None and scale > 0:
        arr = montecarlo(loc, scale, size)
    elif loc is not None and scale == 0:
        arr = np.ones(size) * loc
    else:  # loc is None
        arr = np.ones(size) * -1
    return arr
