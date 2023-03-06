"""Main driver."""
import argparse

from .utils import create_dir, param_array
from .juliet_utils import *
from .error import StellarParamsError


def arg_parse():
    """Parse command line arguments."""
    desc = """
    Compute planetary and orbital parameters from Juliet posteriors.
    """
    p = argparse.ArgumentParser('Juliet Extractor', description=desc)
    p.add_argument('input', help='The juliet posterior pkl file.', type=str)
    p.add_argument('output',
                   help='The output location for the computed parameters.',
                   type=str)
    p.add_argument('-rs', help='The stellar radius.', type=float,
                   required=False)
    p.add_argument('-rse', help='The uncertainty in stellar radius.',
                   type=float, default=0, reqiured=False)
    p.add_argument('-ms', help='The stellar mass.', type=float, required=False)
    p.add_argument('-mse', help='The uncertainty in stellar mass.', type=float,
                   default=0, required=False)
    p.add_argument('-teff', help='The effective temperature.', type=float,
                   required=False)
    p.add_argument('-teffe', help='The uncertainty in effective temperature.',
                   type=float, default=0, required=False)
    p.add_argument('-verbose', help='Set verbosity.', required=False, type=bool)
    p.add_argument('-munit', help='Units for planetary mass.', required=False,
                   type=str, default='jupiter')
    p.add_argument('-runit', help='Units for planetary radius.', required=False,
                   type=str, default='jupiter')
    alb = 'Bond albedo for each planet. If a single number is provided then'
    alb += ' one albedo will be used for all. Structure is -Ab p1 0 p2 0.5'
    p.add_argument('-Ab', help='Bond albedo.', required=False, type=float,
                   nargs='+')
    fixed = 'Fixed parameters. Must be param value param value.'
    fixed += ' For example -fixed P_p1 5 ecc_p2 0 omega_p2 90'
    p.add_argument('-fixed', nargs='+', help=fixed, required=False)

    return p.parse_args()


def main():
    args = arg_parse()  # Read command line args
    verbose = args.verbose
    create_dir(args.output, verbose=verbose)  # Create output directory
    # Compile fixed parameters...
    fixed = args.fixed
    fixed_dict = dict()
    if fixed is not None:
        for par, val in zip(fixed[::2], fixed[1::2]):
            fixed_dict[par] = val
    # Get stellar params.
    rs, ers = args.rs, args.rse
    if rs is None:
        StellarParamsError('rs').warn()
    ms, ems = args.ms, args.mse
    if ms is None:
        StellarParamsError('ms').warn()
    teff, eteff = args.teff, args.teffe
    if teff is None:
        StellarParamsError('teff').warn()
    # Read Juliet posteriors.
    post = read_juliet_posteriors(args.input)
    # Fill post with fixed params to make my life easier.
    nsamp = len(post[list(post.keys())[0]])  # Get number of samples.
    for par in fixed_dict.keys():
        post[par] = np.ones(nsamp) * fixed_dict[par]
    # Fill post with stellar params to make my life easier.
    # -1 means no param provided.
    post['Ms'] = param_array(ms, ems, nsamp)
    post['Rs'] = param_array(rs, ers, nsamp)
    post['Teff'] = param_array(teff, eteff, nsamp)
    nplans = get_number_of_planets(post)
