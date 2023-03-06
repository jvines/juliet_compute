"""Error raising and handling."""
import sys

from termcolor import colored


class Error(Exception):
    """Base class for exceptions in this module."""

    def __init__(self, *args):
        if args:
            self.errorname = args[0]
            self.message = args[1]
        else:
            self.message = None

    def __repr__(self):
        """Error identification for logging."""
        return self.errorname

    def __str__(self):
        """Error identification for logging."""
        return self.errorname

    def __raise__(self):
        """Raise an exception and print the error message."""
        self.warn(color='red')
        sys.exit()

    def warn(self, color='yellow'):
        """Print error message."""
        print(colored('An exception was caught!', color), end=': ')
        print(colored(self, color), end='\nError message: ')
        print(colored(self.message, color))


class UnitError(Error):
    def __init__(self, name):
        self.errorname = f'{name} Error'
        self.message = 'Unknown unit selected. Allowed ones are '
        self.message += '`jupiter` or `earth`.'


class StellarParamsError(Error):
    params = {
        'rs': 'Stellar Radius',
        'ms': 'Stellar Mass',
        'teff': 'Effective Temperature'
    }

    def __init__(self, param):
        self.errorname = 'Stellar Param Error'
        self.message = f'Parameter `{self.params[param]}` is missing.'
        self.message += ' Some calculations will be unavailable.'