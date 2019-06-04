import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from scipy.signal import convolve, convolve2d, correlate2d
from scipy.optimize import curve_fit
from astropy.table import Table, Column
import time
import pdb

class Grid:
    """
    N-dimensional grid of two-state spins (1 = up, -1 = down).

    Parameters
    ----------
    Each argument is its own parameter

    shape : int or tuple of ints
        if int (>0), initialize spin chain with length shape
        if tuple (of ints all > 0), initialize

    init : float or str
        if float, probability that a given spin is initialized in up state. Number-like
             types will be cast to float
        if str, will attempt to look for a preset grid initialization
    """
    def __init__(self, shape, init=0.5):
        # Check user input
        shape_type_error = ValueError('shape must be either a positive integer or tuple of positive integers.')
        cast_warning = False
        if type(shape) == int:
            if int <= 0:
                raise shape_type_error
            shape = (shape)
        elif (type(shape) == list) or (type(shape) == np.ndarray):
            # Check that all items in list are positive integers
            for item in shape:
                if type(item) != int:
                    try:
                        item = int(item)
                        if cast_warning == False:
                            print('Warning: One or more dimensions in shape has been cast to int.')
                    except:
                        raise shape_type_error
                if int <= 0:
                    raise shape_type_error
            shape = tuple(shape)

        presets = []
        pup_given = True # remains true if user directly gives a probability for init
        init_type_error = ValueError('init must be either a float in [0, 1] or recognized str preset.')
        if type(init) == str:
            if init not in presets:
                raise init_type_error
            pup_given = False
        else:
            try:
                init = float(init)
            except:
                raise init_type_error
            if (init > 1) or (init < 0):
                raise init_type_error

        # Initialize grid
        if pup_given:
            self.grid = np.random.choice([1, -1], shape, p=[init, 1-init])
        else:
            pass # Populate this section when the first preset is added

        assert np.shape(self.grid) == shape

class Interaction:
    """
    Object representing terms in the Hamiltonian of the Ising model.

    Parameters
    ----------
    grid : Grid
        Ising grid described by interactions
    strength : float
        Coupling strength
    couplings : array-like (stored as list) with tuple elements
        First N dimensions of couplings must be the same as grid. Last dimension
        has length equal to number of spins coupled together, i.e., K such that
        the interaction enters into H as J*s1*s2*...*sK
    """
    def __init__(self, grid, strength, couplings):
        ## TODO: Write some user checks to make sure they enter couplings correctly

        self.strength = strength
        self.couplings = couplings

class FieldCoupling(Interaction):
    """
    Term in the Hamiltonian which enters as h*Si for specified spins Si.

    Parameters
    ----------
    grid : Grid
        Ising grid described by interactions
    strength : float
        Coupling strength
    preset : str
        Setting describing which spins couple to the field. Valid presets:
          'all', all spins couple to field (default)
          'probabilistic', spins randomly couple to field with probability randomp
    randomp : float in [0,1], inclusive
        if preset == 'probabilistic', couple spins to field with this probability.
        otherwise, this parameter is ignored.
    """
    def __init__(self, strength, preset='all', randomp=None):
        # Check user input
        strength_type_error = ValueError('strength must be a float.')
        if type(strength) != float:
            try:
                strength = float(strength)
            except:
                raise strength_type_error

        preset_type_error = ValueError('preset must be a str.')
        preset_sel_error = ValueError('A valid preset was not given.')
        presets = ['all', 'probabilistic']
        if type(preset) != str:
            raise preset_type_error
        elif preset not in presets:
            raise preset_sel_error

        # Warn user if randomp is ignored
        if (preset != 'probabilistic') and (randomp != None):
            print('Warning: randomp parameter was ignored.')

        # Build coupling list
        shape = list(grid.grid.shape)
        Ndims = len(shape)

        # Populate couplings list
        if preset == 'all':
            # Create array where last index is coordinate number, second to last is point number (there is only one),
            # and other dimensions are just grid dimensions
            couplings = np.indices(shape)
            couplings = np.array([couplings])

            couplings = np.rollaxis(couplings, 0, start=len(coupling.shape))
            couplings = np.rollaxis(couplings, 0, start=len(coupling.shape))

            self.couplings = couplings

            # Fuck this format can't take into account that one interaction can optionally couple fields
            # more generally, it can't handle if different points are involved in different numbers of couplings fuck















        # Define strength
        self.strength = strength



class DipoleCoupling(Interaction):
    """
    Term in the Hamiltonian which enters as J*Si*Sj for specified pairs of Si, Sj.

    Parameters
    ----------
    grid : Grid
        Ising grid described by interactions
    strength : float
        Coupling strength
    preset : str
        Setting describing which spins are coupled to each other. Valid presets:
          'prismnn', spins are coupled to their neighbors
          'prismtri', spins are coupled to six neighbors in triangular fashion in first two dimensions,
            nn interactions in the others
          'prismhex', spins are coupled to three neighbors in hexagonal fashion,
            nn interactions in the others
          'glassnn', randomly initialize prismnn bonds with probability glassp
          'glasstri', randomly initialize prismnn bonds with probability glassp
          'glasshex', randomly initialize prismnn bonds with probability glassp
    glassp : float in [0,1], inclusive
        if the preset is a glass, this describes the bond probability.
        Otherwise, this parameter is ignored.
    """
    ...








    # Define strength
    self.strength = strength

































#
