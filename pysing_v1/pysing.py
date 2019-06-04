import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve, convolve2d, correlate2d
from scipy.optimize import curve_fit
import matplotlib.animation as ani
from astropy.table import Table, Column
import time
import pdb

method = 'metropolis'

def run_ising(h, J, T, N=30, f=0.1, steps=100, plot=True, grid='square', param1=0):
    """ Hamiltonian of the form sum_i -h*Sz(i) - J*Sz(i)Sz(i+1) following canonical ensemble with temperature T (k=1)
    """
    start = time.time()

    # Initialize board with dimension NxN
    X = [np.arange(0,N) + 1 - int(np.ceil(N/2))]
    X = np.array(N * X)
    Y = X.T

    board = np.random.choice([-1,1], (N,N))

    if grid == 'hexagon': # If hexagon, get rid of extra spaces
        board[np.where((X + Y) % 3 == 0)] = 0
    if grid == 'glass8': # If glass8, draw random bonds (param1 in this case is bonding probability)
        bond_ru = np.random.choice([False, True], (N,N), p=[1-param1, param1])
        bond_r0 = np.random.choice([False, True], (N,N), p=[1-param1, param1])
        bond_rd = np.random.choice([False, True], (N,N), p=[1-param1, param1])
        bond_0d = np.random.choice([False, True], (N,N), p=[1-param1, param1])

    # Initialize figure
    if plot:
        plt.close()
        fig, axes = plt.subplots()

        im = axes.imshow(board, cmap='Wistia', origin='lower')

    # Define update function by using Boltzmann factors to calculate the probability of being in a given state
    def update(input):
        """ If plot==True, input is just a dummy variable. If plot==False, input=the current state of the board. """
        # Pull board
        if plot:
            board = im.get_array()
        else:
            board = input

        # Update
        if method == 'metropolis':
            for ii in range(int(np.round(f*N**2))):
                x = np.random.randint(N)
                y = np.random.randint(N)
                dE = 2*h*board[x,y]

                if grid == 'square':
                    dE += 2*J*board[x,y]*board[x, (y-1)%N]
                    dE += 2*J*board[x,y]*board[x, (y+1)%N]
                    dE += 2*J*board[x,y]*board[(x-1)%N, y]
                    dE += 2*J*board[x,y]*board[(x+1)%N, y]
                elif grid == 'triangle': # For some reason, there's a bug that breaks the X, Y symmetry in this simulation
                    dE += 2*J*board[x,y]*board[x, (y-1)%N]
                    dE += 2*J*board[x,y]*board[x, (y+1)%N]
                    dE += 2*J*board[x,y]*board[(x-1)%N, y]
                    dE += 2*J*board[x,y]*board[(x+1)%N, y]
                    dE += 2*J*board[x,y]*board[(x+1)%N, (y+1)%N]
                    dE += 2*J*board[x,y]*board[(x-1)%N, (x-1)%N]
                elif grid == 'hexagon': # This is wrong
                    if (x + y - 1) % 3 == 0:
                        dE += 2*J*board[x,y]*board[x, (y+1)%N]
                        dE += 2*J*board[x,y]*board[(x+1)%N, y]
                        dE += 2*J*board[x,y]*board[(x-1)%N, (x-1)%N]
                    elif (x + y + 1) % 3 == 0:
                        dE += 2*J*board[x,y]*board[x, (y-1)%N]
                        dE += 2*J*board[x,y]*board[(x-1)%N, y]
                        dE += 2*J*board[x,y]*board[(x+1)%N, (x+1)%N]
                    else:
                        dE = -999 # Flip it, who cares (it's ==0)
                elif grid == 'glass8':
                    if bond_ru[x, y]:
                        dE += 2*J*board[x,y]*board[(x+1)%N, (y+1)%N]
                    if bond_r0[x, y]:
                        dE += 2*J*board[x,y]*board[(x+1)%N, y]
                    if bond_rd[x, y]:
                        dE += 2*J*board[x,y]*board[(x+1)%N, (y-1)%N]
                    if bond_0d[x, y]:
                        dE += 2*J*board[x,y]*board[x, (y-1)%N]

                    if bond_ru[(x-1)%N, (y-1)%N]:
                        dE += 2*J*board[x,y]*board[(x-1)%N, (y-1)%N]
                    if bond_r0[(x-1)%N, y]:
                        dE += 2*J*board[x,y]*board[(x-1)%N, y]
                    if bond_rd[(x-1)%N, (y+1)%N]:
                        dE += 2*J*board[x,y]*board[(x-1)%N, (y+1)%N]
                    if bond_0d[x, (y+1)%N]:
                        dE += 2*J*board[x,y]*board[x, (y+1)%N]
                else:
                    raise ValueError('Invalid grid.')

                if dE < 0:
                    board[x,y] *= -1
                elif T != 0:
                    if np.random.uniform(0, 1) < np.exp(-dE / T):
                        board[x,y] *= -1
        else:
            raise ValueError('No such method.')

        if plot:
            im.set_array(board)
            return im
        else:
            return board

    if plot:
        #anim = ani.FuncAnimation(fig, update, interval=1)
        anim = ani.FuncAnimation(fig, update, interval=1, frames=300)
        anim.save('fm1.gif', writer='imagemagick', fps=30)
        print('Time elapsed: {0} s'.format(time.time() - start))
        plt.show()
    else:
        for ii in range(steps):
            board = update(board)

        print('Time elapsed: {0} s'.format(time.time() - start))
        return board

def calculate_correlation_function(board, stdev=0):
    """ Function which takes a board and returns correlation function.
    Assumes board is homogenous and isotropic so that correlation function is purely a function of correlation length. """
    correlation = correlate2d(board, board, mode='same', boundary='wrap')
    X = [np.arange(0,len(board)) + 1 - int(np.ceil(len(board)/2))]
    X = np.array(len(board) * X)
    Y = X.T
    dist = np.hypot(X, Y)

    plt.imshow(correlation)

    dist = dist.flatten()
    correlation = correlation.flatten()

    # Convert 2d correlation array into correlation function evaluated at every available point

    s_arr = np.unique(dist)
    s_arr.sort()

    c_arr = np.array([])

    for s in s_arr:
        c_arr = np.append(c_arr, np.mean(correlation[dist == s]))

    # Apply a smoothing kernel with standard deviation stdev
    if stdev != 0:
        x = np.arange(-3 * np.ceil(stdev), 3 * np.ceil(stdev) + 1)
        ker = np.exp(-0.5 * x ** 2 / stdev ** 2) / np.sqrt(2 * np.pi * stdev ** 2)
        c_arr = convolve(c_arr, ker, mode='same')

    # Normalize and chop off all values past board length / 2 (since the BC's wrap)
    c_arr /= c_arr[0]
    c_arr = c_arr[s_arr < len(board) / 2]
    s_arr = s_arr[s_arr < len(board) / 2]

    return s_arr, c_arr

# Calculate order parameters
def order_parameters(board):
    """
    Given a board, calculates the following two order parameters:
    * M: Sum of all ups minus sum of all downs divided by Nspins (magnetization, will be nonzero in FM case)
    * m: Average absolute value of local magnetization
    * q: Same as local magnetization but weighted by position parity (+1/-1), absolute valued (will be nonzero in AFM case)
    """
    # Get M
    M = np.sum(board) / len(board) ** 2

    # Get m
    fm_ker = np.array([[+1, +1, +1],
                       [+1, +1, +1],
                       [+1, +1, +1]]) / 9
    m_board = convolve2d(board, fm_ker, boundary='wrap', mode='same')
    m = np.mean(np.abs(m_board))

    afm_ker = np.array([[+1, -1, +1],
                        [-1, +1, -1],
                        [+1, -1, +1]]) / 9
    q_board = convolve2d(board, afm_ker, boundary='wrap', mode='same')
    q = np.mean(np.abs(q_board))

    return M, m, q

def ising_grid(grid='square', param1=0):
    """ Run a 30x30x30 (h, J, T) grid on a 30x30 grid and extract M, m, and q for each.

    ~~ 5/30 ~~ 1559216678.0852456
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 30, 30)
    time: 2536.3413004875183 s # T = 0 was wrong

    ~~ 5/30 ~~ 1559217117.5116973
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = [0]

    ~~ 5/30 ~~ 1559257188.867237
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    time: 2499.3187720775604 s

    ~~ 5/30 ~~ 1559269113.4033473
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    time: 5595.280873298645 s # Attempted hexagonal implementation which was wrong

    ~~ 5/30 ~~ 1559269814.3782125 (triangle)
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    time: 6790.1584849357605 s # There seems to be a problem here with preferentially horizontal domains, don't know why

    ~~ 5/30 ~~ 1559273938.5639591 (glass8, p=0.5)
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    time: 4393.840525865555 s

    ~~ 5/31 ~~ 1559288446.0055742 (hexagon)
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    time: 2112.42906665802 s

    The last two are combined into a correct file 'pysing_run_5-30f.fits', but this wasn't good
    So I just sliced the T=0 case off of the first one and called it 'pysing_run_5-30c.fits'

    """
    #h_grid = np.linspace(-2, 2, 30)
    #J_grid = np.linspace(-2, 2, 30)
    #T_grid = np.linspace(0, 30, 30)
    h_grid = np.linspace(-2, 2, 30)
    J_grid = np.linspace(-2, 2, 30)
    T_grid = np.linspace(0, 10, 31)
    T_grid = T_grid[1:]
    #T_grid = np.linspace(0, 30, 30)

    start = time.time()

    h_arr = np.array([])
    J_arr = np.array([])
    T_arr = np.array([])
    M_arr = np.array([])
    m_arr = np.array([])
    q_arr = np.array([])

    for h in h_grid:
        for J in J_grid:
            for T in T_grid:
                board = run_ising(h, J, T, N=30, f=0.1, steps=100, plot=False, grid=grid, param1=param1)
                print('========')
                print('h = {0}'.format(h))
                print('J = {0}'.format(J))
                print('T = {0}'.format(T))

                M, m, q = order_parameters(board)

                h_arr = np.append(h_arr, h)
                J_arr = np.append(J_arr, J)
                T_arr = np.append(T_arr, T)
                M_arr = np.append(M_arr, M)
                m_arr = np.append(m_arr, m)
                q_arr = np.append(q_arr, q)

    h_col = Column(h_arr, name='h')
    J_col = Column(J_arr, name='J')
    T_col = Column(T_arr, name='T')
    M_col = Column(M_arr, name='M')
    m_col = Column(m_arr, name='m')
    q_col = Column(q_arr, name='q')

    tab = Table([h_col, J_col, T_col, M_col, m_col, q_col])
    Table.write(tab, 'pysing-grid_{0}.fits'.format(time.time()))

    print('Grid time elapsed: {0} s'.format(time.time() - start))

def build_phase_diagram(fname, fix='h'):
    """ Function to read ising_grid() output and make phase diagrams """
    data = Table.read(fname)
    cmap1 = 'seismic'
    cmap2 = 'Oranges'
    cmap3 = 'Greens'

    if fix == 'h':
        col1 = 'T'
        col2 = 'J'
        aspect = 'auto'
    elif fix == 'J':
        col1 = 'h'
        col2 = 'T'
        aspect = 'auto'
    elif fix == 'T':
        col1 = 'h'
        col2 = 'J'
        aspect = None

    s_arr = np.unique(data[fix])
    s_arr.sort()

    # Figure out all unique values for unfixed parameters
    col1_vals = np.unique(data[col1])
    col1_vals.sort()
    col2_vals = np.unique(data[col2])
    col2_vals.sort()

    # Figure out how to pad edges of imshow plot
    col1_halfspace = (col1_vals[1] - col1_vals[0]) / 2
    col2_halfspace = (col2_vals[1] - col2_vals[0]) / 2
    extent = [np.min(col1_vals)-col1_halfspace,
              np.max(col1_vals)+col1_halfspace,
              np.min(col2_vals)-col2_halfspace,
              np.max(col2_vals)+col2_halfspace]

    M_grid_list = []
    m_grid_list = []
    q_grid_list = []

    for s in s_arr:
        s_data = data[data[fix] == s]

        # Create parameter grids (assume equally spaced square grid same across all fixed values)
        M_grid = np.zeros(shape=(len(col1_vals), len(col2_vals)))
        m_grid = np.zeros(shape=(len(col1_vals), len(col2_vals)))
        q_grid = np.zeros(shape=(len(col1_vals), len(col2_vals)))

        for ii in range(len(col1_vals)):
            for jj in range(len(col2_vals)):
                idx = np.where( (s_data[col1] == col1_vals[ii]) & (s_data[col2] == col2_vals[jj]) )
                M_grid[jj,ii] = s_data['M'][idx]
                m_grid[jj,ii] = s_data['m'][idx]
                q_grid[jj,ii] = s_data['q'][idx]

        M_grid_list.append(M_grid)
        m_grid_list.append(m_grid)
        q_grid_list.append(q_grid)

    # Construct plot
    plt.close()
    fig = plt.figure(figsize=(10,3))

    plt.subplot(131)
    M_im = plt.imshow(M_grid_list[0], vmin=-1, vmax=1, extent=extent, cmap=cmap1, origin='lower', aspect=aspect)
    plt.xlabel(col1)
    plt.ylabel(col2)
    cbar = plt.colorbar()
    cbar.set_label('M')

    plt.subplot(132)
    m_im = plt.imshow(m_grid_list[0], vmin=0, vmax=1, extent=extent, cmap=cmap2, origin='lower', aspect=aspect)
    plt.xlabel(col1)
    plt.ylabel(col2)
    cbar = plt.colorbar()
    cbar.set_label('m')

    plt.subplot(133)
    q_im = plt.imshow(q_grid_list[0], vmin=0, vmax=1, extent=extent, cmap=cmap3, origin='lower', aspect=aspect)
    plt.xlabel(col1)
    plt.ylabel(col2)
    cbar = plt.colorbar()
    cbar.set_label('q')

    plt.tight_layout()

    def update(kk):
        M_im.set_data(M_grid_list[kk])
        m_im.set_data(m_grid_list[kk])
        q_im.set_data(q_grid_list[kk])

        # Change titles
        plt.subplot(131)
        plt.title('{0} = {1}'.format(fix, np.round(s_arr[kk], 2)))
        plt.subplot(132)
        plt.title('{0} = {1}'.format(fix, np.round(s_arr[kk], 2)))
        plt.subplot(133)
        plt.title('{0} = {1}'.format(fix, np.round(s_arr[kk], 2)))

        return M_im, m_im, q_im

    # Make the animation
    anim = ani.FuncAnimation(fig, update, np.arange(len(s_arr)),
                            interval=100, blit=True)

    anim.save('phase_diagram_fixed-{0}.gif'.format(fix), writer='imagemagick', fps=10)



#def extract_correlation_parameters(s_arr, c_arr): # This function is broken
#    """ Fit a function of the form exp(-gamma *  s) * cos(pi * s) """
#    # First take the smoothed absolute value of the data and fit a straight exponential
#    c_abs = np.abs(c_arr)
#
#    stdev = 3
#    x = np.arange(-3 * np.ceil(stdev), 3 * np.ceil(stdev) + 1)
#    ker = np.exp(-0.5 * x ** 2 / stdev ** 2) / np.sqrt(2 * np.pi * stdev ** 2)
#    c_abs = convolve(c_abs, ker, mode='same')
#
#    def abs_curve(x_arr, gamma, const1, const2):
#        y_model = const1 * np.exp(-gamma * x_arr) + const2
#
#        # Return chi square
#        return np.sum((s_arr - y_model) ** 2)
#
#    pdb.set_trace()
#
#    init_guess = np.array([3, c_abs[0], np.average(c_arr)])
#    out, cov = curve_fit(abs_curve, s_arr, c_abs, p0=init_guess)
#
#    gamma_abs, const1_abs, const2_abs = out
#    y_abs = const1_abs * np.exp(-gamma_abs * x_arr) + const2_abs
#
#    plt.close()
#    plt.plot(s_arr, c_abs)
#    plt.plot(s_arr, y_abs)
#    plt.show()
#
#    # Define functional form
#    def curve(x_arr, gamma, alpha, const):
#        y_model = (1 - const) * np.exp(-gamma * x_arr) * np.cos(alpha * np.pi * x_arr) + const
#
#        # Return chi square
#        return np.sum((s_arr - y_model) ** 2)
#
#    init_guess = np.array([3, np.pi, np.average(c_arr)])
#    out, cov = curve_fit(curve, s_arr, c_arr, p0=init_guess)
#
#    gamma_out, alpha_out, const_out = out
#    gamma_err, alpha_err, const_err = np.sqrt(np.diag(cov))
#    gamma_out, alpha_out, const_out = init_guess
#    y_out = (1 - const_out) * np.exp(-gamma_out * s_arr) * np.cos(alpha_out * np.pi * s_arr) + const_out
#
#    plt.close()
#    plt.plot(s_arr, c_arr)
#    plt.plot(s_arr, y_out)
#    plt.show()





















































#
