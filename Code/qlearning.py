# --- Packages 
import numpy as np
import matplotlib.pyplot as plt

from numba import jit, float64, int64, njit, prange
from numba.experimental import jitclass

#------------------------------------------------------------
# 1. DYNAMIC PROGRAMMING APPROACH
#------------------------------------------------------------

# --- Define RBC class

rbc_data = [
    ('α', float64),          # Production parameter
    ('β', float64),          # Discount factor
    ('δ', float64),          # Depreciation rate
    ('σ', float64),          # Risk aversion parameter
    ('grid', float64[:]),    # Grid (array)
    ('kstar', float64),      # Steady-state capital
    ('kmin', float64),       # Minimum grid value
    ('kmax', float64)        # Maximum grid value
]

@jitclass(rbc_data)
class RBC:
    def __init__(self,
        α=0.33,
        β=0.95,
        δ=0.10,
        σ=2.00,
        grid_size=100):

        self.α, self.β, self.δ, self.σ = α, β, δ, σ

        # Set up grid
        self.kstar = (α / (1 / β - (1 - δ))) ** (1 / (1 - α))
        self.kmin = 0.25 * self.kstar
        self.kmax = 1.75 * self.kstar
        self.grid = np.linspace(self.kmin, self.kmax, grid_size)

    def f(self, k):
        "The production function"
        return k**self.α

    def u(self, c):
        "The utility function"
        if self.σ == 1:
            return np.log(c)
        else:
            return (c**(1 - self.σ) - 1) / (1 - self.σ)

#--- Value of consumption choice c given capital k

@jit
def state_action_value(i, j, v, rbc):
    '''
    Right hand side of the Bellman equation for given state and action indices
    * rbc is an instance of the RBC model
    * c is consumption
    * k is current capital
    * k_next is next period capital
    '''

    α, β, δ, σ, k_grid = rbc.α, rbc.β, rbc.δ, rbc.σ, rbc.grid
    k, k_next = k_grid[i], k_grid[j]
    c = rbc.f(k) + (1 - δ) * k - k_next
    value = -np.inf
    if c > 0:
        value = rbc.u(c) + β * v[j]
    return value

# --- Bellman operator

@jit(parallel=True)
def T(v, rbc):
    '''
    The Bellman operator.
     * rbc is an instance of RBCModel
     * v is an array representing a guess of the value function
     * v_new is the updated estimate of the value function
    '''

    α, β, δ, σ, k_grid = rbc.α, rbc.β, rbc.δ, rbc.σ, rbc.grid
    v_new = np.empty_like(v)
    for i in prange(len(k_grid)):
        x_tmp = np.array([state_action_value(i, j, v, rbc) for j in np.arange(len(k_grid))])
        v_new[i] = np.max(x_tmp)
    return v_new

#--- Value function iteration

def vfi(rbc, tol=1e-5, max_iter=1e5, print_step=25):
    '''
    Value function iteration
    '''
    v = np.zeros(len(rbc.grid))
    error = tol + 1
    iter = 1
    while (error > tol ) & (iter < max_iter):
        v_new = T(v, rbc)
        error = np.max(np.abs(v_new - v))
        if iter % print_step == 0:
            print(f'Error at iteration {iter} is {error}')
        v = v_new
        iter += 1
    if error < tol:
        print(f'Converged in {iter} iterations')
    else:
        print('Failed to converge')
    return v


#------------------------------------------------------------
# 2. Q-LEARNING APPROACH
#------------------------------------------------------------

# --- Define QRBC class

qrbc_data = [
    ('α', float64),          # Production parameter
    ('β', float64),          # Discount factor
    ('δ', float64),          # Depreciation rate
    ('σ', float64),          # Risk aversion parameter
    ('grid_size', int64),    # Grid size (number of grid points)
    ('grid', float64[:]),    # Grid (array)
    ('kstar', float64),      # Steady-state capital
    ('kmin', float64),       # Minimum grid value
    ('kmax', float64),       # Maximum grid value
    ('ε', float64),          # Exploration rate
    ('δ_tol', float64),      # Tolerance for convergence
    ('N', int64),            # Number of episodes
    ('qtable', float64[:, :]),  # Q-table
]

@jitclass(qrbc_data)
class QRBC:
    def __init__(self,
                 α=0.33, 
                 β=0.95, 
                 δ=0.10, 
                 σ=2.0, 
                 grid_size=10, 
                 ε=0.1, 
                 δ_tol=1e-5, 
                 N=20000,
                 seed=None):
        
        self.α, self.β, self.δ, self.σ = α, β, δ, σ
        self.grid_size = grid_size
        self.ε = ε
        self.δ_tol = δ_tol
        self.N = N

        # Set up capital grid
        self.kstar = (α / (1 / β - (1 - δ))) ** (1 / (1 - α))
        self.grid = np.linspace(0.25 * self.kstar, 1.75 * self.kstar, grid_size)

        if seed is not None:
            np.random.seed(seed)

    def f(self, k):
        '''Define the production function'''
        return k**self.α

    def u(self, c):
        '''Define the utility function'''
        if self.σ == 1:
            return np.log(c)
        else:
            return (c**(1 - self.σ) - 1) / (1 - self.σ)

    def select_state(self):
        return np.random.randint(0, len(self.grid))

    def temp_diff(self, state, action, t):
        '''
        Compute the temporal difference for a given state-action pair
        '''
        state_next = np.argmin(np.abs(self.grid - self.grid[action]))
        c = max(self.f(self.grid[state]) + (1 - self.δ) * self.grid[state] - self.grid[action], 1e-8)
        TD = self.u(c) + self.β * np.max(self.qtable[state_next, :]) - self.qtable[state, action]
        return TD, state_next

    def run_one_step(self, max_times=15000):
        '''Run a single step in an episode'''
        s = self.select_state()
        for t in range(1, max_times + 1):
            # Choose action
            if np.random.rand() < self.ε:
                action = np.random.randint(0, self.grid_size)  # Explore
            else:
                action = np.argmax(self.qtable[s, :])  # Exploit

            TD, s_next = self.temp_diff(s, action, t)

            # Update Q-table with decaying learning rate
            qtable_new = self.qtable.copy()
            lr_t = 1 / t
            qtable_new[s, action] = self.qtable[s, action] + lr_t * TD

            # Check convergence
            error = np.abs(qtable_new - self.qtable).max()
            if error < self.δ_tol:
                break
            s, self.qtable = s_next, qtable_new

        return self.qtable

@jit
def run_episodes(qrbc):
    '''
    Run N episodes of Q-learning with qtable of last iteration each time
    '''
    N = qrbc.N
    for i in range(N):
        if i % (N / 10) == 0:
            print(f"Progress: EPOCHS = {i}")
        new_qtable = qrbc.run_one_step()

    return new_qtable

def vf(qtable):
    '''Compute the value function from the Q-table'''
    v = np.max(qtable, axis=1)
    return v