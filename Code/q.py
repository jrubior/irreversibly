# --- Packages
import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe
import os

from numba import jit, float64, int64, prange
from numba.experimental import jitclass

os.chdir(os.path.dirname(os.path.abspath(__file__)))
print(f"Current working directory: {os.getcwd()}") 
CODE_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(CODE_DIR, '..', 'Output')

# ============================================================
# TWO-SECTOR MODEL WITH ONE SHOCK (z1)
# Q-Learning
#
# State space  : (k1, k2, z1)             -- 3-dimensional, same as VFI
# Actions      : (k1_next, k2_next)       -- next-period capitals (indices)
# Q-table      : Q[i1, i2, iz, j1, j2]   -- value of action (j1,j2) in state (i1,i2,iz)
#
# Key difference from VFI
# -----------------------
# VFI computes the expected continuation value analytically:
#     ev = sum_jz  P_z1[iz, jz] * V[j1, j2, jz]
#
# Q-learning instead draws a *realised* next shock iz_next by sampling
# from the transition row P_z1[iz, :] and uses:
#     TD = r + β * max_{j1',j2'} Q[j1, j2, iz_next, j1', j2'] - Q[i1, i2, iz, j1, j2]
#
# No matrix multiply is needed; the agent learns the transition structure
# implicitly through repeated sampling.
#
# The rest of the economics (preferences, production, resource constraints)
# is identical to vfi.py.  P_z1 is still passed in, but only used as a
# sampler, not as a weighted sum.
# ============================================================


# ------------------------------------------------------------
# Model class  (identical to vfi.py; reproduced for self-containedness)
# ------------------------------------------------------------

model_data = [
    # Preferences
    ('β',       float64),
    ('σ',       float64),
    ('θ1',      float64),
    ('θ2',      float64),
    # Production
    ('α1',      float64),
    ('α2',      float64),
    ('γ1',      float64),
    ('γ2',      float64),
    # Depreciation
    ('δ1',      float64),
    ('δ2',      float64),
    # Materials shares
    ('φ11',     float64),
    ('φ21',     float64),
    ('φ12',     float64),
    ('φ22',     float64),
    # Investment shares
    ('ω11',     float64),
    ('ω21',     float64),
    ('ω12',     float64),
    ('ω22',     float64),
    # Shock process
    ('ρ_z1',    float64),
    ('σ_z1',    float64),
    # Grids
    ('k1_grid', float64[:]),
    ('k2_grid', float64[:]),
    ('z1_grid', float64[:]),
    ('P_z1',    float64[:, :]),   # Markov transition matrix for z1 (used as sampler)
    # Steady-state TFP scalars (z2, l1, l2 fixed at SS=0)
    ('A1_ss',   float64),
    ('A2_ss',   float64),
]

@jitclass(model_data)
class TwoSector:
    def __init__(self,
                 β=0.95,#0.99,
                 σ=2.0,
                 θ1=0.5,   θ2=0.5,
                 α1=0.33,  α2=0.33,
                 γ1=0.5,   γ2=0.5,
                 δ1=0.025, δ2=0.025,
                 φ11=0.5,  φ21=0.5,
                 φ12=0.5,  φ22=0.5,
                 ω11=0.5,  ω21=0.5,
                 ω12=0.5,  ω22=0.5,
                 ρ_z1=0.90,
                 σ_z1=0.01,
                 nk1=20,
                 nk2=20,
                 k_width=0.10,
                 z1_grid=np.zeros(7),
                 P_z1=np.zeros((7, 7))):

        self.β  = β;  self.σ  = σ
        self.θ1 = θ1; self.θ2 = θ2
        self.α1 = α1; self.α2 = α2
        self.γ1 = γ1; self.γ2 = γ2
        self.δ1 = δ1; self.δ2 = δ2
        self.φ11 = φ11; self.φ21 = φ21
        self.φ12 = φ12; self.φ22 = φ22
        self.ω11 = ω11; self.ω21 = ω21
        self.ω12 = ω12; self.ω22 = ω22
        self.ρ_z1 = ρ_z1; self.σ_z1 = σ_z1

        self.A1_ss = 1.0
        self.A2_ss = 1.0

        k1_ss = self._kstar(α1, δ1, γ1, self.A1_ss)
        k2_ss = self._kstar(α2, δ2, γ2, self.A2_ss)

        self.k1_grid = np.linspace((1 - k_width) * k1_ss,
                                   (1 + k_width) * k1_ss, nk1)
        self.k2_grid = np.linspace((1 - k_width) * k2_ss,
                                   (1 + k_width) * k2_ss, nk2)

        self.z1_grid = z1_grid
        self.P_z1    = P_z1

    def _kstar(self, α, δ, γ, A):
        kv = α * self.β / (1 - self.β * (1 - δ))
        k  = (kv * A * α**(-α))**(1 / (1 - α))
        return k

    def value_added(self, k, α, A):
        return A * (k / α)**α

    def output(self, v, γ):
        return v / γ

    def u(self, C):
        if self.σ == 1.0:
            return np.log(C)
        else:
            return (C**(1 - self.σ) - 1) / (1 - self.σ)

    def composite_C(self, c1, c2):
        return (c1 / self.θ1)**self.θ1 * (c2 / self.θ2)**self.θ2


# ------------------------------------------------------------
# Reward: flow utility for a given state-action pair
# ------------------------------------------------------------

@jit(nopython=True)
def reward(i1, i2, iz, j1, j2, mod):
    '''
    Compute the immediate reward r(s, a) for state (i1, i2, iz) and
    action (j1, j2).  Returns -inf for infeasible choices.

    This is the same feasibility and production block as bellman_rhs
    in vfi.py, but WITHOUT the continuation value term -- that is
    handled separately in the TD update using the sampled next state.
    '''
    k1      = mod.k1_grid[i1]
    k2      = mod.k2_grid[i2]
    z1      = mod.z1_grid[iz]
    k1_next = mod.k1_grid[j1]
    k2_next = mod.k2_grid[j2]

    x1 = k1_next - (1 - mod.δ1) * k1
    x2 = k2_next - (1 - mod.δ2) * k2

    # Irreversibility: uncomment when ready
    # if x1 < 0.0 or x2 < 0.0:
    #     return -np.inf

    # TFP -- note exp(z1): z1_grid holds log-levels
    A1 = mod.A1_ss * np.exp(z1)
    A2 = mod.A2_ss

    v1 = mod.value_added(k1, mod.α1, A1)
    v2 = mod.value_added(k2, mod.α2, A2)
    y1 = mod.output(v1, mod.γ1)
    y2 = mod.output(v2, mod.γ2)

    m1  = (1 - mod.γ1) * y1;  m2  = (1 - mod.γ2) * y2
    m11 = mod.φ11 * m1;       m21 = mod.φ21 * m1
    m12 = mod.φ12 * m2;       m22 = mod.φ22 * m2

    x11 = mod.ω11 * x1;  x21 = mod.ω21 * x1
    x12 = mod.ω12 * x2;  x22 = mod.ω22 * x2

    c1 = y1 - m11 - m12 - x11 - x12
    c2 = y2 - m21 - m22 - x21 - x22

    if c1 <= 0.0 or c2 <= 0.0:
        return -np.inf

    C = mod.composite_C(c1, c2)
    return mod.u(C)


# ------------------------------------------------------------
# Precompute reward table: R[i1, i2, iz, j1, j2]
# ------------------------------------------------------------

@jit(nopython=True, parallel=True)
def compute_reward_table(mod):
    '''
    Build the full reward table once before Q-learning starts.
    Subsequent training steps do O(1) lookups instead of recomputing
    exp, pow, and feasibility checks each time.
    Shape mirrors the Q-table: (nk1, nk2, nz, nk1, nk2).
    Memory note: same size as Q (~800 MB for nk1=nk2=100, nz=1).
    '''
    nk1 = len(mod.k1_grid)
    nk2 = len(mod.k2_grid)
    nz  = len(mod.z1_grid)
    R = np.full((nk1, nk2, nz, nk1, nk2), -np.inf)
    for i1 in prange(nk1):
        for i2 in range(nk2):
            for iz in range(nz):
                for j1 in range(nk1):
                    for j2 in range(nk2):
                        R[i1, i2, iz, j1, j2] = reward(i1, i2, iz, j1, j2, mod)
    return R


# ------------------------------------------------------------
# Sample next shock index from transition row (used instead of P_z1 weighted sum)
# ------------------------------------------------------------

@jit(nopython=True)
def sample_next_z(iz, mod):
    '''
    Draw iz_next ~ Multinomial(P_z1[iz, :]).
    This replaces the weighted sum  ev = sum_jz P_z1[iz,jz]*V[j1,j2,jz]
    used in VFI.  The agent never sees the probabilities explicitly;
    it only observes the realised next state.
    '''
    u   = np.random.rand()
    cdf = np.cumsum(mod.P_z1[iz])
    jz  = np.searchsorted(cdf, u)
    return min(jz, len(mod.z1_grid) - 1)   # clamp for floating-point safety


# ------------------------------------------------------------
# update_step: one episode of Q-learning
# ------------------------------------------------------------

@jit(nopython=True)
def update_step(Q, R, mod, ε=0.1, lr=0.1):
    '''
    Run one episode of Q-learning: randomly pick a starting state,
    then take steps until the within-episode Q-table change is below
    a tight threshold or a step limit is reached.

    This is the Q-learning analogue of the Bellman sweep in vfi.py.
    Rather than updating every state simultaneously (as VFI does),
    Q-learning updates one (state, action) cell at a time using the
    TD error computed from the *sampled* next state.

    Parameters
    ----------
    Q   : current Q-table, shape (nk1, nk2, nz, nk1, nk2)
              Q[i1, i2, iz, j1, j2] = value of choosing (j1,j2) in state (i1,i2,iz)
    R   : precomputed reward table, same shape as Q
    mod : TwoSector instance
    ε   : exploration probability (epsilon-greedy)
    lr  : learning rate α_t (held fixed per episode; caller may decay across episodes)

    Returns
    -------
    Q      : updated Q-table (in-place modification + return for clarity)
    max_td : max |TD error| seen this episode (used as convergence signal)
    '''
    nk1 = len(mod.k1_grid)
    nk2 = len(mod.k2_grid)
    nz  = len(mod.z1_grid)

    # Suggestion 3: Increased steps from 20,000 to 30,000
    max_steps = 30000
    max_td    = 0.0

    i1 = np.random.randint(0, nk1)
    i2 = np.random.randint(0, nk2)
    iz = np.random.randint(0, nz)

    for _ in range(max_steps):
        # ε-greedy selection
        if np.random.rand() < ε:
            j1 = np.random.randint(0, nk1)
            j2 = np.random.randint(0, nk2)
        else:
            best_val = -np.inf
            j1 = 0; j2 = 0
            for a1 in range(nk1):
                for a2 in range(nk2):
                    if Q[i1, i2, iz, a1, a2] > best_val:
                        best_val = Q[i1, i2, iz, a1, a2]
                        j1 = a1; j2 = a2

        r = R[i1, i2, iz, j1, j2]

        if r == -np.inf:
            j1 = np.random.randint(0, nk1)
            j2 = np.random.randint(0, nk2)
            r  = R[i1, i2, iz, j1, j2]
            if r == -np.inf:
                continue

        iz_next = sample_next_z(iz, mod)

        max_q_next = -np.inf
        for a1 in range(nk1):
            for a2 in range(nk2):
                q_val = Q[j1, j2, iz_next, a1, a2]
                if q_val > max_q_next:
                    max_q_next = q_val

        TD = r + mod.β * max_q_next - Q[i1, i2, iz, j1, j2]
        Q[i1, i2, iz, j1, j2] += lr * TD

        abs_td = TD if TD >= 0.0 else -TD
        if abs_td > max_td:
            max_td = abs_td

        i1, i2, iz = j1, j2, iz_next

    return Q, max_td



# ------------------------------------------------------------
# qlearning: run episodes until the Q-table converges
# ------------------------------------------------------------

def qlearning(mod,
              N=500000,
              eps=None,
              eps_0=0.50,      # Start with more exploration
              callback=None,
              eps_min=0.01,
              eps_decay=1e-3,  # Slower decay over 500k episodes
              lr_0=0.50,
              lr_min=0.001,    # Suggestion 1: Lowered floor from 0.01
              lr_decay=1e-3,   # Slower decay
              tol=1e-5,
              print_step=1000):
    '''
    Run N episodes of Q-learning, decaying the learning rate across
    episodes, and stop early if the max |TD error| within an episode
    falls below tol.

    Convergence criterion: max |TD| across steps in the episode, which
    measures Bellman residuals directly and avoids copying the Q-table.

    This is the Q-learning analogue of vfi() in vfi.py.

    Parameters
    ----------
    mod        : TwoSector instance
    N          : maximum number of episodes
    ε          : exploration rate (fixed across episodes)
    lr_0       : initial learning rate
    lr_min     : floor for the decayed learning rate
    lr_decay   : exponential decay constant: lr_t = max(lr_0*exp(-decay*t), lr_min)
    tol        : convergence tolerance on max |TD error| within an episode
    print_step : print progress every this many episodes

    Returns
    -------
    Q      : converged Q-table, shape (nk1, nk2, nz, nk1, nk2)
    pol_k1 : greedy policy for k1_next (index), shape (nk1, nk2, nz)
    pol_k2 : greedy policy for k2_next (index), shape (nk1, nk2, nz)
    '''

    if eps is not None:
        eps_0 = eps
    nk1, nk2, nz = len(mod.k1_grid), len(mod.k2_grid), len(mod.z1_grid)
    Q = np.zeros((nk1, nk2, nz, nk1, nk2))

    print('Precomputing reward table ...')
    R = compute_reward_table(mod)
    print('Done.\n')

    for episode in range(1, N + 1):
        # Suggestion 2: Decaying learning rate and epsilon
        lr  = max(lr_0  * np.exp(-lr_decay  * episode), lr_min)
        eps = max(eps_0 * np.exp(-eps_decay * episode), eps_min)

        Q, error = update_step(Q, R, mod, ε=eps, lr=lr)

        if episode % print_step == 0:
            print(f'Ep {episode:6d} | lr: {lr:.4f} | eps: {eps:.4f} | Error: {error:.2e}')
            if callback is not None:
                callback(episode, Q, mod)

        if error < tol:
            print(f'Converged in {episode} episodes')
            break
    else:
        print('Warning: reached maximum episodes without convergence')

    flat   = Q.reshape(nk1, nk2, nz, nk1 * nk2)
    best   = flat.argmax(axis=3)
    pol_k1 = (best // nk2).astype(np.int64)
    pol_k2 = (best  % nk2).astype(np.int64)

    return Q, pol_k1, pol_k2


def vf_from_Q(Q):
    '''
    Recover the value function from the Q-table by taking the max over actions.
    Analogous to vf() in qlearning.py.
    Shape: (nk1, nk2, nz, nk1, nk2) -> (nk1, nk2, nz)
    '''
    return Q.max(axis=(3, 4))

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------

if __name__ == '__main__':

    # --- Tauchen discretisation (identical to vfi.py)
    ρ_z1, σ_z1, nz = 0.90, 0.01, 2
    mc      = qe.markov.approximation.tauchen(nz, ρ_z1, σ_z1)
    z1_grid = mc.state_values.astype(np.float64)
    P_z1    = np.asarray(mc.P, dtype=np.float64)

    # # --- Single z=0 (SS), single k2 (SS), large k1 grid to test convergence
    # z1_grid = np.array([0.0])
    # P_z1    = np.array([[1.0]])

    # --- Instantiate model (same parameters as vfi.py)
    mod = TwoSector(nk1=30, nk2=30, z1_grid=z1_grid, P_z1=P_z1)

    print('Grid sizes:')
    print(f'  k1: {len(mod.k1_grid)} points in '
          f'[{mod.k1_grid[0]:.3f}, {mod.k1_grid[-1]:.3f}]')
    print(f'  k2: {len(mod.k2_grid)} points in '
          f'[{mod.k2_grid[0]:.3f}, {mod.k2_grid[-1]:.3f}]')
    print(f'  z1: {len(mod.z1_grid)} points in '
          f'[{mod.z1_grid[0]:.4f}, {mod.z1_grid[-1]:.4f}]')

    # --- Live plot setup
    iz_mid      = len(mod.z1_grid) // 2
    nk2_n       = len(mod.k2_grid)
    ik2_indices = [nk2_n // 4, nk2_n // 2, 3 * nk2_n // 4]
    colors      = ['steelblue', 'darkorange', 'forestgreen']

    def live_plot(episode, Q, mod):
        v      = vf_from_Q(Q)
        flat   = Q.reshape(len(mod.k1_grid), len(mod.k2_grid), len(mod.z1_grid), -1)
        pol_k1 = flat.argmax(axis=3) // len(mod.k2_grid)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        for ik2, col in zip(ik2_indices, colors):
            axes[0].plot(mod.k1_grid, v[:, ik2, iz_mid], color=col,
                         label=f'k2={mod.k2_grid[ik2]:.2f}')
            axes[1].plot(mod.k1_grid, mod.k1_grid[pol_k1[:, ik2, iz_mid]], color=col,
                         label=f'k2={mod.k2_grid[ik2]:.2f}')
        axes[0].set(xlabel='k1', ylabel='V',
                    title=f'Value function  (ep {episode}, z1={mod.z1_grid[iz_mid]:.3f})')
        axes[1].set(xlabel='k1', ylabel="k1'",
                    title=f"Policy k1'  (ep {episode}, z1={mod.z1_grid[iz_mid]:.3f})")
        for ax in axes:
            ax.legend()
        fig.tight_layout()
        plt.show()

    # --- Solve
    print('\nStarting Q-learning ...\n')
    Q, pol_k1, pol_k2 = qlearning(mod, N=500000, eps=0.50, eps_decay=1e-5,
                                   lr_0=0.50, lr_decay=1e-5,
                                   lr_min=0.0001,
                                   tol=1e-5, print_step=500,
                                   callback=live_plot)

    plt.savefig('qlearning_two_sector.png', dpi=150)
    plt.show()
    print('Figure saved to qlearning_two_sector.png')