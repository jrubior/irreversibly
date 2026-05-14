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
# Value Function Iteration
#
# State space : (k1, k2, z1)   -- 3-dimensional
# Actions     : (k1_next, k2_next) -- next-period capital in each sector
#
# Model equations (logs dropped; all variables in levels):
#
#   Preferences
#     U(C) = C^(1-σ)/(1-σ),  C = (c1/θ1)^θ1 * (c2/θ2)^θ2
#
#   Production (value-added only, materials = (1-γ)*y)
#     v1 = A1 * (k1/α1)^α1,   A1 = exp(z1) * (l1_ss/(1-α1))^(1-α1)
#     v2 = A2 * (k2/α2)^α2,   A2 =            (l2_ss/(1-α2))^(1-α2)
#     y_i = v_i / γ_i          (from FOC v_i = γ_i * y_i)
#
#   Materials (from FOC m_i = (1-γ_i)*y_i)
#     m_i = (1-γ_i) * y_i
#     m_ij = φ_ij * m_j        (sectoral allocation)
#
#   Investment aggregates
#     x_ij = ω_ij * x_j        (sectoral allocation)
#
#   Resource constraints
#     c1 = y1 - m11 - m12 - x11 - x12
#     c2 = y2 - m21 - m22 - x21 - x22
#
#   Capital accumulation
#     k1' = x1 + (1-δ1)*k1
#     k2' = x2 + (1-δ2)*k2
#     => x_i = k_i' - (1-δ_i)*k_i   (action choice pins down investment)
#
#   Shock process (Tauchen discretisation)
#     z1' = ρ_z1 * z1 + ε,  ε ~ N(0, σ_z1^2)
# ============================================================


# ------------------------------------------------------------
# Model class
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
    ('P_z1',    float64[:, :]),   # Markov transition matrix for z1
    # Steady-state TFP scalars (z2, l1, l2 fixed at SS=0)
    ('A1_ss',   float64),         # TFP factor for sector 1 excluding z1
    ('A2_ss',   float64),         # TFP factor for sector 2 (constant)
]

@jitclass(model_data)
class TwoSector:
    def __init__(self,
                 β=0.95,
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
                 k_width=0.05,     # grid half-width around SS as a fraction
                 z1_grid=np.zeros(7),   # shock grid  (pass from qe.tauchen)
                 P_z1=np.zeros((7, 7))): # transition matrix (pass from qe.tauchen)
        
        # --- Store parameters
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

        # --- Steady-state TFP factors (z2=l1=l2=0 => log-level = 0)
        # A_i = exp(z_i) * (l_i/(1-α_i))^(1-α_i)
        # At SS: exp(z_i)=1, l_i_ss=1-α_i => (l_i_ss/(1-α_i))^(1-α_i)=1
        self.A1_ss = 1.0   # multiplied by exp(z1) at runtime
        self.A2_ss = 1.0

        # --- Steady-state capital (from Euler + production, symmetric SS)
        # r_i = α_i * A_i * (k_i/α_i)^(α_i-1) => k_i/v_i = α_i*β/(1-β*(1-δ_i))
        # Solve for k_i numerically via fixed point (simple here due to symmetry)
        k1_ss = self._kstar(α1, δ1, γ1, self.A1_ss)
        k2_ss = self._kstar(α2, δ2, γ2, self.A2_ss)

        # --- Capital grids (symmetric around SS)
        self.k1_grid = np.linspace((1 - k_width) * k1_ss,
                                   (1 + k_width) * k1_ss, nk1)
        self.k2_grid = np.linspace((1 - k_width) * k2_ss,
                                   (1 + k_width) * k2_ss, nk2)

        # --- Shock grid and transition matrix (supplied by caller via qe.tauchen)
        self.z1_grid = z1_grid
        self.P_z1    = P_z1

    # ----------------------------------------------------------
    # Helpers used inside __init__ (must be defined as methods)
    # ----------------------------------------------------------

    def _kstar(self, α, δ, γ, A):
        '''
        Steady-state capital for sector i.
        From Euler:  px_i = α_i * pv_i * v_i / k_i  =>  k_i/v_i = α_i*β/(1-β*(1-δ_i))
        From production: v_i = A_i*(k_i/α_i)^α_i
        Combine to get an implicit equation solved by iteration.
        '''
        kv = α * self.β / (1 - self.β * (1 - δ))   # k/v ratio
        # v = A*(k/α)^α  and  k = kv*v  =>  k = kv*A*(k/α)^α
        # k^(1-α) = kv * A * α^(-α)
        k = (kv * A * α**(-α))**(1 / (1 - α))
        return k

    # ----------------------------------------------------------
    # Production and utility
    # ----------------------------------------------------------

    def value_added(self, k, α, A):
        '''v_i = A_i * (k_i / α_i)^α_i'''
        return A * (k / α)**α

    def output(self, v, γ):
        '''y_i = v_i / γ_i  (from FOC: v_i = γ_i * y_i at interior optimum)'''
        return v / γ

    def u(self, C):
        '''Aggregate CRRA utility over composite consumption C.'''
        if self.σ == 1.0:
            return np.log(C)
        else:
            return (C**(1 - self.σ) - 1) / (1 - self.σ)

    def composite_C(self, c1, c2):
        '''C = (c1/θ1)^θ1 * (c2/θ2)^θ2'''
        return (c1 / self.θ1)**self.θ1 * (c2 / self.θ2)**self.θ2


# ------------------------------------------------------------
# Bellman right-hand side for a single state-action tuple
# ------------------------------------------------------------

@jit(nopython=True)
def bellman_rhs(i1, i2, iz, j1, j2, v, mod):
    '''
    Compute the right-hand side of the Bellman equation for:
      state  = (k1_grid[i1], k2_grid[i2], z1_grid[iz])
      action = (k1_grid[j1], k2_grid[j2])  -- next-period capitals

    Returns -inf for infeasible (negative consumption) choices.

    Parameters
    ----------
    i1, i2 : current capital indices
    iz     : current z1 index
    j1, j2 : next-period capital indices (action)
    v      : current value function array, shape (nk1, nk2, nz)
    mod    : TwoSector instance
    '''
    k1      = mod.k1_grid[i1]
    k2      = mod.k2_grid[i2]
    z1      = mod.z1_grid[iz]
    k1_next = mod.k1_grid[j1]
    k2_next = mod.k2_grid[j2]

    # --- Investment implied by the action
    x1 = k1_next - (1 - mod.δ1) * k1
    x2 = k2_next - (1 - mod.δ2) * k2

    # Irreversibility would go here:  if x1 < 0 or x2 < 0: return -inf
    # (left as a comment for the next stage of the project)

    # --- TFP (z2, l1, l2 fixed at 0 => A2_ss already constant)
    A1 = mod.A1_ss * z1
    A2 = mod.A2_ss            # z2=l2=0 throughout

    # --- Value added and output
    v1 = mod.value_added(k1, mod.α1, A1)
    v2 = mod.value_added(k2, mod.α2, A2)
    y1 = mod.output(v1, mod.γ1)
    y2 = mod.output(v2, mod.γ2)

    # --- Materials aggregates: m_j = (1-γ_j)*y_j, then m_ij = φ_ij*m_j
    m1 = (1 - mod.γ1) * y1
    m2 = (1 - mod.γ2) * y2
    m11 = mod.φ11 * m1;  m21 = mod.φ21 * m1
    m12 = mod.φ12 * m2;  m22 = mod.φ22 * m2

    # --- Sectoral investment: x_ij = ω_ij * x_j
    x11 = mod.ω11 * x1;  x21 = mod.ω21 * x1
    x12 = mod.ω12 * x2;  x22 = mod.ω22 * x2

    # --- Resource constraints
    c1 = y1 - m11 - m12 - x11 - x12
    c2 = y2 - m21 - m22 - x21 - x22

    if c1 <= 0.0 or c2 <= 0.0:
        return -np.inf

    # --- Composite consumption and flow utility
    C      = mod.composite_C(c1, c2)
    flow_u = mod.u(C)

    # --- Expected continuation value: E[V(k1',k2',z1') | z1]
    ev = 0.0
    for jz in range(len(mod.z1_grid)):
        ev += mod.P_z1[iz, jz] * v[j1, j2, jz]

    return flow_u + mod.β * ev


# ------------------------------------------------------------
# update_step: one full sweep of the Bellman operator
# ------------------------------------------------------------

@jit(nopython=True, parallel=True)
def update_step(v, mod):
    '''
    Apply one Bellman update to the entire state space.

    For each state (i1, i2, iz) maximise over all feasible actions
    (j1, j2) -- the next-period capital grid points -- and return
    the updated value function v_new and the greedy policy arrays.

    Parameters
    ----------
    v   : current value function, shape (nk1, nk2, nz)
    mod : TwoSector instance

    Returns
    -------
    v_new      : updated value function, same shape as v
    pol_k1     : greedy policy for k1_next (index into k1_grid)
    pol_k2     : greedy policy for k2_next (index into k2_grid)
    '''
    nk1 = len(mod.k1_grid)
    nk2 = len(mod.k2_grid)
    nz  = len(mod.z1_grid)

    v_new  = np.empty((nk1, nk2, nz))
    pol_k1 = np.empty((nk1, nk2, nz), dtype=int64)
    pol_k2 = np.empty((nk1, nk2, nz), dtype=int64)

    # Parallelise over the outer (k1, z1) loop; inner k2 loop is serial
    for i1 in prange(nk1):
        for iz in range(nz):
            for i2 in range(nk2):
                best_val = -np.inf
                best_j1  = 0
                best_j2  = 0
                for j1 in range(nk1):
                    for j2 in range(nk2):
                        val = bellman_rhs(i1, i2, iz, j1, j2, v, mod)
                        if val > best_val:
                            best_val = val
                            best_j1  = j1
                            best_j2  = j2
                v_new[i1, i2, iz]  = best_val
                pol_k1[i1, i2, iz] = best_j1
                pol_k2[i1, i2, iz] = best_j2

    return v_new, pol_k1, pol_k2


# ------------------------------------------------------------
# vfi: iterate update_step until convergence
# ------------------------------------------------------------

def vfi(mod, tol=1e-5, max_iter=10_000, print_step=25):
    '''
    Value function iteration for the two-sector model.

    Repeatedly calls update_step until the sup-norm distance between
    successive value functions falls below tol.

    Parameters
    ----------
    mod        : TwoSector instance
    tol        : convergence tolerance (sup-norm), default 1e-5
    max_iter   : maximum number of iterations, default 10_000
    print_step : print progress every this many iterations, default 25

    Returns
    -------
    v      : converged value function, shape (nk1, nk2, nz)
    pol_k1 : optimal next-period k1 index, shape (nk1, nk2, nz)
    pol_k2 : optimal next-period k2 index, shape (nk1, nk2, nz)
    '''
    nk1 = len(mod.k1_grid)
    nk2 = len(mod.k2_grid)
    nz  = len(mod.z1_grid)

    # Initial guess: zero value function
    v = np.zeros((nk1, nk2, nz))

    error    = tol + 1.0
    iteration = 1

    while error > tol and iteration <= max_iter:
        v_new, pol_k1, pol_k2 = update_step(v, mod)
        error = np.max(np.abs(v_new - v))

        if iteration % print_step == 0:
            print(f'Iteration {iteration:5d} | Error: {error:.2e}')

        v = v_new
        iteration += 1

    if error < tol:
        print(f'Converged in {iteration - 1} iterations')
    else:
        print('Warning: failed to converge')

    return v, pol_k1, pol_k2


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------

if __name__ == '__main__':

    # --- Tauchen discretisation of z1 (done outside jitclass via quantecon)
    # ρ_z1, σ_z1, nz = 0.90, 0.01, 7
    # mc = qe.markov.approximation.tauchen(nz, ρ_z1, σ_z1)
    # z1_grid = np.exp(mc.state_values.astype(np.float64))
    # P_z1    = np.asarray(mc.P, dtype=np.float64)

    # --- Single z=1 (SS level), single k2 (SS), large k1 grid
    z1_grid = np.array([1.0])
    P_z1    = np.array([[1.0]])

    # --- Instantiate model
    mod = TwoSector(nk1=100, nk2=10, z1_grid=z1_grid, P_z1=P_z1)

    print('Grid sizes:')
    print(f'  k1: {len(mod.k1_grid)} points in '
          f'[{mod.k1_grid[0]:.3f}, {mod.k1_grid[-1]:.3f}]')
    print(f'  k2: {len(mod.k2_grid)} points in '
          f'[{mod.k2_grid[0]:.3f}, {mod.k2_grid[-1]:.3f}]')
    print(f'  z1: {len(mod.z1_grid)} points in '
          f'[{mod.z1_grid[0]:.4f}, {mod.z1_grid[-1]:.4f}]')

    # --- Solve
    print('\nStarting VFI ...\n')
    v, pol_k1, pol_k2 = vfi(mod, tol=1e-5, print_step=10)

    # --- Plot value function at median z1
    iz_mid = len(mod.z1_grid) // 2
    nk2 = len(mod.k2_grid)
    ik2_low  = nk2 // 4
    ik2_mid  = nk2 // 2
    ik2_high = 3 * nk2 // 4
    ik2_indices = [ik2_low, ik2_mid, ik2_high]
    colors = ['steelblue', 'darkorange', 'forestgreen']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Value function lines
    ax = axes[0]
    for ik2, col in zip(ik2_indices, colors):
        ax.plot(mod.k1_grid, v[:, ik2, iz_mid], color=col,
                label=f'k2 = {mod.k2_grid[ik2]:.3f}')
    ax.set_xlabel('k1'); ax.set_ylabel('V')
    ax.set_title(f'Value function  (z1 = {mod.z1_grid[iz_mid]:.3f})')
    ax.legend()

    # Policy function for k1_next lines
    ax = axes[1]
    for ik2, col in zip(ik2_indices, colors):
        ax.plot(mod.k1_grid, mod.k1_grid[pol_k1[:, ik2, iz_mid]], color=col,
                label=f'k2 = {mod.k2_grid[ik2]:.3f}')
    ax.set_xlabel('k1'); ax.set_ylabel("k1'")
    ax.set_title(f"Policy: k1'  (z1 = {mod.z1_grid[iz_mid]:.3f})")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'vfi_two_sector.png'), dpi=150)
    plt.show()
    print('Figure saved to vfi_two_sector.png')