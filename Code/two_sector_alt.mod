% This file is just an illustration of how to use Dynare's macro variable language,
% the documentation is available at https://www.dynare.org/manual/the-model-file.html#macro-processing-language
% To run this file, use the following command: dynare two_sector_alt savemacro=two_sector_mod.mod onlymacro
% The preprocessor then creates a new mod file called "two_sector_mod.mod" which looks nice and clean
% Note that the mod file here is just an illustration and incomplete;
% with a full file you don't have to use savemacro and onlymacro, but can simply run it: dynare illustration

//***
// 2-SECTOR MODEL OF STRUCTURAL CHANGE
// Refactored using Dynare Macro Language
//***

// ------------------------------------------------------------------------
// Macro Definitions
// ------------------------------------------------------------------------
@#define SECTORS = ["1", "2"]

// Grouping variables and parameters by dimensionality
@#define VARS_J = ["c", "y", "v", "m", "k", "x", "A", "py", "pv", "pm", "px", "z", "l"]
@#define VARS_IJ = ["m", "x"]

@#define PARAMS_J = ["theta", "gamma", "alpha", "delta", "rho_z", "rho_l"]
@#define PARAMS_IJ = ["phi", "omega"]

@#define SHOCKS_J = ["eps_z", "eps_l"]

// ------------------------------------------------------------------------
// Declaring variables
// ------------------------------------------------------------------------
var C
@#for j in SECTORS
  @#for v in VARS_J
    @{v}@{j}
  @#endfor
@#endfor
@#for i in SECTORS
  @#for j in SECTORS
    @#for v in VARS_IJ
      @{v}@{i}@{j}
    @#endfor
  @#endfor
@#endfor
;

varexo
@#for j in SECTORS
  @#for s in SHOCKS_J
    @{s}@{j}
  @#endfor
@#endfor
;

predetermined_variables
@#for j in SECTORS
  k@{j}
@#endfor
;

// ------------------------------------------------------------------------
// Parametrization
// ------------------------------------------------------------------------
parameters beta
@#for j in SECTORS
  @#for p in PARAMS_J
    @{p}@{j}
  @#endfor
@#endfor
@#for i in SECTORS
  @#for j in SECTORS
    @#for p in PARAMS_IJ
      @{p}@{i}@{j}
    @#endfor
  @#endfor
@#endfor
;

beta = 0.99;

@#for j in SECTORS
  theta@{j}  = 0.5;   
  gamma@{j}  = 0.5; 
  alpha@{j}  = 0.33;  
  delta@{j}  = 0.025; 
  rho_z@{j}  = 0.90;
  rho_l@{j}  = 0.90;
  
  @#for i in SECTORS
    phi@{i}@{j} = 0.5;
    omega@{i}@{j} = 0.5;
  @#endfor
@#endfor

// ------------------------------------------------------------------------
// Model block
// ------------------------------------------------------------------------
model;

// 1. Aggregate consumption (Constructed via string concatenation for scalability)
@#define C_eq = "C = "
@#for j in SECTORS
  @#if j != SECTORS[1]
    @#define C_eq = C_eq + " * "
  @#endif
  @#define C_eq = C_eq + "(c" + j + "/theta" + j + ")^theta" + j
@#endfor
@{C_eq};

@#for j in SECTORS
  // Production technology in sector j
  y@{j} = (v@{j}/gamma@{j})^gamma@{j} * (m@{j}/(1-gamma@{j}))^(1-gamma@{j});

  // Value added aggregate
  v@{j} = A@{j} * (k@{j}/alpha@{j})^alpha@{j};

  // Capital accumulation
  k@{j}(+1) = x@{j} + (1-delta@{j})*k@{j};

  // Redefined TFP 
  A@{j} = z@{j} * (l@{j} / (1-alpha@{j}))^(1-alpha@{j});

  // FOC: Consumption
  py@{j} = theta@{j} * C / c@{j};

  // FOC: Value-added aggregate
  pv@{j} = gamma@{j} * py@{j} * y@{j} / v@{j};

  // FOC: Materials aggregate
  pm@{j} = (1-gamma@{j}) * py@{j} * y@{j} / m@{j};

  // FOC: Euler equation for investment
  px@{j} = beta * ( alpha@{j} * pv@{j}(+1) * v@{j}(+1) / k@{j}(+1) + px@{j}(+1)*(1-delta@{j}) );

  // AR(1) Shock processes in log-levels
  log(z@{j}) = rho_z@{j} * log(z@{j}(-1)) + eps_z@{j};
  log(l@{j}) = rho_l@{j} * log(l@{j}(-1)) + eps_l@{j};
@#endfor

// Cross-sector aggregations and constraints
@#for j in SECTORS
  // Materials and Investment aggregates
  @#define m_eq = "m" + j + " = "
  @#define x_eq = "x" + j + " = "
  @#for i in SECTORS
    @#if i != SECTORS[1]
      @#define m_eq = m_eq + " * "
      @#define x_eq = x_eq + " * "
    @#endif
    @#define m_eq = m_eq + "(m" + i + j + "/phi" + i + j + ")^phi" + i + j
    @#define x_eq = x_eq + "(x" + i + j + "/omega" + i + j + ")^omega" + i + j
  @#endfor
  @{m_eq};
  @{x_eq};
@#endfor

@#for i in SECTORS
  // Resource constraints
  @#define res_eq = "c" + i
  @#for j in SECTORS
    @#define res_eq = res_eq + " + m" + i + j + " + x" + i + j
  @#endfor
  @#define res_eq = res_eq + " = y" + i
  @{res_eq};

  // FOCs for cross-sector purchases
  @#for j in SECTORS
    py@{i} * m@{i}@{j} = phi@{i}@{j} * pm@{j} * m@{j};
    py@{i} * x@{i}@{j} = omega@{i}@{j} * px@{j} * x@{j};
  @#endfor
@#endfor

end;

// ------------------------------------------------------------------------
// Initial Values for Steady State Solver
// ------------------------------------------------------------------------
initval;
  C = 1;
@#for j in SECTORS
  // Shock variables initialized at 1
  z@{j} = 1; 
  l@{j} = 1;
  
  // Guesses for endogenous variables
  c@{j} = 0.5; 
  y@{j} = 2;
  v@{j} = 1; 
  m@{j} = 1;
  k@{j} = 10; 
  x@{j} = 0.25;
  A@{j} = 1; 
  py@{j} = 1; 
  pv@{j} = 1;
  pm@{j} = 1; 
  px@{j} = 1;

  @#for i in SECTORS
    m@{i}@{j} = 0.5;
    x@{i}@{j} = 0.1;
  @#endfor
@#endfor
end;

steady;
check;

// ------------------------------------------------------------------------
// Shocks
// ------------------------------------------------------------------------
shocks;
@#for j in SECTORS
  var eps_z@{j} = 0.01^2; 
  var eps_l@{j} = 0.01^2;
@#endfor
end;

// ------------------------------------------------------------------------
// Impulse response
// ------------------------------------------------------------------------
stoch_simul(order=1, irf=40) c1;