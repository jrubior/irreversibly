//***
// N-SECTOR MODEL OF STRUCTURAL CHANGE
// Refactored using Dynare Macro Language
//***
// ------------------------------------------------------------------------
// Macro Definitions
// ------------------------------------------------------------------------
// Grouping variables and parameters by dimensionality
// ------------------------------------------------------------------------
// Declaring variables
// ------------------------------------------------------------------------
var C
    c1
    y1
    v1
    m1
    k1
    x1
    A1
    py1
    pv1
    pm1
    px1
    z1
    l1
    c2
    y2
    v2
    m2
    k2
    x2
    A2
    py2
    pv2
    pm2
    px2
    z2
    l2
      m11
      x11
      m12
      x12
      m21
      x21
      m22
      x22
;
varexo
    eps_z1
    eps_l1
    eps_z2
    eps_l2
;
predetermined_variables
  k1
  k2
;
// ------------------------------------------------------------------------
// Parametrization
// ------------------------------------------------------------------------
parameters beta
    theta1
    gamma1
    alpha1
    delta1
    rho_z1
    rho_l1
    theta2
    gamma2
    alpha2
    delta2
    rho_z2
    rho_l2
      phi11
      omega11
      phi12
      omega12
      phi21
      omega21
      phi22
      omega22
;
beta = 0.99;
  theta1  = 0.5;   
  gamma1  = 0.5; 
  alpha1  = 0.33;  
  delta1  = 0.025; 
  rho_z1  = 0.90;
  rho_l1  = 0.90;
  
    phi11 = 0.5;
    omega11 = 0.5;
    phi21 = 0.5;
    omega21 = 0.5;
  theta2  = 0.5;   
  gamma2  = 0.5; 
  alpha2  = 0.33;  
  delta2  = 0.025; 
  rho_z2  = 0.90;
  rho_l2  = 0.90;
  
    phi12 = 0.5;
    omega12 = 0.5;
    phi22 = 0.5;
    omega22 = 0.5;
// ------------------------------------------------------------------------
// Model block
// ------------------------------------------------------------------------
model;
// 1. Aggregate consumption (Constructed via string concatenation for scalability)
C = (c1/theta1)^theta1 * (c2/theta2)^theta2;
  // Production technology in sector j
  y1 = (v1/gamma1)^gamma1 * (m1/(1-gamma1))^(1-gamma1);
  // Value added aggregate
  v1 = A1 * (k1/alpha1)^alpha1;
  // Capital accumulation
  k1(+1) = x1 + (1-delta1)*k1;
  // Redefined TFP 
  A1 = z1 * (l1 / (1-alpha1))^(1-alpha1);
  // FOC: Consumption
  py1 = theta1 * C / c1;
  // FOC: Value-added aggregate
  pv1 = gamma1 * py1 * y1 / v1;
  // FOC: Materials aggregate
  pm1 = (1-gamma1) * py1 * y1 / m1;
  // FOC: Euler equation for investment
  px1 = beta * ( alpha1 * pv1(+1) * v1(+1) / k1(+1) + px1(+1)*(1-delta1) );
  // AR(1) Shock processes in log-levels
  log(z1) = rho_z1 * log(z1(-1)) + eps_z1;
  log(l1) = rho_l1 * log(l1(-1)) + eps_l1;
  // Production technology in sector j
  y2 = (v2/gamma2)^gamma2 * (m2/(1-gamma2))^(1-gamma2);
  // Value added aggregate
  v2 = A2 * (k2/alpha2)^alpha2;
  // Capital accumulation
  k2(+1) = x2 + (1-delta2)*k2;
  // Redefined TFP 
  A2 = z2 * (l2 / (1-alpha2))^(1-alpha2);
  // FOC: Consumption
  py2 = theta2 * C / c2;
  // FOC: Value-added aggregate
  pv2 = gamma2 * py2 * y2 / v2;
  // FOC: Materials aggregate
  pm2 = (1-gamma2) * py2 * y2 / m2;
  // FOC: Euler equation for investment
  px2 = beta * ( alpha2 * pv2(+1) * v2(+1) / k2(+1) + px2(+1)*(1-delta2) );
  // AR(1) Shock processes in log-levels
  log(z2) = rho_z2 * log(z2(-1)) + eps_z2;
  log(l2) = rho_l2 * log(l2(-1)) + eps_l2;
// Cross-sector aggregations and constraints
  // Materials and Investment aggregates
  m1 = (m11/phi11)^phi11 * (m21/phi21)^phi21;
  x1 = (x11/omega11)^omega11 * (x21/omega21)^omega21;
  // Materials and Investment aggregates
  m2 = (m12/phi12)^phi12 * (m22/phi22)^phi22;
  x2 = (x12/omega12)^omega12 * (x22/omega22)^omega22;
  // Resource constraints
  c1 + m11 + x11 + m12 + x12 = y1;
  // FOCs for cross-sector purchases
    py1 * m11 = phi11 * pm1 * m1;
    py1 * x11 = omega11 * px1 * x1;
    py1 * m12 = phi12 * pm2 * m2;
    py1 * x12 = omega12 * px2 * x2;
  // Resource constraints
  c2 + m21 + x21 + m22 + x22 = y2;
  // FOCs for cross-sector purchases
    py2 * m21 = phi21 * pm1 * m1;
    py2 * x21 = omega21 * px1 * x1;
    py2 * m22 = phi22 * pm2 * m2;
    py2 * x22 = omega22 * px2 * x2;
end;
// ------------------------------------------------------------------------
// Initial Values for Steady State Solver
// ------------------------------------------------------------------------
initval;
  C = 1;
  // Shock variables initialized at 1
  z1 = 1; 
  l1 = 1;
  
  // Guesses for endogenous variables
  c1 = 0.5; 
  y1 = 2;
  v1 = 1; 
  m1 = 1;
  k1 = 10; 
  x1 = 0.25;
  A1 = 1; 
  py1 = 1; 
  pv1 = 1;
  pm1 = 1; 
  px1 = 1;
    m11 = 0.5;
    x11 = 0.1;
    m21 = 0.5;
    x21 = 0.1;
  // Shock variables initialized at 1
  z2 = 1; 
  l2 = 1;
  
  // Guesses for endogenous variables
  c2 = 0.5; 
  y2 = 2;
  v2 = 1; 
  m2 = 1;
  k2 = 10; 
  x2 = 0.25;
  A2 = 1; 
  py2 = 1; 
  pv2 = 1;
  pm2 = 1; 
  px2 = 1;
    m12 = 0.5;
    x12 = 0.1;
    m22 = 0.5;
    x22 = 0.1;
end;
steady;
check;
// ------------------------------------------------------------------------
// Shocks
// ------------------------------------------------------------------------
shocks;
  var eps_z1 = 0.01^2; 
  var eps_l1 = 0.01^2;
  var eps_z2 = 0.01^2; 
  var eps_l2 = 0.01^2;
end;
// ------------------------------------------------------------------------
// Impulse response
// ------------------------------------------------------------------------
stoch_simul(order=1, irf=40) c1;
