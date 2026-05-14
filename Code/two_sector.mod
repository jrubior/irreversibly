//***
//TWO SECTOR MODEL OF STRUCTURAL CHANGE
//JUAN RUBIO-RAMIREZ, JAFET BACA
//***

// ------------------------------------------------------------------------
// Declaring variables
// ------------------------------------------------------------------------

var C c1 c2 y1 y2 v1 v2 m1 m2 m11 m12 m21 m22 k1 k2 x1 x2 x11 x12 x21 x22
    A1 A2 py1 py2 pv1 pv2 pm1 pm2 px1 px2 z1 z2 l1 l2;
varexo eps_z1 eps_z2 eps_l1 eps_l2;
predetermined_variables k1 k2;

// ------------------------------------------------------------------------
// Parametrization
// ------------------------------------------------------------------------

parameters beta theta1 theta2 gamma1 gamma2 alpha1 alpha2 delta1 delta2 
           phi11 phi21 phi12 phi22 omega11 omega21 omega12 omega22
           rho_z1 rho_z2 rho_l1 rho_l2;

beta    = 0.99;
theta1  = 0.5;   theta2  = 0.5;   // Must sum to 1
gamma1  = 0.5;   gamma2  = 0.5; 
alpha1  = 0.33;  alpha2  = 0.33;
delta1  = 0.025; delta2  = 0.025;

// Materials parameters (sum_i phi_ij = 1)
phi11 = 0.5; phi21 = 0.5; 
phi12 = 0.5; phi22 = 0.5; 

// Investment parameters (sum_i omega_ij = 1)
omega11 = 0.5; omega21 = 0.5;
omega12 = 0.5; omega22 = 0.5;

// Persistence parameters
rho_z1 = 0.90; rho_z2 = 0.90;
rho_l1 = 0.90; rho_l2 = 0.90;

// ------------------------------------------------------------------------
// Model block
// ------------------------------------------------------------------------

model;
    // 1. Aggregate consumption
    exp(C) = (exp(c1)/theta1)^theta1 * (exp(c2)/theta2)^theta2;

    // 2-3. Production technology
    exp(y1) = (exp(v1)/gamma1)^gamma1 * (exp(m1)/(1-gamma1))^(1-gamma1);
    exp(y2) = (exp(v2)/gamma2)^gamma2 * (exp(m2)/(1-gamma2))^(1-gamma2); 

    // 4-5. Materials aggregate
    exp(m1) = (exp(m11)/phi11)^phi11 * (exp(m21)/phi21)^phi21;
    exp(m2) = (exp(m12)/phi12)^phi12 * (exp(m22)/phi22)^phi22;

    // 6-7. Value added aggregate
    exp(v1) = exp(A1) * (exp(k1)/alpha1)^alpha1;
    exp(v2) = exp(A2) * (exp(k2)/alpha2)^alpha2;

    // 8-9. Capital accumulation
    exp(k1(+1)) = exp(x1) + (1-delta1)*exp(k1);
    exp(k2(+1)) = exp(x2) + (1-delta2)*exp(k2);

    // 10-11. New investment aggregate
    exp(x1) = (exp(x11)/omega11)^omega11 * (exp(x21)/omega21)^omega21;
    exp(x2) = (exp(x12)/omega12)^omega12 * (exp(x22)/omega22)^omega22;

    // 12-13. Resource constraints
    exp(c1) + exp(m11) + exp(m12) + exp(x11) + exp(x12) = exp(y1);
    exp(c2) + exp(m21) + exp(m22) + exp(x21) + exp(x22) = exp(y2); 

    // 14-15. Redefined TFP 
    exp(A1) = exp(z1) * (exp(l1) / (1-alpha1))^(1-alpha1);
    exp(A2) = exp(z2) * (exp(l2) / (1-alpha2))^(1-alpha2);

    // 16-17. FOC: Consumption
    exp(py1) = theta1 * exp(C) / exp(c1); 
    exp(py2) = theta2 * exp(C) / exp(c2);

    // 18-19. FOC: Value-added aggregate
    exp(pv1) = gamma1 * exp(py1) * exp(y1) / exp(v1); 
    exp(pv2) = gamma2 * exp(py2) * exp(y2) / exp(v2); 

    // 20-21. FOC: Materials aggregate
    exp(pm1) = (1-gamma1) * exp(py1) * exp(y1) / exp(m1); 
    exp(pm2) = (1-gamma2) * exp(py2) * exp(y2) / exp(m2); 

    // 22-25. FOC: Materials bought by j from i
    exp(py1) * exp(m11) = phi11 * exp(pm1) * exp(m1); 
    exp(py2) * exp(m21) = phi21 * exp(pm1) * exp(m1);
    exp(py1) * exp(m12) = phi12 * exp(pm2) * exp(m2);
    exp(py2) * exp(m22) = phi22 * exp(pm2) * exp(m2);

    // 26-29. FOC: New capital bought by j from i
    exp(py1) * exp(x11) = omega11 * exp(px1) * exp(x1); 
    exp(py2) * exp(x21) = omega21 * exp(px1) * exp(x1); 
    exp(py1) * exp(x12) = omega12 * exp(px2) * exp(x2); 
    exp(py2) * exp(x22) = omega22 * exp(px2) * exp(x2); 

    // 30-31. FOC: Euler equation for investment
    exp(px1) = beta * ( alpha1 * exp(pv1(+1)) * exp(v1(+1)) / exp(k1(+1)) + exp(px1(+1))*(1-delta1) ); 
    exp(px2) = beta * ( alpha2 * exp(pv2(+1)) * exp(v2(+1)) / exp(k2(+1)) + exp(px2(+1))*(1-delta2) );

    // 32-35. AR(1) Shock processes (already in log-levels)
    z1 = rho_z1 * z1(-1) + eps_z1; 
    z2 = rho_z2 * z2(-1) + eps_z2; 
    l1 = rho_l1 * l1(-1) + eps_l1; 
    l2 = rho_l2 * l2(-1) + eps_l2; 
end;

// ------------------------------------------------------------------------
// Initial Values for Steady State Solver
// ------------------------------------------------------------------------

steady_state_model;
    // 1. Exogenous variables and TFP at steady state
    z1 = 0; z2 = 0; l1 = 0; l2 = 0;
    // log(A) from Eq 14-15
    A1 = z1 + (1-alpha1)*(l1 - log(1-alpha1));
    A2 = z2 + (1-alpha2)*(l2 - log(1-alpha2));

    // 2. Relative Prices (Normalized to 1 in logs for the symmetric case)
    py1 = 0; py2 = 0; pm1 = 0; pm2 = 0; px1 = 0; px2 = 0; pv1 = 0; pv2 = 0;

    // 3. Solve for Real Levels using Euler and Production Eq
    // From Euler: k/v = (alpha * beta) / (1 - beta*(1-delta))
    k_v_ratio1 = (alpha1 * beta) / (1 - beta*(1-delta1));
    k_v_ratio2 = (alpha2 * beta) / (1 - beta*(1-delta2));

    // From Production v = A * (k/alpha)^alpha
    // log(v) = A + alpha*log(k) - alpha*log(alpha)
    // Substituting k = v * k_v_ratio:
    v1 = (A1 + alpha1*log(k_v_ratio1) - alpha1*log(alpha1)) / (1-alpha1);
    v2 = (A2 + alpha2*log(k_v_ratio2) - alpha2*log(alpha2)) / (1-alpha2);
    
    k1 = v1 + log(k_v_ratio1);
    k2 = v2 + log(k_v_ratio2);

    // 4. Solve for Output and Aggregates
    // From FOC VA: v = gamma * y / pv -> log(y) = log(v) - log(gamma) + log(pv)
    y1 = v1 - log(gamma1);
    y2 = v2 - log(gamma2);

    // From FOC Materials: m = (1-gamma) * y
    m1 = log(1-gamma1) + y1;
    m2 = log(1-gamma2) + y2;

    // 5. Solve for Investment Levels
    // In SS: x = delta * k
    x1 = log(delta1) + k1;
    x2 = log(delta2) + k2;

    // 6. Decompose Aggregates into Sectoral Flows
    // Materials: py_i * m_ij = phi_ij * pm_j * m_j
    m11 = log(phi11) + m1;
    m21 = log(phi21) + m1;
    m12 = log(phi12) + m2;
    m22 = log(phi22) + m2;

    // Investment: py_i * x_ij = omega_ij * px_j * x_j
    x11 = log(omega11) + x1;
    x21 = log(omega21) + x1;
    x12 = log(omega12) + x2;
    x22 = log(omega22) + x2;

    // 7. Solve for Consumption from Resource Constraint
    // exp(c1) = exp(y1) - exp(m11) - exp(m12) - exp(x11) - exp(x12)
    c1 = log(exp(y1) - exp(m11) - exp(m12) - exp(x11) - exp(x12));
    c2 = log(exp(y2) - exp(m21) - exp(m22) - exp(x21) - exp(x22));

    // 8. Aggregate Consumption
    C = theta1*(c1 - log(theta1)) + theta2*(c2 - log(theta2));
end;

steady;
check;

// ------------------------------------------------------------------------
// Shocks
// ------------------------------------------------------------------------

shocks;
var eps_z1 = 0.01^2; 
var eps_z2 = 0.01^2;
var eps_l1 = 0.01^2;
var eps_l2 = 0.01^2;
end;

// ------------------------------------------------------------------------
// Impulse response
// ------------------------------------------------------------------------

stoch_simul(order=1, irf=40) k1;