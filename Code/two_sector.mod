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
    C = (c1/theta1)^theta1 * (c2/theta2)^theta2;

    // 2-3. Production technology in sector j
    y1 = (v1/gamma1)^gamma1 * (m1/(1-gamma1))^(1-gamma1);
    y2 = (v2/gamma2)^gamma2 * (m2/(1-gamma2))^(1-gamma2);

    // 4-5. Materials aggregate
    m1 = (m11/phi11)^phi11 * (m21/phi21)^phi21;
    m2 = (m12/phi12)^phi12 * (m22/phi22)^phi22;

    // 6-7. Value added aggregate
    v1 = A1 * (k1/alpha1)^alpha1;
    v2 = A2 * (k2/alpha2)^alpha2;

    // 8-9. Capital accumulation
    k1(+1) = x1 + (1-delta1)*k1;
    k2(+1) = x2 + (1-delta2)*k2;

    // 10-11. New investment aggregate
    x1 = (x11/omega11)^omega11 * (x21/omega21)^omega21;
    x2 = (x12/omega12)^omega12 * (x22/omega22)^omega22;

    // 12-13. Resource constraints
    c1 + m11 + m12 + x11 + x12 = y1;
    c2 + m21 + m22 + x21 + x22 = y2;

    // 14-15. Redefined TFP 
    A1 = z1 * (l1 / (1-alpha1))^(1-alpha1);
    A2 = z2 * (l2 / (1-alpha2))^(1-alpha2);

    // 16-17. FOC: Consumption
    py1 = theta1 * C / c1;
    py2 = theta2 * C / c2;

    // 18-19. FOC: Value-added aggregate
    pv1 = gamma1 * py1 * y1 / v1;
    pv2 = gamma2 * py2 * y2 / v2;

    // 20-21. FOC: Materials aggregate
    pm1 = (1-gamma1) * py1 * y1 / m1;
    pm2 = (1-gamma2) * py2 * y2 / m2;

    // 22-25. FOC: Materials bought by j from i
    py1 * m11 = phi11 * pm1 * m1;
    py2 * m21 = phi21 * pm1 * m1;
    py1 * m12 = phi12 * pm2 * m2;
    py2 * m22 = phi22 * pm2 * m2;

    // 26-29. FOC: New capital bought by j from i
    py1 * x11 = omega11 * px1 * x1;
    py2 * x21 = omega21 * px1 * x1;
    py1 * x12 = omega12 * px2 * x2;
    py2 * x22 = omega22 * px2 * x2;

    // 30-31. FOC: Euler equation for investment
    px1 = beta * ( alpha1 * pv1(+1) * v1(+1) / k1(+1) + px1(+1)*(1-delta1) );
    px2 = beta * ( alpha2 * pv2(+1) * v2(+1) / k2(+1) + px2(+1)*(1-delta2) );

    // 32-35. AR(1) Shock processes in log-levels
    log(z1) = rho_z1 * log(z1(-1)) + eps_z1;
    log(z2) = rho_z2 * log(z2(-1)) + eps_z2;
    log(l1) = rho_l1 * log(l1(-1)) + eps_l1;
    log(l2) = rho_l2 * log(l2(-1)) + eps_l2;
end;

// ------------------------------------------------------------------------
// Initial Values for Steady State Solver
// ------------------------------------------------------------------------

initval;
    // Shock variables initialized at 1
    z1 = 1; z2 = 1; l1 = 1; l2 = 1;
    
    // Guesses for endogenous variables
    C = 1; c1 = 0.5; c2 = 0.5;
    y1 = 2; y2 = 2;
    v1 = 1; v2 = 1;
    m1 = 1; m2 = 1;
    m11 = 0.5; m12 = 0.5; m21 = 0.5; m22 = 0.5;
    k1 = 10; k2 = 10;
    x1 = 0.25; x2 = 0.25;
    x11 = 0.1; x12 = 0.1; x21 = 0.1; x22 = 0.1;
    A1 = 1; A2 = 1;
    py1 = 1; py2 = 1; pv1 = 1; pv2 = 1; pm1 = 1; pm2 = 1; px1 = 1; px2 = 1;
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

stoch_simul(order=1, irf=40) c1;
