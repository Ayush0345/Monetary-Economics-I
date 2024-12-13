// Replication JPT

// Declare endogenous variables (38 variables)
var y, k, L, rho, w, s, ppi, c, lambda, R, u, phi, i, kbar, x, gw, y_star, k_star, L_star, rho_star, w_star, s_star, c_star, lambda_star, R_star, u_star, phi_star, i_star, 
    kbar_star, x_star, gw_star, z, b, g, mu, etamp, lambda_p, lambda_w; 

// Declare exogenous variables (shocks)
varexo eps_z, eps_mu, eps_mp, eps_g, eps_p, eps_w, eps_b;

// Observed variables for estimation
varobs y c i;

// Declare parameters (59 parameters)
parameters aalpha, delta, ip, iw, gammam, h, lambda_p_ss, lambda_w_ss, L_ss, pi_ss, bbeta, nu, zeta_p, zeta_w, xi, s2prime, phi_pi, phi_x, phi_dx, rho_R, rho_mp, rho_z, rho_g, rho_mu, rho_p, rho_w, rho_b, theta_p, theta_w,
           theta_mp, theta_g, theta_mu, theta_b, g_ss, kappa, rho_ss, w_ss, kLratio, FLratio, F, yLratio, iLratio, cLratio, k_ss, y_ss, i_ss, c_ss, R_ss, kLstarratio, FLstarratio, Fstar, yLstarratio,
           iLstarratio, cLstarratio, kstar_ss, ystar_ss, istar_ss, cstar_ss, wstar_ss, kappa_w;

// Calibration values for parameters
delta = 0.025; 
aalpha = 0.3;
bbeta = 0.99;  
gammam = 0.005;
nu = 2;
h = 0.5;  
zeta_p = 0.84;
zeta_w = 0.7;
ip = 0.24;
lambda_w_ss = 0.15;
lambda_p_ss = 0.15;
pi_ss = 1.005; // Steady-state inflation
iw = 0.5; // Wage stickiness parameter
xi = 5; // Elasticity parameter
s2prime = 0.5; // Sensitivity parameter
rho_R = 0.8;
rho_mp = 0.4;
rho_z = 0.9;
rho_g = 0.6;
rho_mu = 0.7;
rho_p = 0.7;
rho_w = 0.7;
rho_b = 0.6;
phi_pi = 1.5;
phi_x = 0.13;
phi_dx = 0.25;
theta_p = 0.75; // Price rigidity parameter
theta_w = 0.75; // Wage rigidity parameter
L_ss = 0.3; // Initial value for steady-state labor
theta_mp = 0.10;
theta_g = 0.5;
theta_mu = 0.5;
theta_b = 0.10;

// Compute parameters that depend on calibrated values (must be before 'model' block)
g_ss = 1 / (1 - 0.22); // Fixed g_ss calculation
rho_ss = (exp(gammam)/bbeta) - (1 - delta);
w_ss = ((1 / (1 + lambda_p_ss)) * (aalpha^aalpha) * ((1 - aalpha)^(1 - aalpha)) * (1 / rho_ss^aalpha))^(1 / (1 - aalpha));
kLratio = (w_ss / rho_ss) * (aalpha / (1 - aalpha));
FLratio = kLratio^aalpha - rho_ss * kLratio - w_ss;
yLratio = kLratio^aalpha - FLratio;
iLratio = (1 - (1 - delta) * exp(-gammam)) * exp(gammam) * kLratio;     
cLratio = yLratio * (1 / g_ss) - iLratio;
lambdaL = cLratio^(-1) * (exp(gammam) - h * bbeta) / (exp(gammam) - h);
// L_ss is already assigned above in the calibration block
k_ss = kLratio * L_ss;
y_ss = yLratio * L_ss;
i_ss = iLratio * L_ss;
c_ss = cLratio * L_ss;
F = FLratio * L_ss;
oomega = ((1 / (1 + lambda_p_ss)) * (aalpha^aalpha) * ((1 - aalpha)^(1 - aalpha)) * (1 / rho_ss)^aalpha)^(1 / (1 - aalpha));
kappa = ((1 - zeta_p * bbeta) * (1 - zeta_p)) / (zeta_p * (1 + ip * bbeta));
kappa_w = ((1 - zeta_w * bbeta) * (1 - zeta_w)) / (zeta_w * (1 + bbeta) * (1 + nu * (1 + (1 / lambda_w_ss))));
R_ss = (1 + gammam) * pi_ss / bbeta;

// Steady-state values for flexible-price economy
wstar_ss = w_ss; // Assuming same as w_ss
kLstarratio = kLratio;
FLstarratio = FLratio;
yLstarratio = yLratio;
iLstarratio = iLratio;
cLstarratio = cLratio;
lambdaLstar = lambdaL;
Lstar_ss = L_ss;
kstar_ss = k_ss;
ystar_ss = y_ss;
istar_ss = i_ss;
cstar_ss = c_ss;
Fstar = F;
// Model equations in log-linearized form
model;

//1. Output equation
  y = ((y_ss + F) / y_ss) * (aalpha * k + (1 - aalpha) * L);
//2. Real wage
  rho = w + L - k;
//3. Real marginal cost
  s = aalpha * rho + (1 - aalpha) * w;
//4. Phillips curve
  ppi = (bbeta / (1 + ip * bbeta)) * ppi(+1) + (ip / (1 + ip * bbeta)) * ppi(-1) + kappa * s + kappa * lambda_p;
//5. Consumption Euler equation
  lambda = ((h * bbeta * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c(+1) 
          - ((exp(2 * gammam) + h^2 * bbeta) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c 
          + ((h * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c(-1) 
          + ((h * bbeta * exp(gammam) * rho_z - h * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * z 
          + ((exp(gammam) - h * bbeta * rho_b) / (exp(gammam) - h * bbeta)) * b;
//6. Euler equation for bonds
  lambda = R + lambda(+1) - z(+1) - ppi(+1);
//7. Capital utilization
  rho = xi * u;
//8. Value of capital
  phi = (1 - delta) * bbeta * exp(-gammam) * (phi(+1) - z(+1)) 
        + (1 - (1 - delta) * bbeta * exp(-gammam)) * (lambda(+1) - z(+1) + rho(+1));
//9. Investment Euler equation
  lambda = phi + mu - exp(2 * gammam) * s2prime * (i - i(-1) + z) 
           + bbeta * exp(2 * gammam) * s2prime * (i(+1) - i + z(+1));
//10. Capital accumulation
  k = u + kbar(-1) - z;
//11. Law of motion for capital stock
  kbar = (1 - delta) * exp(-gammam) * (kbar(-1) - z) 
         + (1 - (1 - delta) * exp(-gammam)) * (mu + i);
//12. Wage equation
  w = (1 / (1 + bbeta)) * w(-1) + (bbeta / (1 + bbeta)) * w(+1) 
      - kappa_w * gw + (iw / (1 + bbeta)) * ppi(-1) 
      - ((1 + bbeta * iw) / (1 + bbeta)) * ppi 
      + (bbeta / (1 + bbeta)) * ppi(+1) 
      + (iw / (1 + bbeta)) * z(-1) 
      - ((1 + bbeta * iw - rho_z * bbeta) / (1 + bbeta)) * z 
      + kappa_w * lambda_w;
//13. Wage gap
  gw = w - (nu * L + b - lambda);
//14. Monetary policy rule
  R = rho_R * R(-1) + (1 - rho_R) * (phi_pi * ppi + phi_x * (x - x_star)) 
      + phi_dx * ((x - x(-1)) - (x_star - x_star(-1))) + etamp;
//15. Output gap
  x = y - y_star;
//16. National income identity
  y = g + (c_ss / y_ss) * c + (i_ss / y_ss) * i + (rho_ss * k_ss / y_ss) * u;

// Equations with flexible prices and wages
//17. Flexible-price output
  y_star = ((ystar_ss + Fstar) / ystar_ss) * (aalpha * k_star + (1 - aalpha) * L_star);
//18. Flexible-price real wage
  rho_star = w_star + L_star - k_star;
//19. Flexible-price marginal cost (zero in flexible-price economy)
  s_star = aalpha * rho_star + (1 - aalpha) * w_star;
//20. Flexible-price marginal cost equals zero
  0 = s_star;
//21. Flexible-price consumption Euler equation
  lambda_star = ((h * bbeta * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c_star(+1) 
                - ((exp(2 * gammam) + h^2 * bbeta) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c_star 
                + ((h * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * c_star(-1) 
                + ((h * bbeta * exp(gammam) * rho_z - h * exp(gammam)) / (exp(gammam) - h * bbeta)) * (exp(gammam) - h) * z;
//22. Flexible-price Euler equation for bonds
  lambda_star = R_star + lambda_star(+1) - z(+1) - ppi(+1);
//23. Flexible-price capital utilization
  rho_star = xi * u_star;
//24. Flexible-price value of capital
  phi_star = (1 - delta) * bbeta * exp(-gammam) * (phi_star(+1) - z(+1)) 
             + (1 - (1 - delta) * bbeta * exp(-gammam)) * (lambda_star(+1) - z(+1) + rho_star(+1));
//25. Flexible-price investment Euler equation
  lambda_star = phi_star + mu - exp(2 * gammam) * s2prime * (i_star - i_star(-1) + z) 
                + bbeta * exp(2 * gammam) * s2prime * (i_star(+1) - i_star + z(+1));
//26. Flexible-price capital accumulation
  k_star = u_star + kbar_star(-1) - z;
//27. Flexible-price law of motion for capital stock
  kbar_star = (1 - delta) * exp(-gammam) * (kbar_star(-1) - z) 
              + (1 - (1 - delta) * exp(-gammam)) * (mu + i_star);
//28. Flexible-price wage equation
  w_star = (1 / (1 + bbeta)) * w_star(-1) + (bbeta / (1 + bbeta)) * w_star(+1) 
           - kappa_w * gw_star + (iw / (1 + bbeta)) * ppi(-1) 
           - ((1 + bbeta * iw) / (1 + bbeta)) * ppi 
           + (bbeta / (1 + bbeta)) * ppi(+1) 
           + (iw / (1 + bbeta)) * z(-1) 
           - ((1 + bbeta * iw - rho_z * bbeta) / (1 + bbeta)) * z;
//29. Flexible-price wage gap
  gw_star = w_star - (nu * L_star + b - lambda_star);
//30. Flexible-price output gap
  x_star = y_star - y_star;
//31. Flexible-price national income identity
  y_star = g + (cstar_ss / ystar_ss) * c_star + (istar_ss / ystar_ss) * i_star + (rho_ss * kstar_ss / ystar_ss) * u_star;

// Exogenous shocks
  z = rho_z * z(-1) + eps_z;             // Technology shock
  mu = rho_mu * mu(-1) + eps_mu;         // Investment-specific shock
  etamp = rho_mp * etamp(-1) + eps_mp;   // Monetary policy shock
  g = rho_g * g(-1) + eps_g;             // Government spending shock
  lambda_p = rho_p * lambda_p(-1) + eps_p; // Price markup shock
  lambda_w = rho_w * lambda_w(-1) + eps_w; // Wage markup shock
  b = rho_b * b(-1) + eps_b;             // Preference (demand) shock

end;

initval;
    y = log(y_ss);
    k = log(k_ss);
    L = log(L_ss);
    rho = log(rho_ss);
    w = log(w_ss);
    s = 0;
    ppi = 0;
    c = log(c_ss);
    lambda = log(lambdaL);
    R = log(R_ss);
    u = 0;
    phi = 0;
    i = log(i_ss);
    kbar = log(k_ss);
    x = 0;
    gw = 0;
    b = 0;
    g = 0;
    mu = 0;
    etamp = 0;
    lambda_p = 0;
    lambda_w = 0;
    z = 0;
    
    // Flexible-price variables
    y_star = log(ystar_ss);
    k_star = log(kstar_ss);
    L_star = log(Lstar_ss);
    rho_star = log(rho_ss);
    w_star = log(wstar_ss);
    s_star = 0;
    c_star = log(cstar_ss);
    lambda_star = log(lambdaLstar);
    R_star = log(R_ss);
    u_star = 0;
    phi_star = 0;
    i_star = log(istar_ss);
    kbar_star = log(kstar_ss);
    x_star = 0;
    gw_star = 0;
end;

// Shock specification (model shocks)
shocks;
    var eps_z; stderr 0.1;
    var eps_p; stderr 0.14;
    var eps_mp; stderr 0.22;
    var eps_mu; stderr 6.03;
    var eps_g; stderr 0.35;
    var eps_w; stderr 0.1;
    var eps_b; stderr 0.1;
end;

// Check residuals
resid;

// Steady-state computation
steady;

// Check Blanchard-Kahn conditions
check;

// Simulation command
//stoch_simul(irf=50, order=1);

// Estimated parameters
estimated_params;
    aalpha, beta_pdf, 0.33, 0.05; 
    bbeta, beta_pdf, 0.99, 0.002;
    nu, gamma_pdf, 2.0, 0.5;
    zeta_p, beta_pdf, 0.75, 0.05;
    zeta_w, beta_pdf, 0.75, 0.05;
    phi_pi, normal_pdf, 1.5, 0.25;
    phi_x, normal_pdf, 0.125, 0.05;
    rho_R, beta_pdf, 0.8, 0.1;
    rho_z, beta_pdf, 0.85, 0.1;
    rho_g, beta_pdf, 0.85, 0.1;
    rho_mu, beta_pdf, 0.7, 0.1;
    rho_mp, beta_pdf, 0.4, 0.1;
    rho_p, beta_pdf, 0.5, 0.1;
    rho_w, beta_pdf, 0.5, 0.1;
    rho_b, beta_pdf, 0.7, 0.1;
    h, beta_pdf, 0.7, 0.1;
    lambda_p_ss, normal_pdf, 0.15, 0.05;
    lambda_w_ss, normal_pdf, 0.15, 0.05;
    iw, beta_pdf, 0.5, 0.1;
    ip, beta_pdf, 0.25, 0.1;
    xi, gamma_pdf, 5.0, 1.0;
    s2prime, gamma_pdf, 4.0, 1.0;
    gammam, normal_pdf, 0.005, 0.003;
    L_ss, normal_pdf, 0.3, 0.1;
    pi_ss, normal_pdf, 1.005, 0.1;
    phi_dx, normal_pdf, 0.13, 0.02;
    theta_p, INV_GAMMA1_PDF, 0.5, 0.1;
    theta_w, INV_GAMMA1_PDF, 0.5, 0.1;
    theta_mp, INV_GAMMA1_PDF, 0.10, 0.1;
    theta_g, INV_GAMMA1_PDF, 0.50, 0.1;
    theta_mu, INV_GAMMA1_PDF, 0.50, 0.1;
    theta_b, INV_GAMMA1_PDF, 0.10, 0.1;
end;


// Fix parameters causing issues


// Provide initial values for estimated parameters within prior bounds
estimated_params_init;
    aalpha, 0.3;
    ip, 0.24;
    iw, 0.5;
    gammam, 0.005;
    h, 0.5;
    lambda_p_ss, 0.15;
    lambda_w_ss, 0.15;
    L_ss, 0.3;
    pi_ss, 1.005;
    bbeta, 0.99;
    nu, 2;
    zeta_p, 0.84;
    zeta_w, 0.7;
    s2prime, 4;
    phi_pi, 1.7;
    phi_x, 0.13; // Fixed to calibrated value
    xi, 5; // Fixed to calibrated value
    phi_dx, 0.13;
    rho_R, 0.6;
    rho_mp, 0.4;
    rho_z, 0.6;
    rho_g, 0.6;
    rho_mu, 0.7;
    rho_p, 0.6;
    rho_w, 0.6;
    rho_b, 0.6;
    theta_p, 0.5;
    theta_w, 0.5;
    theta_mp, 0.10;
    theta_g, 0.50;
    theta_mu, 0.50;
    theta_b, 0.10;
end;


// dynare_sensitivity
dynare_sensitivity;
identification(order=1, advanced=1);




// Estimation command
estimation(datafile='model_data.mat', prior_trunc=0, mode_compute=9, mh_replic=120000, mh_jscale=0.56, mode_file='draft_aj_mh_mode.mat');