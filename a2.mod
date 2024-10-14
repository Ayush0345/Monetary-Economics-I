var c k inv n r w a y infl m;

varexo eps_a eps_m;   // Shocks: technology (eps_a), monetary policy (eps_m)

parameters beta gamma phi delta alpha lambda xi_w xi_p rho_i zeta_pi zeta_y rho_a sigma_a rho_m sigma_m omega eta_w psi s_prime_prime;

beta     = 0.99;      // Discount factor
gamma    = 1;         // Risk aversion parameter
phi      = 1.25;      // Inverse elasticity of labor supply
delta    = 0.025;     // Depreciation rate
alpha    = 0.35;      // Capital share in production
lambda   = 0.75;      // Habit formation parameter
xi_w     = 0.9;       // Wage stickiness
xi_p     = 0.6;       // Price stickiness
rho_i    = 0.78;      // Taylor rule inertia
zeta_pi  = 1.5;       // Taylor rule coefficient on inflation
zeta_y   = 0.1;       // Taylor rule coefficient on output gap
rho_a    = 0.92;      // AR(1) for technology shock
sigma_a  = 0.04;      // Std dev for technology shock
rho_m    = 0.4;       // AR(1) for monetary policy shock
sigma_m  = 0.1;       // Std dev for monetary policy shock
eta_w    = 2.1;       // Elasticity of substitution for labor
psi      = 0.6;       // Capital utilization adjustment
s_prime_prime = 5;    // Sensitivity of investment adjustment

// Steady-state value of labor consistent with 38-hour work week
     n_steady = 38/168;  // 38 hours out of 168 hours in a week

// Assume technology level a = 1 in steady state
     a_steady = 1;       

// Find steady-state capital k_steady
     k_steady = ((n_steady^(1-alpha)) * a_steady / ((1/beta - (1-delta)) / alpha))^(1/(alpha-1));

// Steady-state output
     y_steady = a_steady * k_steady^alpha * n_steady^(1 - alpha);

// Steady-state consumption
     c_steady = y_steady - delta * k_steady;

// Model Local-Variables

// Calibrate Psi using the steady-state labor supply equation
     omega_value = ((1 - alpha) * y_steady) / (c_steady * n_steady^phi); 

      omega = omega_value;  // Assign calibrated Psi

model;
// Euler equation (intertemporal consumption choice)
c = beta * c(+1) + (1 - beta) * r - lambda * (c - c(-1));

// Capital accumulation (law of motion for capital)
k = (1 - delta) * k(-1) + inv - s_prime_prime * (inv - inv(-1));

// Corrected labor supply condition (labor market clearing)
w = phi * n - omega * (c - c(-1));

// Goods market equilibrium (aggregate demand = aggregate supply)
y = c + inv;
y = a + alpha * k(-1) + (1 - alpha) * n;

// Taylor rule (interest rate setting by the central bank)
r = rho_i * r(-1) + (1 - rho_i) * (zeta_pi * infl + zeta_y * (y - y(-1))) + eps_m;

// AR(1) process for technology shock
a = rho_a * a(-1) + eps_a;

// AR(1) process for monetary policy shock
m = rho_m * m(-1) + eps_m;

// Inflation dynamic equation
infl = beta * infl(+1) + xi_p * y;

// Wage inflation dynamic equation
w = beta * w(+1) + xi_w * n;

end;

// Shocks and steady state
shocks;
var eps_a = sigma_a^2;  // Variance of technology shock
var eps_m = sigma_m^2;  // Variance of monetary policy shock
end;

// Initial steady state values
initval;
c = 0; k = 0; inv = 0; n = 0; r = 0; w = 0; a = 0; y = 0; infl = 0; m = 0;
end;

steady;

// Generate impulse response functions (IRFs)
stoch_simul(irf=50, order=1, noprint);

figure('Name', 'IRFs to Monetary Policy Shock');
tiledlayout(3, 3); // Create a 3x3 grid for subplots

% Steady-state value (assumed to be zero in log-linear models)
ss_value = 0;  

// Plot Impulse Response Functions using the new layout with steady-state lines
nexttile;
plot(oo_.irfs.y_eps_m, 'r', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Output (y)');
hold off;

nexttile;
plot(oo_.irfs.c_eps_m, 'r', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Consumption (c)');
hold off;

nexttile;
plot(oo_.irfs.inv_eps_m, 'b', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Investment (inv)');
hold off;

nexttile;
plot(oo_.irfs.n_eps_m, 'b', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Labor (n)');
hold off;

nexttile;
plot(oo_.irfs.r_eps_m, 'g', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Return on Capital (r)');
hold off;

nexttile;
plot(oo_.irfs.w_eps_m, 'g', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Wages (w)');
hold off;

nexttile;
plot(oo_.irfs.a_eps_m, 'k', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Technology (a)');
hold off;

// IRFs to Technology Shock

% Second set of plots: IRFs to Technology Shock
figure('Name', 'IRFs to Technology Shock');
tiledlayout(3, 3); % Create a 3x3 grid for subplots

nexttile;
plot(oo_.irfs.y_eps_a, 'r', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Output (y)');
hold off;

nexttile;
plot(oo_.irfs.c_eps_a, 'r', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Consumption (c)');
hold off;

nexttile;
plot(oo_.irfs.inv_eps_a, 'b', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Investment (inv)');
hold off;

nexttile;
plot(oo_.irfs.n_eps_a, 'b', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Labor (n)');
hold off;

nexttile;
plot(oo_.irfs.r_eps_a, 'g', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Return on Capital (r)');
hold off;

nexttile;
plot(oo_.irfs.w_eps_a, 'g', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Wages (w)');
hold off;

nexttile;
plot(oo_.irfs.a_eps_a, 'k', 'LineWidth', 2); 
hold on; 
yline(ss_value, '--k', 'LineWidth', 1.5); % Add steady-state line
title('Technology (a)');
hold off;

// NK v/s RBC

// Store NK model IRFs for technology shock
irf_y_nk = oo_.irfs.y_eps_a;
irf_c_nk = oo_.irfs.c_eps_a;
irf_inv_nk = oo_.irfs.inv_eps_a;
irf_n_nk = oo_.irfs.n_eps_a;
irf_r_nk = oo_.irfs.r_eps_a;
irf_w_nk = oo_.irfs.w_eps_a;
irf_a_nk = oo_.irfs.a_eps_a;

// Switch to RBC model (shut off price and wage stickiness and monopolistic competition)
xi_w = 0;  // Turn off wage stickiness
xi_p = 0;  // Turn off price stickiness
eta_w = Inf;  // Perfect competition in labor market

// Generate impulse response functions (IRFs) for RBC model
stoch_simul(irf=50, order=1, noprint);

// Store RBC model IRFs for technology shock
irf_y_rbc = oo_.irfs.y_eps_a;
irf_c_rbc = oo_.irfs.c_eps_a;
irf_inv_rbc = oo_.irfs.inv_eps_a;
irf_n_rbc = oo_.irfs.n_eps_a;
irf_r_rbc = oo_.irfs.r_eps_a;
irf_w_rbc = oo_.irfs.w_eps_a;
irf_a_rbc = oo_.irfs.a_eps_a;

// Create a figure comparing IRFs to technology shocks for NK and RBC models
figure('Name', 'NK vs RBC IRFs to Technology Shock');
tiledlayout(3, 2);  // Create a 3x2 grid for comparison plots

% Plot Output (y)
nexttile;
plot(irf_y_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_y_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Output (y)');
legend('NK', 'RBC');
hold off;

% Plot Consumption (c)
nexttile;
plot(irf_c_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_c_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Consumption (c)');
legend('NK', 'RBC');
hold off;

% Plot Investment (inv)
nexttile;
plot(irf_inv_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_inv_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Investment (inv)');
legend('NK', 'RBC');
hold off;

% Plot Labor (n)
nexttile;
plot(irf_n_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_n_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Labor (n)');
legend('NK', 'RBC');
hold off;

% Plot Return on Capital (r)
nexttile;
plot(irf_r_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_r_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Return on Capital (r)');
legend('NK', 'RBC');
hold off;

% Plot Wages (w)
nexttile;
plot(irf_w_nk, 'r', 'LineWidth', 1.5); 
hold on;
plot(irf_w_rbc, 'b--', 'LineWidth', 1);  % RBC response
title('Wages (w)');
legend('NK', 'RBC');
hold off;