// Assignment 1

// Part 1

var c k n y inv r w a pi;  // Endogenous variables
varexo e;  // Exogenous variables (shock)

// Parameters
parameters beta rho alpha delta phi gamma sigma Psi;

// Calibration of Parameters
beta = 0.99;
rho = 0.92;
alpha = 0.35;
delta = 0.025;
phi = 2;
gamma = 1;
sigma = 0.04;

// Variables Definition

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
     Psi_value = ((1 - alpha) * y_steady) / (c_steady * n_steady^phi); 

      Psi = Psi_value;  // Assign calibrated Psi

// Model equations
model;

    // Euler equation
    (1/c) = beta * (1/c(+1)) * (1 + r(+1) - delta);
    
    // Labor supply equation
    Psi * n^phi = (1 - alpha) * y / c;
    
    // Production function
    y = a * k^alpha * n^(1 - alpha);
    
    // Capital accumulation equation
    k = (1 - delta) * k(-1) + inv;
    
    // Resource constraint
    y = c + inv;
    
    // Wage equation
    w = (1 - alpha) * y / n;
    
    // Return on capital
    r = alpha * y / k;
    
    // Profit equation
    pi = y - w * n - r * k;
    
    // Technology shock process
    log(a) = rho * log(a(-1)) + e;
end;

// Initial values
initval;
    c = c_steady;
    k = k_steady;
    n = n_steady; 
    y = y_steady;
    inv = delta * k_steady;
    r = alpha * y_steady / k_steady;
    w = (1 - alpha) * y_steady / n_steady;
    pi = 0;
    a = 1;
    e = 0;
end;

// Shocks
shocks;
    var e; stderr sigma;
end;

// Steady-state computation
steady;

// Impulse response functions (IRFs)
stoch_simul(order=1, irf=50, periods=200);


figure(1);

subplot(3, 2, 1); plot(oo_.irfs.y_e, 'r--', 'LineWidth', 2); title('Output (y)');
subplot(3, 2, 2); plot(oo_.irfs.c_e, 'r--', 'LineWidth', 2); title('Consumption (c)');
subplot(3, 2, 3); plot(oo_.irfs.inv_e, 'b--', 'LineWidth', 2); title('Investment (inv)');
subplot(3, 2, 4); plot(oo_.irfs.n_e, 'b--', 'LineWidth', 2); title('Labor (n)');
subplot(3, 2, 5); plot(oo_.irfs.r_e, 'g--', 'LineWidth', 2); title('Return on Capital (r)');
subplot(3, 2, 6); plot(oo_.irfs.w_e, 'g--', 'LineWidth', 2); title('Wages (w)');

// Part 4

// Autocorrelation for all variables declared

figure(2);

// Autocorrelation for each variable at 6 lags

subplot(3, 2, 1); autocorr(oo_.endo_simul(strmatch('y', M_.endo_names, 'exact'),:), 6); title('Output (y)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'r', 'LineStyle', '--');

subplot(3, 2, 2); autocorr(oo_.endo_simul(strmatch('c', M_.endo_names, 'exact'),:), 6); title('Consumption (c)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'r', 'LineStyle', '--');

subplot(3, 2, 3); autocorr(oo_.endo_simul(strmatch('inv', M_.endo_names, 'exact'),:), 6); title('Investment (inv)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'b', 'LineStyle', '--');

subplot(3, 2, 4); autocorr(oo_.endo_simul(strmatch('n', M_.endo_names, 'exact'),:), 6); title('Labor (n)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'b', 'LineStyle', '--');

subplot(3, 2, 5); autocorr(oo_.endo_simul(strmatch('r', M_.endo_names, 'exact'),:), 6); title('Return on Capital (r)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'g', 'LineStyle', '--');

subplot(3, 2, 6); autocorr(oo_.endo_simul(strmatch('w', M_.endo_names, 'exact'),:), 6); title('Wages (w)');
h = findobj(gca, 'Type', 'Line');
set(h, 'Color', 'g', 'LineStyle', '--');

// END

