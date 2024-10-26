%% Assignment 3 - Monetary Economics

%% Loading all data variables

% 1. Output

file1 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\GDPC1 (1).xls";
rgdp = readmatrix(file1);

output = rgdp(:,2);

lag_output = circshift(output, 1);

output_growth = output - lag_output ./ lag_output;

% 2. Consumption

file2 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\PCEND.xls";

c = readmatrix(file2);

consump = c(:,2);

lag_consump = circshift(consump, 1);

consump_growth = consump - lag_consump ./ consump;

% 3. Investment

file3 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\PCEDG.xls";
file4 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\GPDI (1).xls";

durable = readmatrix(file3);
private_inv = readmatrix(file4);

dg = durable(:,2);
dgs = rmmissing(dg);
pinv = private_inv(:,2);
pinvv = rmmissing(pinv);


investment = dgs + pinvv;

lag_inv = circshift(investment, 1);

inv_growth = investment - lag_inv ./ lag_inv;

% 4. Wages

file5 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\COMPRNFB.xls";

w = readmatrix(file5);

wages = w(:,2);

lag_wages = circshift(wages, 1);

wage_growth = wages - lag_wages ./ lag_wages;

% 5. Labor 

file6 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\PRS85006032.xls";
file66 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\POP.xls";

n = readmatrix(file6);
pp = readmatrix(file66);

labor = n(:,2);
pop = pp(:,2);

labor_input = labor ./ pop;

% 6. Inflation

file7 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\USAGDPDEFQISMEI.xls";
inf = readmatrix(file7);

infl = inf(:,2);

lag_infl = circshift(infl, 1);

inff = infl - lag_infl ./ lag_infl;

inflation = log(inff);

% Interest Rate

file8 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\EFFR.xls";

int = readmatrix(file8);

ffr = int(:,2);

%% Demeaning all data variables

demeaned_output = output_growth - mean(output_growth);

demeaned_consumption = consump_growth - mean(consump_growth);

demeaned_investment = inv_growth - mean(inv_growth);

demeaned_wages = wage_growth - mean(wage_growth);

demeaned_labor = labor_input - mean(labor_input);

demeaned_inf = inflation - mean(inflation);

demeaned_ir = ffr - mean(ffr);

%% Detrending all data variables

[o_trend, o_cycle] = hpfilter(output_growth, 1600);

[c_trend, c_cycle] = hpfilter(consump_growth, 1600);

[n_trend, n_cycle] = hpfilter(labor_input, 1600);

[inv_trend, inv_cycle] = hpfilter(inv_growth, 1600);

[w_trend, w_cycle] = hpfilter(wage_growth, 1600);

[inf_trend, inf_cycle] = hpfilter(inflation, 1600);

[i_trend, i_cycle] = hpfilter(ffr, 1600);

%% Plotting Demeaned and Detrended variables with NBER declared recession periods

file9 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\GDPC1 (1).xls";

dd = readtable(file9);
date = dd.Date;
dates = datetime(date);

file10 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\EFFR.xls";

ffrdd = readtable(file10);
ffr_date = ffrdd.Date;
ffr_dates = datetime(ffr_date);

% Define NBER Recession Periods (Start and End Dates)

recession_periods = [
    datetime(1960,4,1), datetime(1961,1,1);
    datetime(1970,1,1), datetime(1970,10,1);
    datetime(1974,1,1), datetime(1975,4,1);
    datetime(1980,1,1), datetime(1980,7,1);
    datetime(1981,7,1), datetime(1982,10,1);
    datetime(2001,4,1), datetime(2001,10,1);
    datetime(2008,1,1), datetime(2009,7,1);
    datetime(2020,1,1), datetime(2020,4,1)
];

recession_periods_ffr = [
    datetime(2001,4,1), datetime(2001,10,1);
    datetime(2008,1,1), datetime(2009,7,1);
    datetime(2020,1,1), datetime(2020,4,1)
];

figure(1);

% 1. Output Growth: Left (Demeaned), Right (Cycle)
subplot(7, 2, 1);
hold on;
plot(dates, demeaned_output, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned Output Growth');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 2);
hold on;
plot(dates, o_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Output Growth)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 2. Consumption Growth: Left (Demeaned), Right (Cycle)
subplot(7, 2, 3);
hold on;
plot(dates, demeaned_consumption, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned Consumption Growth');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 4);
hold on;
plot(dates, c_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Consumption Growth)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 3. Investment Growth
subplot(7, 2, 5);
hold on;
plot(dates, demeaned_investment, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned Investment Growth');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 6);
hold on;
plot(dates, inv_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Investment Growth)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 4. Wage Growth

subplot(7, 2, 7);
hold on;
plot(dates, demeaned_wages, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Wage Growth)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 8);
hold on;
plot(dates, w_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Wage Growth)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 5. Labor Hours

subplot(7, 2, 9);
hold on;
plot(dates, demeaned_labor, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 10);
hold on;
plot(dates, n_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 6. Inflation

subplot(7, 2, 11);
hold on;
plot(dates(1:255,:), demeaned_inf, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Inflation)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 12);
hold on;
plot(dates(1:255,:), inf_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods, 1)
    fill([recession_periods(j, 1), recession_periods(j, 2), ...
          recession_periods(j, 2), recession_periods(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Inflation)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

% 7. Nominal Interest Rate

subplot(7, 2, 13);
hold on;
plot(ffr_dates, demeaned_ir, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods_ffr, 1)
    fill([recession_periods_ffr(j, 1), recession_periods_ffr(j, 2), ...
          recession_periods_ffr(j, 2), recession_periods_ffr(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Nominal Interest Rate)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(7, 2, 14);
hold on;
plot(ffr_dates, i_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(recession_periods_ffr, 1)
    fill([recession_periods_ffr(j, 1), recession_periods_ffr(j, 2), ...
          recession_periods_ffr(j, 2), recession_periods_ffr(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Nominal Interest Rate)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

%% Alternate Time Series

% Inflation 

% I use Trimmed Mean PCE Inflation Rate (measure of core inflation) as my
% alternate time series for Inflation

file11 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\PCETRIM12M159SFRBDAL.xls";

pce = readmatrix(file11);

alt_inf = pce(:,2);
    
[alt_inf_trend, alt_inf_cycle] = hpfilter(alt_inf, 1600);

% Labor Hours

% I use Average Weekly Hours of All Employees, Total Private as my
% alternate time series for measuring Labor Hours

file12 = "C:\Users\ayush\OneDrive\Desktop\Economics PhD\6.1\JPT Data and Codes\AWHAETP.xls";

ahw = readmatrix(file12);

alt_labor = ahw(:,2);

lag_labor = circshift(alt_labor, 1);

alt_labor_growth = alt_labor - lag_labor ./ lag_labor;

demeaned_alt_labor = alt_labor_growth - mean(alt_labor_growth);

[l_trend, l_cycle] = hpfilter(alt_labor_growth, 1600);


dpce = readtable(file11);
inf_date = dpce.Date;
inf_dates = datetime(inf_date);

alb = readtable(file12);
labor_date = alb.Date;
labor_dates = datetime(labor_date);

% Plotting 
inf_rp = [
    datetime(1980,1,1), datetime(1980,7,1);
    datetime(1981,7,1), datetime(1982,10,1);
    datetime(2001,4,1), datetime(2001,10,1);
    datetime(2008,1,1), datetime(2009,7,1);
    datetime(2020,1,1), datetime(2020,4,1)
];


labor_rp = [
    datetime(2008,1,1), datetime(2009,7,1);
    datetime(2020,1,1), datetime(2020,4,1)
];

figure(2);
subplot(4, 2, 1);
hold on;
plot(dates(73:255,:), demeaned_inf(73:255,:), 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(inf_rp, 1)
    fill([inf_rp(j, 1), inf_rp(j, 2), ...
          inf_rp(j, 2), inf_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Inflation)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 2);
hold on;
plot(dates(73:255,:), inf_cycle(73:255,:), 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(inf_rp, 1)
    fill([inf_rp(j, 1), inf_rp(j, 2), ...
          inf_rp(j, 2), inf_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Inflation)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 3);
hold on;
plot(inf_dates(1:186,:), alt_inf(1:186,:), 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(inf_rp, 1)
    fill([inf_rp(j, 1), inf_rp(j, 2), ...
          inf_rp(j, 2), inf_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Trimmed Mean PCE Inflation Rate');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 4);
hold on;
plot(inf_dates(1:186,:), alt_inf_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(inf_rp, 1)
    fill([inf_rp(j, 1), inf_rp(j, 2), ...
          inf_rp(j, 2), inf_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Trimmed Mean PCE Inflation Rate)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 5);
hold on;
plot(dates(186:258,:), demeaned_labor(186:258,:), 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(labor_rp, 1)
    fill([labor_rp(j, 1), labor_rp(j, 2), ...
          labor_rp(j, 2), labor_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 6);
hold on;
plot(dates(186:258,:), n_cycle(186:258,:), 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(labor_rp, 1)
    fill([labor_rp(j, 1), labor_rp(j, 2), ...
          labor_rp(j, 2), labor_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 7);
hold on;
plot(labor_dates, demeaned_alt_labor, 'b', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(labor_rp, 1)
    fill([labor_rp(j, 1), labor_rp(j, 2), ...
          labor_rp(j, 2), labor_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Demeaned (Average Weekly Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

subplot(4, 2, 8);
hold on;
plot(labor_dates, l_cycle, 'k', 'LineWidth', 2);
ylim_vals = ylim;
for j = 1:size(labor_rp, 1)
    fill([labor_rp(j, 1), labor_rp(j, 2), ...
          labor_rp(j, 2), labor_rp(j, 1)], ...
         [ylim_vals(1), ylim_vals(1), ylim_vals(2), ylim_vals(2)], ...
         [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
title('Cycle (Average Weekly Labor Hours)');
xlabel('Date'); ylabel('Value');
grid on;
hold off;

