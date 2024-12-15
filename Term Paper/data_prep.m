%% Preparing Data for Bayesian Estimation
data = readtable("data.xlsx");
vix = data.vix;
earnings = data.earnings;
savings = data.savings;
credit = data.credit;
housing = data.housing;
durables = data.durables;


save('model_data.mat', 'vix', 'earnings', 'savings', 'credit', 'housing', 'durables');