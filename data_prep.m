%% Preparing Data for Bayesian Estimation
JPT = readtable("JPT.xlsx");
y = JPT.y;
c = JPT.c;
i = JPT.i;
ppi = JPT.ppi; 
R = JPT.R;  
L = JPT.L;  
w = JPT.w;

save('model_data.mat', 'y', 'c', 'i', 'L', 'ppi', 'w', 'R');



