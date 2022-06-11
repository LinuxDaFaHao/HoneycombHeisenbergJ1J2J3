e0_3data = [];
for i = 0:55
    e0_3data = [e0_3data;load(['../data/energy',num2str(i),'beta','0.300000'])];
end

e0_4data = [];
for i = 0:55
    e0_4data = [e0_4data;load(['../data/energy',num2str(i),'beta','0.400000'])];
end

ecdata = [];
for i = 0:55
    ecdata = [ecdata;load(['../data/energy',num2str(i),'beta','0.440687'])];
end

e0_5data = [];
for i = 0:55
    e0_5data = [e0_5data;load(['../data/energy',num2str(i),'beta','0.500000'])];
end

e0_6data = [];
for i = 0:55
    e0_6data = [e0_6data;load(['../data/energy',num2str(i),'beta','0.600000'])];
end

energy_list = [mean(e0_3data),mean(e0_4data),mean(ecdata),mean(e0_5data),mean(e0_6data)];
C_list = [var(e0_3data),var(e0_4data),var(ecdata),var(e0_5data),var(e0_6data)];
beta_list = [0.3,0.4,0.440687, 0.5, 0.6];
%plot(beta_list, energy_list,'-o');
 plot(beta_list, C_list,'-o');