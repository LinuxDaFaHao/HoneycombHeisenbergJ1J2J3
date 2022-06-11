L = 64;
sx_data=[];
for i = 0:4
     sx_data = [sx_data,load(['../data/sum_spin0',num2str(i),'beta','0.440687'])];
%     sx_data = [sx_data,load(['../data/cluster_size',num2str(i),'beta','0.440687'])];
     data = sx_data(:,i+1);
     data = reshape(data,10000,[]);
     data = data.*data/L^2/L^2;
     plot(mean(data),'-o');hold on;
end



sx_data = sx_data(end/2+1:end,:);
sx2_data = sx_data.*sx_data;
chi = mean(sx2_data)/L^2;
normchi = chi/L^2;
fprintf("chi/L^2 = %.12f\n", mean(normchi));
fprintf("var chi/L^2 = %.12f\n", var(normchi));


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('step','Interpreter','latex');
ylabel('$\chi$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 