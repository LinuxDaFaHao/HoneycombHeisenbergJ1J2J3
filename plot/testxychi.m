L = 64;
beta = 1.12;
sx_data=[];
sy_data=[];
num_chain = 5;
for i = 0:num_chain-1
    %      raw_data = load(['../data/sum_spin0',num2str(i),'beta',num2str(beta,'%.6f')]);
    raw_data = load(['../data/sum_spin0',num2str(i)]);
    plot(mean(reshape(raw_data.*raw_data,100000,[]))/L^4,'-o');hold on;
    sx_data = [sx_data,raw_data];
%     sy_data = [sy_data,load(['../data/sum_spin1',num2str(i),'beta',num2str(beta,'%.6f')])];
    sy_data = [sy_data,load(['../data/sum_spin1',num2str(i)])];
end
sx_data = sx_data(100000:end,:);
sx2_data = sx_data.*sx_data;
chix = mean(sx2_data)/L^2;
normchix = chix/L^2;
fprintf("chi_x/L^2 = %.12f\n", mean(normchix));
fprintf("var chi_x/L^2 = %.12f\n", var(normchix));


sy_data = sy_data(end/2:end,:);
sy2_data = sy_data.*sy_data;
chiy = mean(sy2_data)/L^2;
normchiy = chiy/L^2;
fprintf("chi_y/L^2 = %.12f\n", mean(normchiy));
fprintf("var chi_y/L^2 = %.12f\n", var(normchiy));


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('step','Interpreter','latex');
ylabel('$\chi$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);