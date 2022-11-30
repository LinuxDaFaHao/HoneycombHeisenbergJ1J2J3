L_set = [32,48,64,80];

% beta_set = [ 0.2041  0.1538     0.1235 ...
%     0.1031    0.0952      0.0885      0.0826  0.0810 ...
%     0.0800 0.0790 0.0780 ...
%     0.0775  ...
%     0.077 0.0765 0.076 0.0755 ...
%     0.0752  0.0745 0.074 0.0735  0.0730  0.0709   0.0690 ...
%     0.0671           0.0606      0.0578];
beta_set = [ 0.2041      0.1235 ...
    0.1031    0.0952      0.0885      0.0826  0.0810 ...
    0.0800 0.0790 0.0780 ...
    0.0775  ...
    0.077 0.0765 0.076 0.0755 ...
    0.0752  0.0745 0.074 0.0735  0.0730  0.0709   0.0690 ...
    0.0671           0.0606      0.0578];



J1zz = 5.31;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.113;
num_chain = 12;

prefix = '../../data/';
eVtoK_const = 11.606;


correlation2_set = zeros(numel(L_set), numel(beta_set));
correlation4_set = zeros(numel(L_set), numel(beta_set));
ratio_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'correlation_summary', postfix]);
        data_type_size = numel(test_load);%
        averaged_data=zeros(data_type_size, num_chain);
        
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            data = load([prefix, 'correlation_summary', postfix]);
            averaged_data(:,i + 1) = data(1:data_type_size);
        end
      
      
        correlation2 = averaged_data(1,:);
        correlation4 = averaged_data(2,:);
        ratio = averaged_data(3,:);
        
        correlation2_set(system_ind, beta_ind) = mean(correlation2);
        correlation4_set(system_ind, beta_ind) = mean(correlation4);
       
        ratio_set(system_ind, beta_ind) = mean(ratio);
    end
end

T_set = eVtoK_const./beta_set;

% plot(T_set, energy_set,'-^');hold on;
%     plot(T_set, specific_heat_set, '-o');hold on;
%  plot(T_set(10:end), (binder_ratio1_set(:,10:end) + binder_ratio2_set(:,10:end) + binder_ratio3_set(:,10:end))/3,'-o');hold on;
% start_num = 6;
%  plot(T_set(start_num:end), (binder_ratio1_set(:,start_num:end) + binder_ratio2_set(:,start_num:end) + binder_ratio3_set(:,start_num:end))/3,'-o');hold on;
% plot(T_set(start_num:end), 5/3*ones(size(T_set(start_num:end))),'-.');hold on;
plot(T_set,ratio_set,'-o');hold on;
% plot(T_set, binder_ratio2_set,'-o');hold on;
% plot(T_set, binder_ratio3_set,'-o');hold on;
%  plot(T_set, msquare_set,'-o');hold on;
%   plot(T_set, mxy2_set,'-o'); hold on;

% plot(T_set, chix_set, '-o');hold on;
% plot(T_set, chiy_set, '-o');hold on;
% plot(T_set, chiz_set, '-o');hold on;


% Tc_set = [158.5,157.5,157];%K


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
%  ylabel('specific heat','Interpreter','latex');
%  ylabel('$M^2(Q)$','Interpreter','latex');
%  ylabel('complex order parameter','Interpreter','latex');
% ylabel('binder ratio $\langle M^4\rangle / \langle M^2\rangle^2$','Interpreter','latex');
ylabel('$g(L)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

l=legend('$L=32$','$48$','$80$');%
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');
%  set(gca, 'Ylim',[0,inf]);
%
%
% axes('Position',[0.6,0.7,0.25,0.2]);
% fit_x = 1./(log(L_set).^2);
% plot(fit_x,Tc_set,'o'); hold on;
%
% p = fit(fit_x',Tc_set','poly1');
% T_BKT = p.p2;
% fprintf('T_BKT=%.5f\n',T_BKT);
% plot_x = 0:max(fit_x)/100:max(fit_x);
% plot(plot_x, p.p1*plot_x + p.p2,'-.');
%
%
%
% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% % ylabel('$T_{\mathrm{BKT}}/K$','Interpreter','latex');
% ylabel('$T^*(K)$','Interpreter','latex');
% xlabel('$(\ln L)^2$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);
%
% set(gca, 'Xlim',[0,inf]);
%
