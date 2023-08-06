L_set = [32, 48, 64, 128];
% L_set = [32, 48, 64, 72,100];
%
% beta_set = [ 0.2857 0.2326   0.1887      0.1639       0.1370        0.1176      0.1031      0.0952    0.0917    0.0885    0.0855    0.0826    0.0800   0.0775   0.0752   0.0730  0.0709    0.0690    0.0671    0.0654    0.0637    0.0621    0.0606    0.0592    0.0578   0.0565    0.0552];
% beta_set = [  beta_set,  0.0770 0.0760  0.0745 0.0740  ]; % L = ,
% . don't delete
% beta_set = [0.2326  0.1887   0.1639   0.1449     0.1299    0.1176    0.1075    0.0990  0.0917   0.0855   0.0800    0.0775  0.0770    0.0760    0.0752    0.0745    0.0740    0.0730    0.0709    0.0690    0.0671    0.0654    0.0637];
beta_set = [   0.1639    0.1299    0.1176    0.1075    0.0990  0.0917   0.0855   0.0800    0.0775  0.0770    0.0760    0.0752    0.0745    0.0740    0.0730    0.0709    0.0690    0.0671    0.0654    0.0637];

beta_set = sort(beta_set);
J1zz = 5.3;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.113;
num_chain = 5;

prefix = '../../data/';
save_data_prefix = '../plot_data/';
eVtoK_const = 11.606;

stiffness_set = zeros(numel(L_set), numel(beta_set));
energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
chiz_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));


for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);
        averaged_data=zeros(data_type_size, num_chain);
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            temp = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = temp(1:data_type_size);
        end
        energy = averaged_data(1,:);
        %             fprintf("energy = %.12f\n", mean(energy));
        %             fprintf("delta e = %.12f\n", sqrt(var(energy)/num_chain));
        %
        c =  averaged_data(2,:);
        %             fprintf("C = %.12f\n", mean(c));
        %             fprintf("delta C = %.12f\n", sqrt(var(c)/num_chain));
        chix = averaged_data(3,:);
        chiy = averaged_data(4,:);
        chiz = averaged_data(5,:);
        stiffness = averaged_data(6,:);
        %stiffness = stiffness(stiffness>0);
        fprintf("stiffness = %.12f\n", mean(stiffness));
        fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        
        msquare = sum(averaged_data([7,8,10,11,13,14],:))/3 ; %xy-component
        % msquare = sum(averaged_data([9,12,15],:))/3 ; %z-component
        % msquare = sum(averaged_data(7:15,:))/3 ;
        energy_set(system_ind, beta_ind) = mean(energy);
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
    end
end
T_set = eVtoK_const./beta_set;
%    plot(T_set, specific_heat_set, '-o');hold on;
% plot(T_set, stiffness_set,'-o');hold on;

% plot(T_set, energy_set,'-o');hold on;
% plot(T_set, chix_set+ chiy_set + chiz_set,'-o'); hold on;
% h = plot(T_set,   msquare_set, '-o');hold on;
% h2 = plot(T_set, 2./beta_set/pi,'-.'); hold on;
% % l=legend([h0,h1,h2],{'$\rho_s = 8T/\pi$','$\rho_s = 4T/\pi$','$\rho_s = 2T/\pi$'});
% l=legend([h2],{'$\rho_s = \frac{2T}{\pi}$'});
% l=legend('32','48','64','128');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

% T=text(50,10,['$L=32, 48, 64, 72, 100$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);


plot(T_set, 2*sqrt(msquare_set),'-x');hold on;
% minf_set = size(1, numel(T_set));
% for i = 1:size(msquare_set, 2) % for every temperature
%      y_data = sqrt(msquare_set(:,i));
%      x_data = 1./L_set';
%      p = fit(x_data,y_data,'poly2');
%      minf_set(i) = max(p.p3,0.0);
% end
% plot(T_set, minf_set,'-o');hold on;

l=legend('$L=32$','$48$','$64$','$128$');%,'extrapolate to $L=\infty$'
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');
% 
% Tc = 152;
% T_subset = T_set(T_set<Tc);
% minf_subset = msquare_set(end,T_set<Tc);
% loglog( (Tc-T_subset), (minf_subset),'-o');hold on;
% p = fit(log(Tc-T_subset(1:4))',log(minf_subset(1:4))','poly1');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
%  ylabel('spin stiffness $\rho_s$(meV)','Interpreter','latex');
% ylabel('AF magnetization $M_{x}^2(Q)+M_{y}^2(Q)$','Interpreter','latex');
% ylabel('specific heat $C (k_B)$','Interpreter','latex');
ylabel('inplane AF magnetization $M(\mu_B)$','Interpreter','latex');

set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
% set(gca, 'Ylim',[0,inf]);
% set(gca, 'Xlim',[100,inf]);


% Tc_set = [161.98, 157.42, 155.24, 154.46,153.06];%K
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
% ylabel('$T_{\mathrm{BKT}}/K$','Interpreter','latex');
% ylabel('$T^*(K)$','Interpreter','latex');
% xlabel('$(\ln L)^2$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);
%
% set(gca, 'Xlim',[0,inf]);
%
