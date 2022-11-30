L_set = [64,128, 192];

beta_set = [0.1031  0.0855  0.0800 0.0752 0.075 0.0745 0.074 0.0735 ...
    0.0734 0.07335 0.0733 0.07325 0.0732 0.0731 ...
    0.073 0.0725 0.072 0.0715 0.071 0.0709 0.0671 ];
%
J1zz = 5.3;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.3;
num_chain = 12;

prefix = '../data/';
save_data_prefix = './plot_data/';
eVtoK_const = 11.606;

stiffness_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
chiz_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
mxy2_set  = zeros(numel(L_set), numel(beta_set));

for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    chi_x_set = zeros(1, numel(beta_set));
    chi_y_set = zeros(1, numel(beta_set));
    chi_z_set = zeros(1, numel(beta_set));
    c_set = zeros(1, numel(beta_set));
    energy_set = zeros(1, numel(beta_set));
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f, T = %.6f\n', beta, eVtoK_const/beta);
        
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
            data = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = data(1:data_type_size);
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

        %
        stiffness = averaged_data(6,:);
        %             fprintf("stiffness = %.12f\n", mean(stiffness));
        %             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        msquare = averaged_data(8,:);
%         m = max(msquare);
%         msquare = msquare(msquare > m/2)/3;
        fprintf("msquare = %.12f\n", mean(msquare));
%         fprintf("n = %d\n", numel(msquare));
        fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        
        mxy2 = averaged_data(11,:);
        
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        mxy2_set(system_ind, beta_ind) = mean(mxy2);
    end
end
T_set = eVtoK_const./beta_set;


% h = plot(T_set, stiffness_set, '-o');hold on;
% h0 = plot(T_set, 8./beta_set/pi,'-.'); hold on;
% h1 = plot(T_set, 4./beta_set/pi,'-.'); hold on;
% h2 = plot(T_set, 2./beta_set/pi,'-.'); hold on;
% l=legend([h0,h1,h2],{'$\rho_s = 8T/\pi$','$\rho_s = 4T/\pi$','$\rho_s = 2T/\pi$'});
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',32);
% set(l,'Location','SouthWest');



%plot(T_set, specific_heat_set, '-o');hold on;

plot(T_set, msquare_set,'-o');hold on;
plot(T_set, mxy2_set,'-o'); hold on;
% plot(T_set, chix_set, '-o');hold on;
% plot(T_set, chiy_set, '-o');hold on;
% plot(T_set, chiz_set, '-o');hold on;


% Tc_set = [158.5,157.5,157];%K


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('specific heat','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
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
