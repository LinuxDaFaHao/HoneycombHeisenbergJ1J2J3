L_set = [32, 48, 64,80];

beta_set = [0.2041  0.1538     0.1235 ...
    0.1031    0.0952      0.0885      0.0826  0.0810 ...
    0.0800 0.0790 0.0780 ...
    0.0775  ...
    0.077 0.0765 0.076 0.0755 ...
    0.0752  0.0745 0.074 0.0735  0.0730 0.0709    0.0690 ...
    0.0671      0.0637      0.0606      0.0578];


% beta_set = [  0.077 0.0765 0.076 0.0755 ...
%     0.0752  0.0745 0.074 0.0735  0.0730   0.0690 ...
%     0.0671      0.0637      0.0606      0.0578];
J1zz = 5.343;
J2zz = 0.2;
J3zz = -28.166;
Dzz = -0.113;
num_chain = 12;

prefix = '../data/';
eVtoK_const = 11.606;

extrapolation_selected_index = 2:4;

stiffness_set = zeros(numel(L_set), numel(beta_set));
energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
chiz_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
mxy2_set  = zeros(numel(L_set), numel(beta_set));
binder_ratio1_set = zeros(numel(L_set), numel(beta_set));
binder_ratio2_set = zeros(numel(L_set), numel(beta_set));
binder_ratio3_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);%15
        averaged_data=zeros(data_type_size, num_chain);
        
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
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
        %         stiffness = averaged_data(6,:);
        %             fprintf("stiffness = %.12f\n", mean(stiffness));
        %             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        msquare = sum(averaged_data([6,7,9,10,12,13],:))/3 ;
        %           msquare = sum(averaged_data(6:14,:))/3 ;
        fprintf("msquare = %.12f\n", mean(msquare));
        fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        %    fprintf("M_x2/M_y2 = %.12f\n", mean(averaged_data(8,:))/mean(averaged_data(9,:)));
        mxy2 = averaged_data(15,:);
        
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        mxy2_set(system_ind, beta_ind) = mean(mxy2);
        if strcmp(prefix, '../data2/')
            binder_ratio1 = averaged_data(16,:);
            binder_ratio2 = averaged_data(17,:);
            binder_ratio3 = averaged_data(18,:);
            binder_ratio1_set(system_ind, beta_ind) = mean(binder_ratio1);
            binder_ratio2_set(system_ind, beta_ind) = mean(binder_ratio2);
            binder_ratio3_set(system_ind, beta_ind) = mean(binder_ratio3);
        end
        
    end
end

T_set = eVtoK_const./beta_set;

%   plot(T_set, energy_set,'-^');hold on;
%   plot(T_set, specific_heat_set, '-o');hold on;
% plot(T_set, (binder_ratio1_set+binder_ratio2_set+binder_ratio3_set)/3,'-o');hold on;

   
%  plot(T_set, mxy2_set,'-o'); hold on;
% l=legend('$L=32$','$48$','$64$','$80$');%
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

plot(T_set, mxy2_set,'-x');hold on;
minf_set = size(1, numel(T_set));
for i = 1:size(mxy2_set, 2) % for every temperature
     y_data = mxy2_set(extrapolation_selected_index,i);
     x_data = 1./L_set(extrapolation_selected_index)';
     p = fit(x_data,y_data,'poly2');
     minf_set(i) = p.p3;
end
plot(T_set, minf_set,'-o');hold on;
l=legend('$L=32$','$48$','$64$','$80$','extrapolate to $L=\infty$');%
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
%  ylabel('specific heat $C$','Interpreter','latex');
% ylabel('$M^2(Q)$','Interpreter','latex');
 ylabel('complex order parameter','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Ylim',[0,inf]);


figure;
% minf_set = minf_set(T_set < 154);
minf_set = mxy2_set(end, T_set < 155 & T_set >110);
T_set = T_set(T_set<155&T_set >110);
loglog(156-T_set, minf_set,'-x');hold on;
 p = fit(log(156-T_set)',log(minf_set)','poly1');
 p


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$156 - T(K)$','Interpreter','latex');
%  ylabel('specific heat $C$','Interpreter','latex');
% ylabel('$M^2(Q)$','Interpreter','latex');
 ylabel('complex order parameter','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

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
