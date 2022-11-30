L_set = [8,8];

% beta_set = [  0.2041   0.1538 0.1235     ...
%     0.1031    0.0952      0.0885      0.0826  0.0810 ...
%     0.0800 0.0790 0.0780 ...
%     0.0775  ...
%     0.077 0.0765 0.076 0.0755 ...
%     0.0752  0.0745  0.074   0.0735  0.0730     0.0690 ...
%     0.0671  0.0637  0.0606  0.0578];
% 0.0709
beta_set = [  0.2041  ];
J1zz = 5.29;
J2zz = 0.2;
J3zz = -28;
Dzz = 0;
num_chain = 12;

prefix = '../data/';
eVtoK_const = 11.606;

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
        
        %         fprintf('beta=%.6f\n', beta);
        if(system_ind==1)
            postfix = ['hei-rank',num2str(0), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        else
            postfix = ['hei-rank',num2str(0), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        end
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);%15
        averaged_data=zeros(data_type_size, num_chain);
        
        
        for i = 0:num_chain-1
            if(system_ind==1)
                postfix = ['hei-rank',num2str(i), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
                    'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                    'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            else
                postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
                    'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                    'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            end
            data = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = data(1:data_type_size);
        end
        disp(mean(averaged_data'));
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
        %          msquare = sum(averaged_data([6,7,9,10,12,13]+1,:))/3 ;
        
        msquare = sum(averaged_data(7:15,:))/3 ;
        %          msquare = sum(averaged_data(10:12,:)) ;
        fprintf("msquare = %.12f\n", mean(msquare));
%         fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        %    fprintf("M_x2/M_y2 = %.12f\n", mean(averaged_data(8,:))/mean(averaged_data(9,:)));
        %         mxy2 = averaged_data(15,:);
        
        %         stiffness_set(system_ind, beta_ind) = mean(stiffness);
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        %         mxy2_set(system_ind, beta_ind) = mean(mxy2);
        binder_ratio1 = averaged_data(16,:);
        binder_ratio2 = averaged_data(17,:);
        binder_ratio3 = averaged_data(18,:);
        binder_ratio1_set(system_ind, beta_ind) = mean(binder_ratio1);
        binder_ratio2_set(system_ind, beta_ind) = mean(binder_ratio2);
        binder_ratio3_set(system_ind, beta_ind) = mean(binder_ratio3);
        
    end
end

T_set = eVtoK_const./beta_set;

% plot(T_set, energy_set,'-^');hold on;
%      plot(T_set, specific_heat_set, '-o');hold on;
% start_num = 1;
%  plot(T_set(start_num:end), (binder_ratio1_set(:,start_num:end) + binder_ratio2_set(:,start_num:end) + binder_ratio3_set(:,start_num:end))/3,'-o');hold on;

%  plot(T_set, binder_ratio1_set,'-o');hold on;
% plot(T_set, binder_ratio2_set,'-o');hold on;
% plot(T_set, binder_ratio3_set,'-o');hold on;
%   plot(T_set(start_num:end), msquare_set(:,start_num:end),'-o');hold on;
%   plot(T_set, mxy2_set,'-o'); hold on;

% plot(T_set, chix_set, '-o');hold on;
% plot(T_set, chiy_set, '-o');hold on;
% plot(T_set, chiz_set, '-o');hold on;
% eta_set = zeros(1, numel(T_set));
% % loglog(L_set, mxy2_set','x');hold on;
% for i = 1:(numel(T_set))
%     fit_x = msquare_set(:,i)';
%     p = fit(log(L_set)',log(fit_x'),'poly1');
%     fprintf('T = %.5f, \t beta = %.5f,\t eta=%.5f\n',T_set(i), eVtoK_const/T_set(i), -p.p1);
%     x = L_set(1):0.5:L_set(end);
% %     fl=loglog(x,exp(p.p2)*x.^p.p1,'-.'); hold on;
%     eta_set(i) = -p.p1;
% end
% plot(T_set, eta_set,'-x');hold on;
% plot(T_set, 0.25*ones(1, numel(T_set)),'-.');

% Tc_set = [158.5,157.5,157];%K


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
%  ylabel('specific heat','Interpreter','latex');
%   ylabel('$M^2(Q)$','Interpreter','latex');
% ylabel('complex order parameter','Interpreter','latex');
ylabel('binder ratio $\langle M^4\rangle / \langle M^2\rangle^2$','Interpreter','latex');
%  ylabel('binder ratio $\langle m^4\rangle / \langle m^2\rangle^2$','Interpreter','latex');
%  ylabel('complex order parameter','Interpreter','latex');
% ylabel('$\eta$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Ylim',[0,inf]);
