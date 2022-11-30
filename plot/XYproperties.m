L_set = [ 32, 64];

%  beta_set = [ 4 3 2 1.8 1.5 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5];
beta_set = [ 6 5 4.5 4 3.5 3 2.7 2.3 2 1.8 1.5 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5];
J1 = 1;
Hp = 0.01;
num_chain = 12;

prefix = '../data/';
eVtoK_const = 1;% 11.606;

stiffness_set = zeros(numel(L_set), numel(beta_set));
energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
binder_ratio_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        %summaryxy-rank0SquareJ11.000000Hp0.200000beta1.100000L128
        postfix = ['xy-rank',num2str(0), 'Square', 'J1', num2str(J1,'%.6f'),...
            'Hp', num2str(Hp,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);%6
        averaged_data=zeros(data_type_size, num_chain);
        
        for i = 0:num_chain-1
             postfix = ['xy-rank',num2str(i), 'Square', 'J1', num2str(J1,'%.6f'),...
            'Hp', num2str(Hp,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
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
        
        %
        %         stiffness = averaged_data(6,:);
        %             fprintf("stiffness = %.12f\n", mean(stiffness));
        %             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));

         msquare = averaged_data(5,:) ;
%         fprintf("msquare = %.12f\n", mean(msquare));
%         fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        %    fprintf("M_x2/M_y2 = %.12f\n", mean(averaged_data(8,:))/mean(averaged_data(9,:)));
       
        binder_ratio = averaged_data(6,:);

        
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
      
        msquare_set(system_ind, beta_ind) = mean(msquare);
      
        binder_ratio_set(system_ind, beta_ind) = mean(binder_ratio);

    end
end

T_set = eVtoK_const./beta_set;
%  plot(T_set, msquare_set,'-o');hold on;
% plot(T_set, binder_ratio_set ,'-o');hold on;
% plot(T_set, specific_heat_set ,'-o');hold on;
eta_set = zeros(1, numel(T_set));
% loglog(L_set, mxy2_set','x');hold on;
for i = 1:(numel(T_set))
    fit_x = msquare_set(:,i)';
    p = fit(log(L_set)',log(fit_x'),'poly1');
    fprintf('T = %.5f, \t beta = %.5f,\t eta=%.5f\n',T_set(i), eVtoK_const/T_set(i), -p.p1);
    x = L_set(1):0.5:L_set(end);
%     fl=loglog(x,exp(p.p2)*x.^p.p1,'-.'); hold on;
    eta_set(i) = -p.p1;
end
plot(T_set, eta_set,'-x');hold on;
plot(T_set, 0.25*ones(1, numel(T_set)),'-.');
plot(T_set, 1/9*ones(1, numel(T_set)),'-.');

%   plot(T_set, energy_set,'-^');hold on;
%   plot(T_set, specific_heat_set, '-o');hold on;
%  plot(T_set, binder_ratio1_set,'-o');hold on;
%   plot(T_set, binder_ratio2_set,'-o');hold on;
%    plot(T_set, binder_ratio3_set,'-o');hold on;
%    plot(T_set, msquare_set,'-o');hold on;
%     plot(T_set, mxy2_set,'-o'); hold on;

% plot(T_set, chix_set, '-o');hold on;
% plot(T_set, chiy_set, '-o');hold on;
% plot(T_set, chiz_set, '-o');hold on;


% Tc_set = [158.5,157.5,157];%K


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');
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
