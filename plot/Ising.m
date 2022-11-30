L_set = [ 128,256];
%32, 64,

beta_set = [ 0.021 0.0207 0.0206 0.0205 0.02048 0.02045 0.02042 0.0204 0.02038 0.02036 0.02034 0.02032 0.0203 0.0202 0.0201 0.02  ];
% beta_set = [   0.021 0.0207 0.0206 0.0205  0.0204  0.0203 0.0202 0.0201 0.02  ];

J1 = 5.3;
J2 = 0.2;
J3 = -28;
D = 0.0;
num_chain = 12;

prefix = '../data/';
save_data_prefix = './plot_data/';
eVtoK_const = 1;

energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chi_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
% C2_set = zeros(numel(L_set), numel(beta_set));
% C4_set = zeros(numel(L_set), numel(beta_set));
mxy2_set  = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['ising-rank',num2str(0), 'Honeycomb', 'J1', num2str(J1,'%.6f'),...
            'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
            'D', num2str(D,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);
        averaged_data=zeros(data_type_size, num_chain);
        
        for i = 0:num_chain-1
            postfix = ['ising-rank',num2str(i), 'Honeycomb', 'J1', num2str(J1,'%.6f'),...
                'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
                'D', num2str(D,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            data = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = data(1:data_type_size);
        end
        
        chi = averaged_data(3,:);
        %             fprintf("chi_x = %.12f\n", mean(chix));
        %             fprintf("delta chi_x = %.12f\n", sqrt(var(chix)/num_chain));
        
        e = averaged_data(1,:);
        c =  averaged_data(2,:);
        %             fprintf("c = %.12f\n", mean(c));
        %             fprintf("delta c = %.12f\n", sqrt(var(c)/num_chain));
        %
        msquare = averaged_data(4,:);
%         m = max(msquare);
%         msquare = msquare(msquare > m/2)/3;
        fprintf("msquare = %.12f\n", mean(msquare));
        fprintf("n = %d\n", numel(msquare));
        fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        
        mxy2 = averaged_data(5,:);
        %         correlation2 = averaged_data(5,:);
        %         correlation4 = averaged_data(6,:);
        %         fprintf("C2_x = %.12f\n", mean(correlation2));
        %         fprintf("delta C2_x = %.12f\n", sqrt(var(correlation2)/num_chain));
        
        energy_set(system_ind, beta_ind) = mean(e);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chi_set(system_ind, beta_ind) = mean(chi);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        mxy2_set(system_ind, beta_ind) = mean(mxy2);
        %         C2_set(system_ind, beta_ind) = mean(correlation2);
        %         C4_set(system_ind, beta_ind) = mean(correlation4);
        
    end
end
T_set = eVtoK_const./beta_set;


% l=legend([h0,h1,h2],{'$\rho_s = 8T/\pi$','$\rho_s = 4T/\pi$','$\rho_s = 2T/\pi$'});
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',32);
% set(l,'Location','SouthWest');



% plot(T_set, specific_heat_set, '-o');hold on;

% plot(T_set, energy_set,'-o');hold on;
% plot(T_set, msquare_set,'-o');hold on;
% plot(T_set, mxy2_set,'-o'); hold on;
% plot(T_set, chi_set, '-o');hold on;
% plot(T_set, C4_set, '-x');hold on;
% plot(T_set, C2_set./C4_set, '-x');hold on;



eta_set = zeros(1, numel(T_set));
for i = 1:(numel(T_set))
    fit_x = mxy2_set(:,i)';
    p = fit(log(L_set)',log(fit_x'),'poly1');
    fprintf('T = %.5f, \t beta = %.5f,\t eta=%.5f\n',T_set(i), eVtoK_const/T_set(i), -p.p1);
    x = L_set(1):0.5:L_set(end);
    %     fl=loglog(x,exp(p.p2)*x.^p.p1,'-.'); hold on;
    eta_set(i) = -p.p1;
end
plot(T_set, eta_set,'-x');hold on;
plot(T_set, 0.25*ones(1, numel(T_set)),'-.');
plot(T_set, 1/9*ones(1, numel(T_set)),'-.');

% Tc = 48.98;
% a = 0.5;
% eta = 1/9; b = eta/2;
% for l = 1: numel(L_set)
%     L = L_set(l);
%     Tsub_set = T_set(T_set < Tc);
%     mxy_subset = sqrt(mxy2_set(l, T_set < Tc));
%     x_data = 1/L .* exp( a./ sqrt((Tc - Tsub_set)/Tc));
%     y_data = mxy_subset * L^b;
%     semilogx(x_data,y_data,'x');hold on;
% end
% set(gca, 'Xlim',[0,1e6]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('specific heat','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Xlim',[0,inf]);


