%{
Oct 12, ../../data
1: energy; 2: specific heat; 3~5: anti-ferromagnetic
susceptibilit                    y(S[Q1,Q2,Q3]);
6~8: AF magnetization (Q1, Q2, Q3);
9~11: binder ratio of magnetization
12: binder ratio of magnetization(without specify momentum)
13: C3 symmetry order parameter
14: C3 symmetry order parameter susceptibilities
15: complex order parameter
16: binder ratio of complex order parameter
%}
L_set = [64,80];
%

% beta_set = [0.22 0.218 0.215 0.213 0.21 0.208 0.205 0.203 0.2 0.18 0.16 0.14 0.12 0.1 0.09 0.08 0.07];
beta_set = [0.16  0.14 0.13 0.128 0.125 0.123 0.12 0.118 0.115 0.113 0.11 0.1 0.09 ];
beta_set = [beta_set,0.119  0.117 0.116 0.114  0.1138    0.1137    0.1136  0.1135 0.1134    0.1133    0.1132    0.1131 0.1125 0.112 0.1115 0.1105];
beta_set = sort(beta_set);
%
J1 = 1;
J2 = 0;
J3 = -5;
D = 0.0;
num_chain = 12;

prefix = '../../data/';
eVtoK_const = 1;

specific_heat_set = zeros(numel(L_set), numel(beta_set));
specific_heat_error_set = zeros(numel(L_set), numel(beta_set));
chi_set = zeros(numel(L_set), numel(beta_set));
energy_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
mxy2_set  = zeros(numel(L_set), numel(beta_set));
binder_ratio1_set = zeros(numel(L_set), numel(beta_set));
binder_ratio2_set = zeros(numel(L_set), numel(beta_set));
binder_ratio3_set = zeros(numel(L_set), numel(beta_set));
binder_ratio_set = zeros(numel(L_set), numel(beta_set));
binder_ratiom_set = zeros(numel(L_set), numel(beta_set));
bc3_set  = zeros(numel(L_set), numel(beta_set));
chic3_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    
    for beta_ind = 1:numel(beta_set)
        
        
        beta = beta_set(beta_ind);
        geometry = 'Honeycomb';
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['ising-rank',num2str(0), geometry, 'J1', num2str(J1,'%.6f'),...
            'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
            'D', num2str(D,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        %       test_load = load([prefix, 'summary', postfix]);
        %       data_type_size = numel(test_load);%15
        data_type_size = 16;
        averaged_data=zeros(data_type_size, num_chain);
        
        for i = 0:num_chain-1
            postfix = ['ising-rank',num2str(i), geometry, 'J1', num2str(J1,'%.6f'),...
                'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
                'D', num2str(D,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            file_name = [prefix, 'summary', postfix];
            if(exist(file_name,'file'))
                data = load(file_name);
                averaged_data(:,i + 1) = data(1:data_type_size);
            else
                averaged_data(:,i + 1) = NaN;
            end
        end
        
        energy = averaged_data(1,:);
        
        c =  averaged_data(2,:);
        c = c(~isnan(c));
        
        
        chi = sum(averaged_data([3,4,5]))/3 /N *beta;
        chi = chi(~isnan(chi));
        
        
        %         msquare = sum(averaged_data([6,7,9,10,12,13],:))/3 ;
        %         msquare = sum(averaged_data(6:14,:))/3 ;
        msquare = sum(averaged_data([6,7,8],:))/3 ;
        msquare = msquare(~isnan(msquare));
        %         fprintf("msquare = %.12f\n", mean(msquare));
        %         fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        
        binder_ratio1 = averaged_data(9,:);
        binder_ratio2 = averaged_data(10,:);
        binder_ratio3 = averaged_data(11,:);
        binder_ratio = averaged_data(12,:);
        bc3 = averaged_data(13,:);
        chic3 = averaged_data(14,:);
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        specific_heat_error_set(system_ind, beta_ind) = std(c);
        chi_set(system_ind, beta_ind) = mean(chi);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        
        
        binder_ratio1_set(system_ind, beta_ind) = mean(binder_ratio1);
        binder_ratio2_set(system_ind, beta_ind) = mean(binder_ratio2);
        binder_ratio3_set(system_ind, beta_ind) = mean(binder_ratio3);
        binder_ratio_set(system_ind, beta_ind) = mean(binder_ratio);
        bc3_set(system_ind, beta_ind) = mean(bc3);
        chic3_set(system_ind, beta_ind) = mean(chic3);
        mxy2 = averaged_data(15,:);
        mxy2_set(system_ind, beta_ind) = mean(mxy2);
        binder_ratiom = averaged_data(16,:);
        binder_ratiom_set(system_ind, beta_ind) = mean(binder_ratiom);
        
    end
end
T_set = eVtoK_const./beta_set;


% l=legend([h0,h1,h2],{'$\rho_s = 8T/\pi$','$\rho_s = 4T/\pi$','$\rho_s = 2T/\pi$'});
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',32);
% set(l,'Location','SouthWest');


%  beta2 = beta_set(1:end-1) + diff(beta_set)/2;
%     C = -beta2.^2 .* (diff(energy_set)./diff(beta_set));
%     plot(1./beta2, C,'-x');hold on;

% plot(T_set, specific_heat_set, '-o');hold on;
% errorbar(T_set, specific_heat_set, specific_heat_error_set, '-o');hold on;

%      plot(T_set, energy_set,'-.');hold on;
%   plot(T_set, msquare_set ,'-.');hold on;
%  plot(T_set, mxy2_set,'-o'); hold on;
%     plot(T_set, chi_set, '-o');hold on;


% start_num = 1;
% plot(T_set(start_num:end), (binder_ratio1_set(:,start_num:end) + binder_ratio2_set(:,start_num:end) + binder_ratio3_set(:,start_num:end))/3,'-o');hold on;
% plot(T_set, binder_ratio3_set,'-o');

%  plot(T_set, binder_ratio_set,'-o');hold on;

% plot(T_set, binder_ratiom_set,'-o');hold on;

%  plot(T_set, bc3_set, '-x');hold on;
% plot(T_set, chic3_set, '-o');hold on;


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


