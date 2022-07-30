L_set = [64, 72, 100,128, 144];
beta_set = [ 0.2041    0.1887    0.1754   0.1639    0.1538  0.1449    0.1370    0.1299    0.1235    0.1176    0.1124    0.1075    0.1031    0.0990  0.0952 ...
    0.0917    0.0885    0.0855    0.0826    0.0800    0.0775    0.0752    0.0730  0.0709    0.0690    0.0671    0.0654    0.0637    0.0621    0.0606    0.0592    0.0578   0.0565    0.0552];
J1zz = 5.3;
J2zz = 0.2;
J3zz = -38;
Dzz = -0.113;
num_chain = 5;

prefix = '../data2/';
save_data_prefix = './plot_data2/';
eVtoK_const = 11.606;

% stiffness_set = zeros(numel(L_set), numel(beta_set));
% for system_ind = 1:numel(L_set)
%     L = L_set(system_ind);
%     N = 2 * L^2;
%     
%     T_set = 1./beta_set * eVtoK_const;
%     chi_x_set = zeros(1, numel(beta_set));
%     chi_y_set = zeros(1, numel(beta_set));
%     chi_z_set = zeros(1, numel(beta_set));
%     c_set = zeros(1, numel(beta_set));
%     energy_set = zeros(1, numel(beta_set));
%     
%     for beta_ind = 1:numel(beta_set)
%         beta = beta_set(beta_ind);
%         
%         fprintf('beta=%.6f\n', beta);
%         
%         postfix = ['hei-rank',num2str(0), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
%             'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
%             'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
%         data_example = load([prefix, 'sum_spin0', postfix]);
%         chain_length = numel(data_example);
%         sx_data=zeros(chain_length, num_chain);
%         sy_data=zeros(chain_length, num_chain);
%         sz_data=zeros(chain_length, num_chain);
%         energy_data=zeros(chain_length, num_chain);
%         stiffness_data = zeros(chain_length, num_chain);
%         
%         direct_data_file = [save_data_prefix,'/hei','Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
%             'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
%             'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L),'.mat'];
%         if exist(direct_data_file,'file')
%             load(direct_data_file);
%         else
%             
%             for i = 0:num_chain-1
%                 postfix = ['hei-rank',num2str(i), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
%                     'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
%                     'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
%                 sx_data(:, i+1) = load([prefix, 'sum_spin0', postfix]);
%                 sy_data(:, i+1) = load([prefix, 'sum_spin1', postfix]);
%                 sz_data(:, i+1) = load([prefix, 'sum_spin2', postfix]);
%                 energy_data(:, i+1) = load([prefix, 'energy', postfix]);
%                 stiffness_data(:,i+1) = load([prefix, 'stiffness', postfix]);
%             end
%             
%             
%             [chix, delta_chix] = CalMagSusp(sx_data(end/2:end,:), N);
%             
%             fprintf("chi_x = %.12f\n", mean(chix));
%             fprintf("delta chi_x = %.12f\n", sqrt(var(chix)/num_chain));
%             
%             
%             [chiy, delta_chiy] = CalMagSusp(sy_data(end/2:end,:), N);
%             fprintf("chi_y = %.12f\n", mean(chiy));
%             fprintf("delta chi_y = %.12f\n", sqrt(var(chiy)/num_chain));
%             
%             
%             [chiz, delta_chiz] = CalMagSusp(sz_data(end/2:end,:), N);
%             fprintf("chi_z = %.12f\n", mean(chiz));
%             fprintf("delta chi_z = %.12f\n", sqrt(var(chiz)/num_chain));
%             
%             chi=chix+chiy+chiz;
%             fprintf("chi(x+y+z) = %.12f\n", mean(chi));
%             
%             energy_data = energy_data(end/2+1:end,:);
%             energy = mean(energy_data)/N;
%             fprintf("e = %.12f\n", mean(energy));
%             fprintf("delta e = %.12f\n", sqrt(var(energy)/num_chain));
%             
%             c = (mean(energy_data.^2) - mean(energy_data).^2) * beta^2 /N;
%             fprintf("delta c = %.12f\n", sqrt(var(c)/num_chain));
%             
%             stiffness = mean(stiffness_data(end/2+1:end,:));
%             fprintf("stiffness = %.12f\n", mean(stiffness));
%             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
%             
%             save(direct_data_file, 'chix','chiy','chiz','energy','c','stiffness');
%         end
%         
%         stiffness_set(system_ind, beta_ind) = mean(stiffness);
%     end 
% end
load('Stiffness.mat','stiffness_set');
h = plot(T_set, stiffness_set, '-o');hold on;
h0 = plot(T_set, 8./beta_set/pi,'-.'); hold on;
l=legend(h0,'$\rho_s = 8T/\pi$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',32);
set(l,'Location','SouthWest');


Tc_set = [152.059154,151.6366,150.251,149.25648,148.6865];%K


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('spin stiffness $\rho_s$(meV)','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
set(gca, 'Ylim',[0,inf]);


axes('Position',[0.6,0.7,0.25,0.2]);
fit_x = 1./(log(L_set).^2);
plot(fit_x,Tc_set,'o'); hold on;

p = fit(fit_x',Tc_set','poly1');
T_BKT = p.p2;
fprintf('T_BKT=%.5f\n',T_BKT);
plot_x = 0:max(fit_x)/100:max(fit_x);
plot(plot_x, p.p1*plot_x + p.p2,'-.');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5); 
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% ylabel('$T_{\mathrm{BKT}}/K$','Interpreter','latex');
ylabel('$T^*(K)$','Interpreter','latex');
xlabel('$(\ln L)^2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

set(gca, 'Xlim',[0,inf]);



function [chi, delta_chi] = CalMagSusp(smu_data, N)
% smu_data is matrix
smu_square = smu_data .* smu_data;
chi_list = (mean(smu_square) - mean(smu_data).^2)/N;
chi = mean(chi_list);
delta_chi = 0;
% delta_chi = sqrt( var(chi_list)/numel(chi_list));
end
