L_set = [144];
for L = L_set
% L = 64;
N = 2 * L^2;
% beta_set = [ 0.05,0.06,0.07,0.077,0.08,0.09,0.10,0.12,0.14,0.16,0.18,0.3];
% beta_set = [3.3333    2.0000    1.4286    1.1111    0.9091    0.7692 ...
% 0.6667    0.5882    0.5263    0.4762  0.4348    0.4000    0.3704 ...
% 0.3448    0.3226    0.3030    0.2857    0.2703    0.2564    0.2439  ...
% 0.2326   0.2222    0.2128    0.2041    0.1961    0.1887    0.1818   ...
% 0.1754    0.1695    0.1639    0.1587    0.1538  0.1493    0.1449    ...
% 0.1408    0.1370    0.1333    0.1299    0.1266    0.1235    0.1205    ...
% 0.1176    0.1149   0.1124    0.1099    0.1075    0.1053    0.1031    ...
% 0.1010    0.0990    0.0971    0.0952    0.0935    0.0917   0.0901    ...
% 0.0885    0.0870    0.0855    0.0840    0.0826    0.0813    0.0800    ...
% 0.0787    0.0775    0.0763  0.0752    0.0741    0.0730    0.0719    ...
% 0.0709    0.0699    0.0690    0.0680    0.0671    0.0662    0.0654 ...
% 0.0645    0.0637    0.0629    0.0621    0.0613    0.0606    0.0599    ...
% 0.0592    0.0585    0.0578    0.0571  0.0565    0.0559    0.0552    ...
% 0.0546    0.0541    0.0535    0.0529    0.0524    0.0518    0.0513    0.0508  0.0503];
%   beta_set = beta_set(beta_set < 0.1754  );
beta_set = [ 0.2041    0.1887    0.1754   0.1639    0.1538  0.1449    0.1370    0.1299    0.1235    0.1176    0.1124    0.1075    0.1031    0.0990  0.0952 ...
  0.0917    0.0885    0.0855    0.0826    0.0800    0.0775    0.0752    0.0730  0.0709    0.0690    0.0671    0.0654    0.0637    0.0621    0.0606    0.0592    0.0578   0.0565    0.0552];
% beta_set = beta_set(beta_set < 0.22 );
% J1zz = 5.3;
% J2zz = 0.2;
% J3zz = -38;
% Dzz = -0.01;
J1zz = -38;
J2zz = 0;
J3zz = 0;
Dzz = -0.113;
num_chain = 5;

% beta_set = [ 0.05,0.06,0.07,0.08,0.09,0.10,0.15,0.2,0.3,0.4,0.5,1,2, 10,100];
% J1zz = 5.343;
% J2zz = 0.190;
% J3zz = -28.166;
% Dzz = -0.113;
% num_chain = 5;
prefix = '../data2/';
save_data_prefix = './plot_data2/';
eVtoK_const = 11.606;
T_set = 1./beta_set * eVtoK_const;
chi_x_set = zeros(1, numel(beta_set));
chi_y_set = zeros(1, numel(beta_set));
chi_z_set = zeros(1, numel(beta_set));
c_set = zeros(1, numel(beta_set));
energy_set = zeros(1, numel(beta_set));
stiffness_set = zeros(1, numel(beta_set));
for j = 1:numel(beta_set)
    beta = beta_set(j);
    
    fprintf('beta=%.6f\n', beta);
    
    postfix = ['hei-rank',num2str(0), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
        'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
        'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
    data_example = load([prefix, 'sum_spin0', postfix]);
    chain_length = numel(data_example);
    sx_data=zeros(chain_length, num_chain);
    sy_data=zeros(chain_length, num_chain);
    sz_data=zeros(chain_length, num_chain);
    energy_data=zeros(chain_length, num_chain);
    stiffness_data = zeros(chain_length, num_chain);
    
    direct_data_file = [save_data_prefix,'/hei','Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
        'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
        'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L),'.mat'];
    if exist(direct_data_file,'file')
        load(direct_data_file);
    else
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            sx_data(:,i+1) = load([prefix, 'sum_spin0', postfix]);
            
            %     data = sx_data(:,i+1);
            %     data = reshape(data,1000,[]);
            %     data = data.*data/L^2;
            %     plot(mean(data),'-o');hold on;
            
            sy_data(:, i+1) = load([prefix, 'sum_spin1', postfix]);
            sz_data(:, i+1) = load([prefix, 'sum_spin2', postfix]);
            energy_data(:, i+1) = load([prefix, 'energy', postfix]);
            stiffness_data(:,i+1) = load([prefix, 'stiffness', postfix]);
        end
        
        
        [chix, delta_chix] = CalMagSusp(sx_data(end/2:end,:), N);
        
        fprintf("chi_x = %.12f\n", mean(chix));
        fprintf("delta chi_x = %.12f\n", sqrt(var(chix)/num_chain));
        
        
        [chiy, delta_chiy] = CalMagSusp(sy_data(end/2:end,:), N);
        fprintf("chi_y = %.12f\n", mean(chiy));
        fprintf("delta chi_y = %.12f\n", sqrt(var(chiy)/num_chain));
        
        
        [chiz, delta_chiz] = CalMagSusp(sz_data(end/2:end,:), N);
        fprintf("chi_z = %.12f\n", mean(chiz));
        fprintf("delta chi_z = %.12f\n", sqrt(var(chiz)/num_chain));
        
        chi=chix+chiy+chiz;
        fprintf("chi(x+y+z) = %.12f\n", mean(chi));
        
        energy_data = energy_data(end/2+1:end,:);
        energy = mean(energy_data)/N;
        fprintf("e = %.12f\n", mean(energy));
        fprintf("delta e = %.12f\n", sqrt(var(energy)/num_chain));
        
        c = (mean(energy_data.^2) - mean(energy_data).^2) * beta^2 /N;
        fprintf("delta c = %.12f\n", sqrt(var(c)/num_chain));
        
        stiffness = mean(stiffness_data(end/2+1:end,:));
        fprintf("stiffness = %.12f\n", mean(stiffness));
        fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        
        save(direct_data_file, 'chix','chiy','chiz','energy','c','stiffness');
    end
    
    chi_x_set(j) = mean(chix);
    chi_y_set(j) = mean(chiy);
    chi_z_set(j) = mean(chiz);
    c_set(j) = mean(c);
    energy_set(j) = mean(energy);
    stiffness_set(j) = mean(stiffness);
end

% h1=plot(T_set, chi_x_set,'-x');hold on;
% h2=plot(T_set, chi_y_set,'-.');hold on;
% h3=plot(T_set, chi_z_set,'-o');hold on;


% l=legend([h1,h2,h3],'$\chi_x$', '$\chi_y$', '$\chi_z$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',36);
% set(l,'Location','SouthWest');
% 
% 
% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% xlabel('$T/K$','Interpreter','latex');
% ylabel('Magnetic susceptibility $\chi$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);
% set(gca, 'Ylim',[0,inf]);


% figure;
% plot(T_set, c_set,'-o');hold on;
% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% xlabel('$T/K$','Interpreter','latex');
% ylabel('specific heat $C$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);
% 
% 
% figure;
if(stiffness_set(1) < 1)
    plot(T_set, stiffness_set*N, '-o');hold on; %for the wrong code 2020 Jun 10
else
    plot(T_set, stiffness_set, '-o');hold on;
end
plot(T_set, 2./beta_set/pi,'-.');hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T/K$','Interpreter','latex');
ylabel('spin stiffness $\rho$(meV)','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
set(gca, 'Ylim',[0,inf]);
end

function [chi, delta_chi] = CalMagSusp(smu_data, N)
% smu_data is matrix
smu_square = smu_data .* smu_data;
chi_list = (mean(smu_square) - mean(smu_data).^2)/N;
chi = mean(chi_list);
delta_chi = 0;
% delta_chi = sqrt( var(chi_list)/numel(chi_list));
end
