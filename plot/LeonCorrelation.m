figure;
L_set = [ 32, 48,  64,80,96];

%beta_set = [    0.1124       0.1031     0.0952          0.0855       0.0800       0.0752     0.0709       0.0671       0.0637       0.0606       0.0578      ];
  beta_set = [0.0800 0.0752 0.075 0.0745 0.074 0.0735 ...
      0.0734 0.07335 0.0733 0.07325 0.0732 0.0731 ...
      0.073 0.0725 0.072 0.0715 0.071 0.0709];
%   beta_set = [ ];
% beta_set = [];
% beta_set = [  0.03 0.02 0.01 0.005 0.002 ];
% beta_set = sort(beta_set);
% beta_set = [0.054 0.053 0.052 0.051 0.050 0.049 0.048 0.047 0.046 0.043 0.041 0.039 0.038 0.037 0.036 0.034];
J1zz = 5.3;
J2zz = 0.2;
J3zz = -28.00002;
Dzz = -0.3;
num_chain = 5;

prefix = '../data/';
save_data_prefix = './plot_data/';
eVtoK_const = 11.606;

stiffness_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
C2x_set = zeros(numel(L_set), numel(beta_set));
C2y_set = zeros(numel(L_set), numel(beta_set));
C2z_set = zeros(numel(L_set), numel(beta_set));
C4x_set = zeros(numel(L_set), numel(beta_set));
C4y_set = zeros(numel(L_set), numel(beta_set));
C4z_set = zeros(numel(L_set), numel(beta_set));

msquare_set = zeros(numel(L_set), numel(beta_set));
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
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
%         test_load = load([prefix, 'summary', postfix]);
        data_type_size = 16;%numel(test_load);
        averaged_data=zeros(data_type_size, num_chain);

        
        direct_data_file = [save_data_prefix,'/hei','Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L),'.mat'];
        if exist(direct_data_file,'file') && 0
            load(direct_data_file);
        else
            
            for i = 0:num_chain-1
                postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
                    'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                    'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
                data = load([prefix, 'summary', postfix]);
                averaged_data(:,i + 1) = data(1:data_type_size);
            end
            
            C2x = averaged_data(11,:);
            C4x = averaged_data(12,:);
            
            fprintf("C2_x = %.12f\n", mean(C2x));
            fprintf("delta C2_x = %.12f\n", sqrt(var(C2x)/num_chain));
            
            C2y = averaged_data(13,:);
            C4y = averaged_data(14,:);
            fprintf("chi_y = %.12f\n", mean(C2y));
            fprintf("delta chi_y = %.12f\n", sqrt(var(C2y)/num_chain));
            
            
            C2z = averaged_data(15,:);
            C4z = averaged_data(16,:);
%             fprintf("chi_z = %.12f\n", mean(chiz));
%             fprintf("delta chi_z = %.12f\n", sqrt(var(chiz)/num_chain));
%             
%             chi=chix+chiy+chiz;
%             fprintf("chi(x+y+z) = %.12f\n", mean(chi));
% 
%             
%             energy = averaged_data(1,:);
%             fprintf("e = %.12f\n", mean(energy));
%             fprintf("delta e = %.12f\n", sqrt(var(energy)/num_chain));
%             
            c =  averaged_data(2,:);
%             fprintf("c = %.12f\n", mean(c));
%             fprintf("delta c = %.12f\n", sqrt(var(c)/num_chain));
%             
            stiffness = averaged_data(6,:);
            stiffness = stiffness(stiffness>0);
%             fprintf("stiffness = %.12f\n", mean(stiffness));
%             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
            
            msquare = averaged_data(8,:);
            m = max(msquare);
            
            msquare = msquare(msquare > m/2)/3;
            fprintf("msquare = %.12f\n", mean(msquare));
            fprintf("n = %d\n", numel(msquare));
            fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
            
            
%             save(direct_data_file, 'chix','chiy','chiz','energy','c','stiffness');
        end
        
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        C2x_set(system_ind, beta_ind) = mean(C2x);
        C2y_set(system_ind, beta_ind) = mean(C2y);
        C2z_set(system_ind, beta_ind) = mean(C2z);
        C4x_set(system_ind, beta_ind) = mean(C4x);
        C4y_set(system_ind, beta_ind) = mean(C4y);
        C4z_set(system_ind, beta_ind) = mean(C4z);
        msquare_set(system_ind, beta_ind) = mean(msquare)  ;
    end
end
% save('StiffnessD3dJ328Dzz0_113.mat','stiffness_set');
% load('StiffnessD3dJ328Dzz0_113.mat','stiffness_set');
T_set = eVtoK_const./beta_set;

% plot(T_set, msquare_set,'-o');hold on;

%  plot(T_set, C2x_set, '-o');hold on;
% plot(T_set, C2y_set, '-o');hold on;
% plot(T_set, C2z_set, '-o');hold on;
% plot(T_set, C4x_set, '-o');hold on;
% plot(T_set, C4y_set, '-o');hold on;
% plot(T_set, C4z_set, '-o');hold on;

%  plot(T_set, (C2x_set+C2y_set+C2z_set)./(C4x_set+C4y_set+C4z_set),'-x');hold on;
 plot(T_set, (C2x_set)./(C4x_set),'-x');hold on;


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('specific heat','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Ylim',[0,inf]);
