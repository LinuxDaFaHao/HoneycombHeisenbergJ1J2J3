%{
Oct 23, ../../data
1: energy; 2: specific heat;
%}
L_set = [16];


beta_set = [ 0.3800    0.3815    0.3830    0.3845    0.3860    0.3874    0.3889    0.3904    0.3919    0.3934    0.3949    0.3964    0.3979    0.3994    0.4009    0.4023    0.4038    0.4053   0.4068    0.4083    0.4098    0.4113    0.4128    0.4143    0.4157    0.4172    0.4187   0.4202    0.4217    0.4232    0.4247    0.4262    0.4277    0.4291    0.4306    0.4321  0.4336    0.4351    0.4366    0.4381    0.4396    0.4411    0.4426    0.4440    0.4455      0.4470  0.4485    0.4500];
%
J1 = 1;
J2 = 0;
J3 = 0;
D = 0.0;
num_chain = 5;

prefix = '../../data/';
eVtoK_const = 1;

specific_heat_set = zeros(numel(L_set), numel(beta_set));
specific_heat_error_set = zeros(numel(L_set), numel(beta_set));
energy_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    
    for beta_ind = 1:numel(beta_set)
        
        
        beta = beta_set(beta_ind);
        geometry = 'Square';
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['ising-rank',num2str(0), geometry, 'J1', num2str(J1,'%.6f'),...
            'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
            'D', num2str(D,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        %       test_load = load([prefix, 'summary', postfix]);
        %       data_type_size = numel(test_load);%15
        data_type_size = 2;
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
        
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        
    end
end
T_set = eVtoK_const./beta_set;



% beta2 = beta_set(1:end-1) + diff(beta_set)/2;
% C = -beta2.^2 .* (diff(energy_set)./diff(beta_set));
% plot(1./beta2, C,'-x');hold on;

plot(T_set, specific_heat_set, '-o');hold on;
errorbar(T_set, specific_heat_set, specific_heat_error_set, '-o');hold on;

% plot(T_set, energy_set,'-.');hold on;

% eta_set = zeros(1, numel(T_set));
% for i = 1:(numel(T_set))
%     fit_x = mxy2_set(:,i)';
%     p = fit(log(L_set)',log(fit_x'),'poly1');
%     fprintf('T = %.5f, \t beta = %.5f,\t eta=%.5f\n',T_set(i), eVtoK_const/T_set(i), -p.p1);
%     x = L_set(1):0.5:L_set(end);
%     %     fl=loglog(x,exp(p.p2)*x.^p.p1,'-.'); hold on;
%     eta_set(i) = -p.p1;
% end
% plot(T_set, eta_set,'-x');hold on;
% plot(T_set, 0.25*ones(1, numel(T_set)),'-.');
% plot(T_set, 1/9*ones(1, numel(T_set)),'-.');

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


