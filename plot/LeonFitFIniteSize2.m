L_set = [32,48,64, 80, 96];

beta_set = [   0.0735   ];
%
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
        %             fprintf("C_y2 = %.12f\n", mean(C2y));
        %             fprintf("delta C_y2 = %.12f\n", sqrt(var(C2y)/num_chain));
        
        
        C2z = averaged_data(15,:);
        C4z = averaged_data(16,:);
        
        msquare = averaged_data(8,:);% + averaged_data(9,:) + averaged_data(10,:);
        %             m = max(msquare);
        %
        %             msquare = msquare(msquare > m/2);
        fprintf("msquare = %.12f\n", mean(msquare));
        fprintf("n = %d\n", numel(msquare));
        fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        
        
        
        
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        C2x_set(system_ind, beta_ind) = mean(C2x);
        C2y_set(system_ind, beta_ind) = mean(C2y);
        C2z_set(system_ind, beta_ind) = mean(C2z);
        C4x_set(system_ind, beta_ind) = mean(C4x);
        C4y_set(system_ind, beta_ind) = mean(C4y);
        C4z_set(system_ind, beta_ind) = mean(C4z);
        msquare_set(system_ind, beta_ind) = mean(msquare) ;%-14.5554/L^2  ;
    end
end
% save('StiffnessD3dJ328Dzz0_113.mat','stiffness_set');
% load('StiffnessD3dJ328Dzz0_113.mat','stiffness_set');
T_set = eVtoK_const./beta_set;

%plot(T_set, specific_heat_set, '-o');hold on;

loglog(L_set, C2x_set','-o');hold on;

p = fit(log(L_set(2:end))',log(C2x_set(2:end)),'poly1');
fprintf('eta=%.5f\n',-p.p1);
x = L_set(1):0.5:L_set(end)+2;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');

T1=text(50,1e-3,['$\eta=$',num2str(-p.p1)]);
set(T1,'Interpreter','latex');set(T1,'Fontsize',32);


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('$M^2(Q)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Ylim',[0,inf]);
