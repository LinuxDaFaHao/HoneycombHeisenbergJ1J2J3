L_set = [  48, 64,128];

% beta_set = [ 0.2857 0.2326   0.1887      0.1639       0.1370        0.1176      0.1031      0.0952    0.0917    0.0885    0.0855    0.0826    0.0800   0.0775   0.0752   0.0730  0.0709    0.0690    0.0671    0.0654    0.0637    0.0621    0.0606    0.0592    0.0578   0.0565    0.0552];
% beta_set = [  beta_set,  0.0770 0.0760  0.0745 0.0740  ]; % L = 32, 48,
% 64, 72,100. don't delete
% beta_set = [0.2326  0.1887   0.1639   0.1449     0.1299    0.1176    0.1075    0.0990  0.0917   0.0855   0.0800    0.0775  0.0770    0.0760    0.0752    0.0745    0.0740    0.0730    0.0709    0.0690    0.0671    0.0654    0.0637];
beta_set = [   0.1639    0.1299    0.1176    0.1075    0.0990  0.0917   0.0855   0.0800    0.0775  0.0770    0.0760    0.0752    0.0745    0.0740    0.0730    0.0709    0.0690    0.0671    0.0654    0.0637];

J1zz = 5.3;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.113;
num_chain = 12;

prefix = '../../data/';
save_data_prefix = './plot_data/';
eVtoK_const = 11.606;

stiffness_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
chiz_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));


for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);
        averaged_data=zeros(data_type_size, num_chain);
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            temp = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = temp(1:data_type_size);
        end
        energy = averaged_data(1,:);
        c =  averaged_data(2,:);
        chix = averaged_data(3,:);
        chiy = averaged_data(4,:);
        chiz = averaged_data(5,:);
        stiffness = averaged_data(6,:);
        fprintf("stiffness = %.12f\n", mean(stiffness));
        fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        
        %         msquare = sum(averaged_data(7:15,:))/3 ;
        % msquare = sum(averaged_data([9,12,15],:))/3 ;
        msquare = sum(averaged_data([7,8,10,11,13,14],:))/3 ;
        stiffness_set(system_ind, beta_ind) = mean(stiffness);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
    end
end
T_set = eVtoK_const./beta_set;
% plot(T_set, specific_heat_set, '-o');hold on;
% loglog(L_set, msquare_set','-o');hold on;
eta_set = zeros(1, numel(T_set));
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

% h = plot(T_set, stiffness_set, '-o');hold on;
% h2 = plot(T_set, 2./beta_set/pi,'-.'); hold on;
% % l=legend([h0,h1,h2],{'$\rho_s = 8T/\pi$','$\rho_s = 4T/\pi$','$\rho_s = 2T/\pi$'});
% l=legend([h2],{'$\rho_s = \frac{2T}{\pi}$'});
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',36);
% set(l,'Location','SouthWest');

% T=text(50,10,['$L=32, 48, 64, 72, 100$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
% ylabel('spin stiffness $\rho_s$(meV)','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');

set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
% set(gca, 'Ylim',[0,inf]);
% set(gca, 'Xlim',[0,200]);

