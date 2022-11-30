L_set = [32,64,128,256];

beta_set = [   0.02042   ];

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
        m = max(msquare);
        msquare = msquare(msquare > m/2)/3;
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
        msquare_set(system_ind, beta_ind) = mean(msquare)  ;
        mxy2_set(system_ind, beta_ind) = mean(mxy2)  ;
        %         C2_set(system_ind, beta_ind) = mean(correlation2);
        %         C4_set(system_ind, beta_ind) = mean(correlation4);
        
    end
end
T_set = eVtoK_const./beta_set;


%plot(T_set, specific_heat_set, '-o');hold on;

loglog(L_set, mxy2_set','-o');hold on;

p = fit(log(L_set(2:end))',log(mxy2_set(2:end)),'poly1');
fprintf('eta=%.5f\n',-p.p1);
x = L_set(1):0.5:L_set(end)+2;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');

T1=text(100,0.38,['$\eta=$',num2str(-p.p1)]);
set(T1,'Interpreter','latex');set(T1,'Fontsize',32);


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('$m^2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%  set(gca, 'Ylim',[0,inf]);
% 