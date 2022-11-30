%{
Sep 28, ../../data
1: energy; 2: specific heat; 3~5: Q=0 magnetic susceptibilities;
6~14: AF magnetization (Q1xyz, Q2xyz, Q3xyz);
15~17: binder ratio of magnetization
%}
%, 48, 64, 80
%{
Sep 28, ../../data2
1: energy; 2: specific heat; 3~11: anti-ferromagnetic
susceptibility(Sx[Q1,Q2,Q3], Sy[Q1,Q2,Q3], Sz[Q1,Q2,Q3]);
12~20: AF magnetization (Q1xyz, Q2xyz, Q3xyz);
21~23: binder ratio of magnetization
24: C3 symmetry order parameter
25: C3 symmetry order parameter susceptibilities
26: complex order parameter
27: binder ratio of complex order parameter
%}
L_set = [ 32, 48,64, 80];

% beta_set = [ 0.2041  0.1538     0.1235 ...
%     0.1031    0.0952      0.0885      0.0826  0.0810 ...
%     0.0800 0.0790 0.0780 ...
%     0.0775  ...
%     0.077 0.0765 0.076 0.0755 ...
%     0.0752  0.0745 0.074 0.0735  0.0730  0.0709   0.0690 ...
%     0.0671           0.0606      0.0578];


beta_set = [ 0.2041  0.1538      ...
    0.0810 ...
    0.0800 0.0790 0.0780 ...
    0.0775  ...
    0.077 0.0765 0.076 0.0755 ...
    0.0752  0.0745 0.074 0.0735  0.0730  0.0709   0.0690 ...
    0.0671           0.0606      0.0578];

% 
% beta_set = [ 0.2041      ...
%     0.0810 ...
%     0.0800 ...
%     0.0775  ...
%     0.077   0.0755 ...
%     0.0752  0.0745 0.074 0.0735  0.0730  0.0709   0.0690 ...
%     0.0671           0.0606      0.0578];


J1zz = 5.31;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.113;
num_chain = 12;

prefix = '../../data2/';
eVtoK_const = 11.606;


energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
chix_set = zeros(numel(L_set), numel(beta_set));
chiy_set = zeros(numel(L_set), numel(beta_set));
chiz_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
mxy2_set  = zeros(numel(L_set), numel(beta_set));
binder_ratio1_set = zeros(numel(L_set), numel(beta_set));
binder_ratio2_set = zeros(numel(L_set), numel(beta_set));
binder_ratio3_set = zeros(numel(L_set), numel(beta_set));
binder_ratiom_set = zeros(numel(L_set), numel(beta_set));
bc3_set  = zeros(numel(L_set), numel(beta_set));
chic3_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = 2 * L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
        postfix = ['hei-rank',num2str(0), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
            'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
            'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
        test_load = load([prefix, 'summary', postfix]);
        data_type_size = numel(test_load);%15
        averaged_data=zeros(data_type_size, num_chain);
        
        
        for i = 0:num_chain-1
            postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
                'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
                'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            data = load([prefix, 'summary', postfix]);
            averaged_data(:,i + 1) = data(1:data_type_size);
        end
        energy = averaged_data(1,:);
        
        c =  averaged_data(2,:);
        
        %         chix = averaged_data(3,:);
        %         chiy = averaged_data(4,:);
        %         chiz = averaged_data(5,:);
        chix = sum(averaged_data([3,6,9]))/N* beta/eVtoK_const;
        chiy = sum(averaged_data([4,7,10]))/N* beta/eVtoK_const;
        chiz = sum(averaged_data([5,8,11]))/N* beta/eVtoK_const;
        
        %         msquare = sum(averaged_data([6,7,9,10,12,13],:))/3 ;
        %         msquare = sum(averaged_data(6:14,:))/3 ;
        msquare = sum(averaged_data([12,13,15,16,18,19],:))/3 ;
        %         fprintf("msquare = %.12f\n", mean(msquare));
        %         fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        mxy2 = averaged_data(26,:);
        binder_ratio1 = averaged_data(21,:);
        binder_ratio2 = averaged_data(22,:);
        binder_ratio3 = averaged_data(23,:);
        bc3 = averaged_data(24,:)/N^3* beta/eVtoK_const;
        chic3 = averaged_data(25,:)/N^3 * beta/eVtoK_const;
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        chix_set(system_ind, beta_ind) = mean(chix);
        chiy_set(system_ind, beta_ind) = mean(chiy);
        chiz_set(system_ind, beta_ind) = mean(chiz);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        mxy2_set(system_ind, beta_ind) = mean(mxy2);
        
        binder_ratio1_set(system_ind, beta_ind) = mean(binder_ratio1);
        binder_ratio2_set(system_ind, beta_ind) = mean(binder_ratio2);
        binder_ratio3_set(system_ind, beta_ind) = mean(binder_ratio3);
        bc3_set(system_ind, beta_ind) = mean(bc3);
        chic3_set(system_ind, beta_ind) = mean(chic3);
        if strcmp(prefix, '../../data2/')
            binder_ratiom = averaged_data(27,:);
            binder_ratiom_set(system_ind, beta_ind) = mean(binder_ratiom);
        end
    end
end

T_set = eVtoK_const./beta_set;

% plot(T_set, energy_set,'-^');hold on;
% plot(T_set, specific_heat_set, '-o');hold on;
% plot(T_set(10:end), (binder_ratio1_set(:,10:end) + binder_ratio2_set(:,10:end) + binder_ratio3_set(:,10:end))/3,'-o');hold on;
% start_num = 8;
% plot(T_set(start_num:end), (binder_ratio1_set(:,start_num:end) + binder_ratio2_set(:,start_num:end) + binder_ratio3_set(:,start_num:end))/3,'-o');hold on;
% plot(T_set(start_num:end), 5/3*ones(size(T_set(start_num:end))),'-.');hold on;
% plot(T_set,binder_ratiom_set,'-o');hold on;
% plot(T_set, binder_ratio2_set,'-o');hold on;
% plot(T_set, binder_ratio3_set,'-o');hold on;
% plot(T_set, msquare_set,'-o');hold on;
% plot(T_set, mxy2_set,'-o'); hold on;
% plot(T_set, bc3_set, '-x');hold on;
% plot(T_set, chic3_set, '-o');hold on;
% plot(T_set, chix_set, '-o');hold on;
% plot(T_set, chiy_set, '-o');hold on;
% plot(T_set, chiz_set, '-o');hold on;
 plot(T_set, chix_set+chiy_set+chiz_set, '-o');hold on;

% eta_set = zeros(1, numel(T_set));
% % loglog(L_set, mxy2_set','x');hold on;
% for i = 1:(numel(T_set))
%     fit_x = msquare_set(:,i)';
%     p = fit(log(L_set)',log(fit_x'),'poly1');
%     fprintf('T = %.5f, \t beta = %.5f,\t eta=%.5f\n',T_set(i), eVtoK_const/T_set(i), -p.p1);
%     x = L_set(1):0.5:L_set(end);
%     %     fl=loglog(x,exp(p.p2)*x.^p.p1,'-.'); hold on;
%     eta_set(i) = -p.p1;
% end
% plot(T_set, eta_set,'-x');hold on;
% plot(T_set, 0.25*ones(1, numel(T_set)),'-.');
% plot(T_set, 1/9*ones(1, numel(T_set)),'-.');
% 

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
% ylabel('specific heat','Interpreter','latex');
% ylabel('$M^2(Q)$','Interpreter','latex');
% ylabel('complex order parameter','Interpreter','latex');
% ylabel('binder ratio $\langle M^4\rangle / \langle M^2\rangle^2$','Interpreter','latex');
ylabel('binder ratio $\langle m^4\rangle / \langle m^2\rangle^2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

l=legend('$L=32$','$48$','$64$','$80$');%
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');
%  set(gca, 'Ylim',[0,inf]);
%
