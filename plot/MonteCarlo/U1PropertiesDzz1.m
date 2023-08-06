L_set = [32,48,64];
beta_min_set = [0.055,0.06,0.06];
beta_max_set = [0.1,0.1,0.1];


geometry = 'Honeycomb';
J1 = 5.3;
J2 = 0.2;
J3 = -28;
Dzz = -1;
thread_num_set = [96,96,96];
% marker_color = ['b','r','g','g','g','y'];
marker_color = {[0.0,0.45,0.74],[0.85,0.33,0.10],[0.47,0.67,0.19], ...
    [0.47,0.67,0.19],...
    [0.47,0.67,0.19],...
    [0.93,0.69,0.13],...
    [0.30,0.75,0.93]};
prefix = '../../data/';

eVtoK_const = 11.606;

for i = 1:numel(L_set)
    L = L_set(i);
    N = L*L*2;
    beta_min = beta_min_set(i);
    beta_max = beta_max_set(i);
    thread_num = thread_num_set(i);
    postfix = ['hei', geometry, 'J1zz', num2str(J1,'%.6f'),...
        'J2zz', num2str(J2,'%.6f'),  'J3zz', num2str(J3,'%.6f'),...
        'Dzz', num2str(Dzz,'%.6f'), ...
        'beta_max',num2str(beta_max,'%.6f'),...
        'beta_min',num2str(beta_min,'%.6f'),'L', num2str(L)];
    file_id = fopen([prefix,'results',postfix],'rb');
    beta_set = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    T_set = eVtoK_const./beta_set;
    energy = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    specific_heat = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    stiffness = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    af_magnetization_1x = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_1y = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_1z = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_2x = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_2y = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_2z = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    af_magnetization_3x = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_3y = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_3z = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    af_magnetization_chix = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_chiy = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_chiz = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    af_magnetization_1_binder_ratio = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_2_binder_ratio = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_3_binder_ratio = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_binder_ratio = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_chi_b = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    c3_order_parameter = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    c3_order_parameter_chi = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    complex_order_parameter = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    complex_order_parameter_binder_ratio = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    binder_ratio_c3 = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    
    fclose(file_id);
%         [T_set,I] = sort(T_set);
%         plot(T_set, af_magnetization_1x(I) ...
%              + af_magnetization_1y(I) ...
%              + af_magnetization_2x(I) ...
%              + af_magnetization_2y(I) ...
%              + af_magnetization_3x(I) ...
%              + af_magnetization_3y(I) ,'^');hold on;
    
    %     plot(T_set, af_magnetization_chix(I) ...
    %         + af_magnetization_chiy(I) ...
    %         + af_magnetization_chiz(I), '-o');hold on;
    %
%     h(i) = semilogy(T_set, af_magnetization_chi_b,'o','color',marker_color{i});hold on;
%     plot(T_set, af_magnetization_binder_ratio(I), '-o');hold on;
%     plot(T_set, specific_heat, 'o');hold on;
%      h(i) = semilogy(T_set, c3_order_parameter_chi, 'x','color',marker_color{i});hold on;
%       h(i) = plot(T_set, c3_order_parameter, 'x','color',marker_color{i});hold on;
%      h(i) = plot(T_set, binder_ratio_c3, 'x','color',marker_color{i});hold on;
%     plot(T_set, energy, 'o');hold on;
 plot(T_set, stiffness, '-o');hold on;
%       plot(T_set, specific_heat, 'o');hold on;
%     beta2 = beta_set(1:end-1) + diff(beta_set)/2;
%     C = -beta2.^2 .* (diff(energy)./diff(beta_set));
%     plot(1./beta2, C,'-x');hold on;
end
plot(T_set, 2./beta_set/pi,'-.');hold on;
% l=legend([h(1),h(2),h(3),h(6),h(7)],{'16','24','32','48','64'});
l=legend('$L=32$','48','64','$\rho = \frac{2T}{\pi}$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% xlabel('$T(K)$','Interpreter','latex');
% ylabel('AF susceptibility','Interpreter','latex');
% ylabel('$C_3$ order parameter','Interpreter','latex');
% ylabel('binder ratio of $C_3$ order parameter','Interpreter','latex');
ylabel('stiffness $\rho$(meV)','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);



Tc_set = [168.61, 164.96, 162.96];%K
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
ylabel('$T_{\mathrm{BKT}}/K$','Interpreter','latex');
ylabel('$T^*(K)$','Interpreter','latex');
xlabel('$(\ln L)^2$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

set(gca, 'Xlim',[0,inf]);

