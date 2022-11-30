L_set = [32,32,32];
beta_min_set = [0.06,0.075,0.075,0.05,0.05,0.05];
beta_max_set = [1.1,0.2,0.1,0.125,0.125,0.125];


geometry = 'Honeycomb';
J1 = 5.3;
J2 = 0.2;
J3 = -28;
Dzz = -0.113;
thread_num_set = [96,96,48];
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
    
%     complex_order_parameter = fread(file_id, thread_num, 'double');
%     fread(file_id, 1, 'char');
%     complex_order_parameter_binder_ratio = fread(file_id, thread_num, 'double');
%     fread(file_id, 1, 'char');
    %
    fclose(file_id);
%     plot(T_set, (af_magnetization_1+af_magnetization_2+af_magnetization_3)/3,'-o');hold on;
    [T_set,I] = sort(T_set);
%     plot(T_set, af_magnetization_1x(I) ...
%          + af_magnetization_1y(I) ...
%          + af_magnetization_2x(I) ...
%          + af_magnetization_2y(I) ...
%          + af_magnetization_3x(I) ...
%          + af_magnetization_3y(I) ,'^');hold on;

%     plot(T_set, af_magnetization_chix(I) ...
%         + af_magnetization_chiy(I) ...
%         + af_magnetization_chiz(I), '-o');hold on;

%     plot(T_set, af_magnetization_chi_b,'o');hold on;
%       plot(T_set, af_magnetization_binder_ratio(I), '-o');hold on;
%        plot(T_set, specific_heat(I), 'o');hold on;
     plot(T_set, c3_order_parameter_chi(I), 'x');hold on;
    %     
%     beta2 = beta_set(1:end-1) + diff(beta_set)/2;
%     C = -beta2.^2 .* (diff(energy)./diff(beta_set));
%     plot(1./beta2, C,'-x');hold on;
end

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T(K)$','Interpreter','latex');
ylabel('specific heat','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);