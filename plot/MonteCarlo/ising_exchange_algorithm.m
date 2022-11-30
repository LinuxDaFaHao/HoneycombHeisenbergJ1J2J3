% L_set = [32,48,48,64,80,128];
% beta_min_set = [0.05,0.05,0.1,0.05,0.05,0.05];
% beta_max_set = [0.125,0.125,0.113,0.125,0.125,0.125];
L_set = [80,128,128];
beta_min_set = [0.05,0.05, 0.109,0.1,0.05,0.05,0.05];
beta_max_set = [0.125,0.125,0.112,0.113,0.125,0.125,0.125];

%Note:
% Tc = 8.88 ~ 8.91, beta=0.1126~0.1122

geometry = 'Honeycomb';
J1 = 1;
J2 = 0;
J3 = -5;
thread_num_set = [96,96,48,96,96,96];
prefix = '../../data/';

eVtoK_const = 1;

for i = 1:numel(L_set)
    L = L_set(i);
    N = L*L*2;
    beta_min = beta_min_set(i);
    beta_max = beta_max_set(i);
    thread_num = thread_num_set(i);
    postfix = ['ising', geometry, 'J1', num2str(J1,'%.6f'),...
        'J2', num2str(J2,'%.6f'),  'J3', num2str(J3,'%.6f'),...
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
    af_magnetization_1 = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_2 = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_3 = fread(file_id, thread_num, 'double');
    fread(file_id, 1, 'char');
    af_magnetization_chi = fread(file_id, thread_num, 'double');
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
    %
    fclose(file_id);
%     plot(T_set, (af_magnetization_1+af_magnetization_2+af_magnetization_3)/3,'-o');hold on;
    [T_set,I] = sort(T_set);
    plot(T_set, af_magnetization_binder_ratio(I), '-o');hold on;

%      plot(T_set, af_magnetization_binder_ratio, '-o');hold on;
%       plot(T_set, specific_heat, '-o');hold on;
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