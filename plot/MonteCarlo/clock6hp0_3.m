L_set = [ 32, 64, 128, 256];
% 
beta_set = [4 3.5 3 2.7 2.3 2 1.8 1.2 1.1 1 0.9 0.8 0.7 0.6 0.5];
% 1.5 

J = 1;
hp = 0.3 *J ;
num_chain = 12;
prefix = '../../data/';
eVtoK_const = 1/J;


energy_set = zeros(numel(L_set), numel(beta_set));
specific_heat_set = zeros(numel(L_set), numel(beta_set));
msquare_set = zeros(numel(L_set), numel(beta_set));
binder_ratio_set = zeros(numel(L_set), numel(beta_set));
for system_ind = 1:numel(L_set)
    L = L_set(system_ind);
    N = L^2;
    
    for beta_ind = 1:numel(beta_set)
        beta = beta_set(beta_ind);
        
        fprintf('beta=%.6f\n', beta);
        
%         postfix = ['clock-rank',num2str(0), 'Square', 'J', num2str(J,'%.6f'),...
%             'hp', num2str(hp,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
%         test_load = load([prefix, 'summary', postfix]);
%         data_type_size = numel(test_load);%4

        data_type_size = 4;
        averaged_data=zeros(data_type_size, num_chain);
        for i = 0:num_chain-1
            postfix = ['clock-rank',num2str(i), 'Square', 'J', num2str(J,'%.6f'),...
                'hp', num2str(hp,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
            file_name = [prefix, 'summary', postfix];
            if(exist(file_name,'file'))
                data = load(file_name);
                averaged_data(:,i + 1) = data(1:data_type_size);
            else
                averaged_data(:,i + 1) = NaN;
            end
        end
        energy = averaged_data(1,:);
        %             fprintf("energy = %.12f\n", mean(energy));
        %             fprintf("delta e = %.12f\n", sqrt(var(energy)/num_chain));
        %
        c =  averaged_data(2,:);
        %             fprintf("C = %.12f\n", mean(c));
        %             fprintf("delta C = %.12f\n", sqrt(var(c)/num_chain));
        
        
        %
        %         stiffness = averaged_data(6,:);
        %             fprintf("stiffness = %.12f\n", mean(stiffness));
        %             fprintf("delta rho = %.12f\n", sqrt(var(stiffness)/num_chain));
        %          msquare = sum(averaged_data([6,7,9,10,12,13],:))/3 ;
        %         msquare = sum(averaged_data([6,7],:)) ;
        msquare = averaged_data(3,:);
        %         fprintf("msquare = %.12f\n", mean(msquare));
        %         fprintf("delta msquare = %.12f\n", sqrt(var(msquare)/numel(msquare)));
        %    fprintf("M_x2/M_y2 = %.12f\n", mean(averaged_data(8,:))/mean(averaged_data(9,:)));
        
        %         stiffness_set(system_ind, beta_ind) = mean(stiffness);
        energy_set(system_ind, beta_ind) = mean(energy);
        specific_heat_set(system_ind, beta_ind) = mean(c);
        msquare_set(system_ind, beta_ind) = mean(msquare);
        %         binder_ratio = averaged_data(4,:);
        binder_ratio = averaged_data(4,:);
        binder_ratio_set(system_ind, beta_ind) = mean(binder_ratio);
        
    end
end

T_set = eVtoK_const./beta_set;

% plot(T_set, energy_set/J,'-^');hold on;
%      plot(T_set, specific_heat_set, '-o');hold on;
% start_num = 1;
% plot(T_set(start_num:end), binder_ratio_set(:,start_num:end),'-o');hold on;

% plot(T_set, binder_ratio_set,'-o');hold on;
%  plot(T_set, msquare_set,'-o');hold on;

eta_set = zeros(1, numel(T_set));
% loglog(L_set, msquare_set','x');hold on;
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
plot(T_set, 1/9*ones(1, numel(T_set)),'-.');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T = 1/\beta$','Interpreter','latex');
%  ylabel('specific heat','Interpreter','latex');
%  ylabel('$m^2$','Interpreter','latex');
  ylabel('$\eta$','Interpreter','latex');
% ylabel('binder ratio $\langle M^4\rangle / \langle M^2\rangle^2$','Interpreter','latex');
%   ylabel('binder ratio $\langle m^4\rangle / \langle m^2\rangle^2$','Interpreter','latex');
%  ylabel('complex order parameter','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
% set(gca, 'Ylim',[0,inf]);

% l=legend('$L=32$','$64$','$128$', '$256$');%
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');
