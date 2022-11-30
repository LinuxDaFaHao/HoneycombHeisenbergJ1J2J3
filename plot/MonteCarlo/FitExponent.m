% Tc = (8.88 + 8.91)/2;
Tc = 8.842;
L_set = [32,48,64,80];
T_max = [ 9.2500, 9.0370, 8.9800,8.93];
t_max = T_max - Tc;
loglog(L_set, t_max,'o');hold on;

fit_x = L_set;
fit_y = t_max;
p = fit(log(fit_x'),log(fit_y'),'poly1');
nu = -1/p.p1;
fprintf('nu=%.5f\n',nu);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('$t_max$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

figure;
C_max = [2.505,3.115,3.81,4.33];
loglog(L_set, C_max,'o');hold on;


fit_x = L_set;
fit_y = C_max;
p = fit(log(fit_x'),log(fit_y'),'poly1');
alpha = p.p1 * nu;
fprintf('alpha/nu=%.5f\n',p.p1);
fprintf('alpha=%.5f\n',alpha);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('$C_max$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);



figure;
chi_max = [5.3,10.95,18.4,27.9];
loglog(L_set, chi_max,'o');hold on;


fit_x = L_set;
fit_y = chi_max;
p = fit(log(fit_x'),log(fit_y'),'poly1');
gamma = p.p1 * nu;
fprintf('gamma/nu=%.5f\n',p.p1);
fprintf('gamma=%.5f\n',gamma);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L$','Interpreter','latex');
ylabel('$\chi_{max}$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
