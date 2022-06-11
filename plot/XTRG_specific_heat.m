Geometry = 'YC';
Ly = 6;
Lx = 6;
J1 = -1;
J2 = -0.04;
J3 = 5.3;
Dzz = 0.02;
D = 1500;
eVtoK_const = 11.606;
J1_in_eV = 5.21;
spline_density = 5; % integer, largerl than or equal to 1. 

% FileNamePostfix=['HoneyHei',Geometry, num2str(Ly),'x',num2str(Lx),'J1',num2str(J1),'J2',num2str(J2),'J3',num2str(J3),'D',num2str(D),'.json'];
FileNamePostfix=['HoneyHei',Geometry, num2str(Ly),'x',num2str(Lx),'J1',num2str(J1),'J2',num2str(J2),'J3',num2str(J3),'Dzz',num2str(Dzz),'D',num2str(D),'.json'];

FreeEnergyData = jsondecode(fileread(['../data/free_energy',FileNamePostfix]));

%  FreeEnergyData = jsondecode(fileread('../data/free_energyXC4x4D600.json'));



beta = FreeEnergyData(:,1);
Fn = FreeEnergyData(:,2);
log_beta = log(beta);

dlog_beta = (log_beta(2)-log_beta(1))/spline_density;
log_beta_set = log_beta(1):dlog_beta:log_beta(end);
Fn_set = spline(log_beta, Fn, log_beta_set);
beta_set = exp(log_beta_set);
T_set = 1./beta_set;
mlnZ = Fn_set.*beta_set;

En = diff(mlnZ)./diff(beta_set);
beta_en = beta_set(1:end-1)+diff(beta_set)/2;
T_en = 1./beta_en;

C = diff(En)./diff(T_en);
T_specific_heat = T_en(1:end-1) + diff(T_en)/2;
% semilogx(T_specific_heat, C,'-o');hold on;

T = 1./beta * eVtoK_const * J1_in_eV;
C_original = spline(T_specific_heat, C, 1./beta);
semilogx( T, C_original,'-x');hold on;



% set(gca, 'Xlim',[0,100]);
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T$','Interpreter','latex');
ylabel('specific heat $C$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


