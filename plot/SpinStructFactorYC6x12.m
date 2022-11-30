geometry = 'YC';
Ly = 6;
Lx = 12;
J1 = -1;
J2 = -0.04;
J3 = 5.3;
Dzz = 0.02;
Db = 15000;
FileNamePrefix = '../data/';
%template: HoneyHeiYC4x8J1-1J2-0.04J35.3Dzz0.02D15000.json
FileNamePostfix = ['HoneyHei',geometry, num2str(Ly), 'x', num2str(Lx),...
    'J1', num2str(J1), 'J2', num2str(J2),  'J3', num2str(J3),...
    'Dzz', num2str(Dzz),'D', num2str(Db),'.json'];

% SzData = jsondecode(fileread([FileNamePrefix,'sz',FileNamePostfix]));
% SzSzCorrelationData = jsondecode(fileread([FileNamePrefix,'zzsf',FileNamePostfix]));
% SpSmCorrelationData = jsondecode(fileread([FileNamePrefix,'pmsf',FileNamePostfix]));
% SmSpCorrelationData = jsondecode(fileread([FileNamePrefix,'mpsf',FileNamePostfix]));

datanum_of_structfactor = numel(SzSzCorrelationData);
sdots_data = zeros(1, datanum_of_structfactor);
szsz_correlation = zeros(1, datanum_of_structfactor);
spsm_correlation = zeros(1, datanum_of_structfactor);
smsp_correlation = zeros(1, datanum_of_structfactor);

kx_set = -4*pi/3:0.02:4*pi/3; % exchange the range of kx and ky if XC cylinder
ky_set = -2*pi/sqrt(3):0.02:2*pi/sqrt(3);
kx_num = numel(kx_set);
ky_num = numel(ky_set);
struct_factor = zeros(ky_num, kx_num);
aa_datanum_of_structfactor = 0;
for i = 1:datanum_of_structfactor
    site1 = SzSzCorrelationData{i}{1}(1);%count from 0
    site2 = SzSzCorrelationData{i}{1}(2);
    [x1,y1] = HoneycombYCCylinderSiteInd2XYCoor(site1, Ly);
    [x2,y2] = HoneycombYCCylinderSiteInd2XYCoor(site2, Ly);
    
    y1_idx = mod(site1, Ly);
    x1_idx = floor(site1 / Ly);
    a1_sublattice = mod(x1_idx + y1_idx, 2);
    
    y2_idx = mod(site2, Ly);
    x2_idx = floor(site2 / Ly);
    a2_sublattice = mod(x2_idx + y2_idx, 2);
    if(a1_sublattice ~= 0 || a2_sublattice ~=0)
        continue
    end
    
    dx = x1 - x2;
    dy = y1 - y2;
    sz_site1 = SzData(site1, 2);
    sz_site2 = SzData(site2, 2);
    szsz_correlation(i) = SzSzCorrelationData{i}{2} ;%- sz_site1 * sz_site2;
    spsm_correlation(i) = SpSmCorrelationData{i}{2};
    smsp_correlation(i) = SmSpCorrelationData{i}{2};
    sdots_data(i) = szsz_correlation(i) + 1/2 * ( spsm_correlation(i) + smsp_correlation(i) );
    
    struct_factor = struct_factor +  sdots_data(i) * cos( (kx_set *dx + ky_set' * dy));
    aa_datanum_of_structfactor = aa_datanum_of_structfactor+1;
    
end
struct_factor = struct_factor / aa_datanum_of_structfactor;
for kx_idx = 1:kx_num
    kx = kx_set(kx_idx);    
    for ky_idx = 1:ky_num
      ky = ky_set(ky_idx);
      if( abs(ky/(4*pi/sqrt(3))) + abs(kx/(4*pi/3)) > 1 )
          struct_factor(ky_idx, kx_idx) = NaN;
      end
    end
end

imagesc(kx_set, ky_set,struct_factor);hold on;
colorbar

T1=text(0,0,'$\Gamma$');
set(T1,'Interpreter','latex');set(T1,'Fontsize',32);

% T2=text(0, 

set(gca,'YDir','normal');
axis equal;
axis off;

set(gca,'fontsize',32);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$M$','Interpreter','latex');
ylabel('$K$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',32); 
set(get(gca,'YLabel'),'FontSize',32); 
set(gca,'xtick',[]); 
set(gca,'ytick',[]); 
set(gca,'linewidth',1.5);
set(gcf,'position',[1000,1000, 2*300, sqrt(3)*300]);