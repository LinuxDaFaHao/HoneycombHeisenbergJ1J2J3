geometry = 'YC';
Ly = 6;
Lx = 12;
J1 = -5.3;
J2 = -0.2;
J3 = 38;
Dzz = 0.113;
Db = 17000;
FileNamePrefix = '../data/';
%template: HoneyHeiYC4x8J1-1J2-0.04J35.3Dzz0.02D15000.json
FileNamePostfix = ['HoneyHei',geometry, num2str(Ly), 'x', num2str(Lx),...
    'J1', num2str(J1), 'J2', num2str(J2),  'J3', num2str(J3),...
    'Dzz', num2str(Dzz),'D', num2str(Db),'.json'];

% draw parameter
lattice_linewidth  =  3;
spin_arrow_scale = 0.6;
reference_point_color = [148 0 211]/256;



SzData = jsondecode(fileread([FileNamePrefix,'sz',FileNamePostfix]));
SzSzCorrelationData = jsondecode(fileread([FileNamePrefix,'zzsf',FileNamePostfix]));
SpSmCorrelationData = jsondecode(fileread([FileNamePrefix,'pmsf',FileNamePostfix]));
SmSpCorrelationData = jsondecode(fileread([FileNamePrefix,'mpsf',FileNamePostfix]));

N = Lx*Ly*2;


reference_site = N/2;
num_of_center_correlation = 0;
for i = 1:numel(SzSzCorrelationData)
    if(SzSzCorrelationData{i}{1}(2) == reference_site)
        num_of_center_correlation = num_of_center_correlation + 1;
    elseif( SzSzCorrelationData{i}{1}(1) == reference_site)
        num_of_center_correlation = num_of_center_correlation + 1;
    end
end

site2x_set = zeros(1, num_of_center_correlation);
site2y_set = zeros(1, num_of_center_correlation);
sxsx_correlation = zeros(1, num_of_center_correlation);
sysy_correlation = zeros(1, num_of_center_correlation);
j = 0;
for i = 1:numel(SzSzCorrelationData)
    
    if(SzSzCorrelationData{i}{1}(2) == reference_site)
        j = j + 1;
        [site2x_set(j), site2y_set(j)] = HoneycombYCCylinderSiteInd2XYCoor(SzSzCorrelationData{i}{1}(1),Ly);
        sxsx_correlation(j) = SpSmCorrelationData{i}{2}/2;
%         sysy_correlation(j) = SpSmCorrelationData{i}{2}/2;
    elseif( SzSzCorrelationData{i}{1}(1) == reference_site)
        j = j + 1;
        [site2x_set(j), site2y_set(j)] = HoneycombYCCylinderSiteInd2XYCoor(SzSzCorrelationData{i}{1}(2),Ly);
        sxsx_correlation(j) = SpSmCorrelationData{i}{2}/2;
%         sysy_correlation(j) = SpSmCorrelationData{i}{2}/2;
    end
end
selectA = sxsx_correlation<0;
q1=quiver(site2x_set(selectA),site2y_set(selectA),sxsx_correlation(selectA),sysy_correlation(selectA),spin_arrow_scale,'linewidth',3);
hold on;
selectB = sxsx_correlation>0;
q2=quiver(site2x_set(selectB),site2y_set(selectB),sxsx_correlation(selectB),sysy_correlation(selectB),spin_arrow_scale,'linewidth',3);

hold on;
[refx,refy] = HoneycombYCCylinderSiteInd2XYCoor(reference_site,Ly);

plot(refx,refy,'x','markersize',25,'linewidth',5,'color',reference_point_color);hold on;



% figure
%plot lattice
x_set = zeros(1,N);
y_set = zeros(1,N);
for site_idx = 0:N-1
    [x_set(site_idx+1), y_set(site_idx+1)] = HoneycombYCCylinderSiteInd2XYCoor(site_idx, Ly);
    %    plot(x,y,'o','Color', 'k', 'MarkerFaceColor','k');hold on;
end
scatter(x_set,y_set,20,'k','o','filled');
for site_idx = 0:N-1
    
    y_idx = mod(site_idx, Ly);
    x_idx = floor(site_idx / Ly);
    x = x_idx / 2;
    a_sublattice = mod(x_idx + y_idx, 2);
    y = y_idx * sqrt(3) / 2 + (1 - a_sublattice) / (2 * sqrt(3));
    if(a_sublattice==0)
       line([x,x],[y,y+1/sqrt(3)],'color','k','linewidth',lattice_linewidth); 
       line([x,x+1/2],[y,y-1/2/sqrt(3)],'color','k','linewidth',lattice_linewidth);
    else
        line([x,x+1/2],[y,y+1/2/sqrt(3)],'color','k','linewidth',lattice_linewidth);
    end
end
% set(gca,);
axis off; % axis tight;
axis equal;