function [] = PlotYCHoneycombLattice(Lx, Ly)
N = 2 * Lx * Ly;
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
       line([x,x],[y,y+1/sqrt(3)]); 
       line([x,x+1/2],[y,y-1/2/sqrt(3)]);
    else
        line([x,x+1/2],[y,y+1/2/sqrt(3)]);
    end
end
axis off;
% axis tight;
axis equal;
end