function [x, y] = HoneycombYCCylinderSiteInd2XYCoor(site_idx, ly)
% for the definition, refer to
% 'https://online.kitp.ucsb.edu/online/fragnets12/honeydmrg/pdf1/Zhu_Honeycomb_Fragnets12_KITP.pdf'
% page 5
% site_idx count from 0
% assume |a1| = |a2| = 1, a1 and a2 is the unit cell vector

y_idx = mod(site_idx, ly);
x_idx = floor(site_idx / ly);

x = x_idx / 2;
a_sublattice = mod(x_idx + y_idx, 2);
y = y_idx * sqrt(3) / 2 + (1 - a_sublattice) / (2 * sqrt(3));
end