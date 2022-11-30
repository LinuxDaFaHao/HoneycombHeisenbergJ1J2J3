L = 32;
beta = 0.2041;
J1zz = 5.31;
J2zz = 0.2;
J3zz = -28;
Dzz = -0.113;
num_chain = 12;
prefix = '../data/';
xy_real = [];
xy_imag = [];
for i = 0:num_chain-1
    postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
        'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
        'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
    xy_real_data = load([prefix, 'xy_order_parameter_real', postfix]);
    xy_imag_data = load([prefix, 'xy_order_parameter_imag', postfix]);
    xy_real = [xy_real; xy_real_data(500001:end)];
    xy_imag = [xy_imag; xy_imag_data(500001:end)];
end
hist3([xy_real, xy_imag],[100,100],'CdataMode','auto');
colorbar
view(2)

% axis equal;