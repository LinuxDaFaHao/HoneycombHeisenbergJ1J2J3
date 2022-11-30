beta = 0.1538;
L = 32;
J1zz = 5.343;
J2zz = 0.2;
J3zz = -28.166;
Dzz = -0.113;
num_chain = 12;
prefix = '../data/';
for i = 0:num_chain-1
    postfix = ['hei-rank',num2str(i), 'Honeycomb3', 'J1zz', num2str(J1zz,'%.6f'),...
        'J2zz', num2str(J2zz,'%.6f'),  'J3zz', num2str(J3zz,'%.6f'),...
        'Dzz', num2str(Dzz,'%.6f'),'beta',num2str(beta,'%.6f'),'L', num2str(L)];
    Sx_data = load([prefix, 'magnetization0', postfix]);
    Sy_data = load([prefix, 'magnetization1', postfix]);
    Sz_data = load([prefix, 'magnetization2', postfix]);
    S_data = [Sx_data(end/2+1:end),Sy_data(end/2+1:end),Sz_data(end/2+1:end)]; 
end
