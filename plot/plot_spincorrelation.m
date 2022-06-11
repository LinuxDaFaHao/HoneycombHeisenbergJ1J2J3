
geometry  = 'XC';
Lx = 12;
Ly = 6;
J1 = -1;
J2 = -0.04;
J3 = 5.3;
Dzz = 0.02;
Db = 10000; %bond dimension
postfix = ['HoneyHei',geometry, num2str(Ly),'x', num2str(Lx), 'J1',num2str(J1), 'J2', num2str(J2), 'J3', num2str(J3), 'Dzz', num2str(Dzz), 'D',num2str(Db),'.json'];
prefix = '../data/';
sz_data = jsondecode(fileread([prefix, 'sz',postfix]));
szsz_data = jsondecode(fileread([prefix, 'zzsf',postfix]));
spsm_data = jsondecode(fileread([prefix, 'pmsf',postfix]));
smsp_data = jsondecode(fileread([prefix, 'mpsf',postfix]));

start_point = szsz_data{1}{1}(1);
for i = 1:numel(szsz_data)
    if(szsz_data{i}{1}(1) ~= start_point)
        correlation_data_size = i;
        break
    end
end

szsz_data = szsz_data(1:correlation_data_size);
