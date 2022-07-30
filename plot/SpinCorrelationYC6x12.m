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
for i = 1:datanum_of_structfactor
    szsz_correlation = SzSzCorrelationData{i}{2};
    spsm_correlation = SpSmCorrelationData{i}{2};
    smsp_correlation = SmSpCorrelationData{i}{2};
    
    sdots_data(i) = szsz_correlation + 1/2 * ( spsm_correlation + smsp_correlation );
end
moment = sqrt(mean(abs(sdots_data)));
fprintf("the magnetic moment = %.12f\n", moment);

reference_site = SzSzCorrelationData{1}{1}(1);
i = 1;
num_of_chaincorrelation = 0;
while(SzSzCorrelationData{i}{1}(1) == reference_site)
    site2 = SzSzCorrelationData{i}{1}(2);
    if(mod(site2 - reference_site, Ly ) == 0 )
        num_of_chaincorrelation = num_of_chaincorrelation + 1;
    end
    i = i + 1;
end
szsz_correlation = zeros(1, num_of_chaincorrelation);
spsm_correlation = zeros(1, num_of_chaincorrelation);
smsp_correlation = zeros(1, num_of_chaincorrelation);
i = 1;
num_of_chaincorrelation = 0;
sz_reference_point = SzData(reference_site, 2);
while(SzSzCorrelationData{i}{1}(1) == reference_site)
    site2 = SzSzCorrelationData{i}{1}(2);
    if(mod(site2 - reference_site, Ly ) == 0 )
        num_of_chaincorrelation = num_of_chaincorrelation + 1;
        sz_site2 = SzData(site2, 2);
        szsz_correlation(num_of_chaincorrelation) = SzSzCorrelationData{i}{2} - sz_reference_point * sz_site2;
        spsm_correlation(num_of_chaincorrelation) = SpSmCorrelationData{i}{2};
        smsp_correlation(num_of_chaincorrelation) = SmSpCorrelationData{i}{2};
    end
    i = i + 1;
end
distance = 1:num_of_chaincorrelation;
plot(distance, szsz_correlation,'x');hold on;
plot(distance, spsm_correlation,'o');hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('distance $x$','Interpreter','latex');
ylabel('correlator $S^{\alpha}(0) S^{\alpha}(x)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
% set(gca, 'Ylim',[0,inf]);