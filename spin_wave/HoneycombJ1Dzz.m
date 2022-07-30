AFMhoneycomb = spinw;
AFMhoneycomb.genlattice('lat_const', [3, 3,10],'angled',[90 90 120],'spgr','P -3');
AFMhoneycomb.addatom('r',[2/3 1/3 0], 'S', 2, 'label', 'MCu1','color','r');
fprintf('Magnetic atom positions:\n');
AFMhoneycomb.table('matom')

AFMhoneycomb.gencoupling('maxDistance',4);

AFMhoneycomb.addmatrix('value',1,'label','J','color','SteelBlue')
AFMhoneycomb.addmatrix('value',diag([0 0 0]),'label','D','color','r')

AFMhoneycomb.addcoupling('mat','J','bond',1)
AFMhoneycomb.addaniso('D')
% AFMhoneycomb.symbolic(true);
  AFMhoneycomb.genmagstr('mode','direct','S',cat(2,[0; 1; 0],[0;-1;0]),'k',[0, 0,  0],'n', [0 0 1]);
% AFMhoneycomb.genmagstr('mode','direct','S',cat(2,[0; 0; 1],[0;0;-1]),'k',[0, 0, 0],'n', [0 0 1]);
%  plot(AFMhoneycomb,'range',[3 3 1],'magColor','purple','baseShift',[0;-1;0],'atomLegend',false)
%  plot(AFMhoneycomb,'range',[3 3 1],'atomColor','gold')
% swplot.zoom(2);

% symSpec = AFMhoneycomb.spinwavesym();
% pretty(symSpec.omega)
honeySpec = AFMhoneycomb.spinwave({[0 0 0] [1/3 1/3 0] [1/2 0 0] [0 0 0] 500});
% honeySpec = sw_neutron(honeySpec);
% 
% figure
% sw_plotspec(honeySpec,'mode','disp','axLim',[0 6],'colormap',[0 0 0],'colorbar',false)
% 
% honeySpec = sw_egrid(honeySpec,'Evect',linspace(0,7,500),'component','Sperp');
% 
% 
% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% ylabel('$\omega$(meV)','Interpreter','latex');
% xlabel('Momentum $k$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);

% figure
% sw_plotspec(honeySpec,'mode','color','axLim',[0 2],'dE',0.4)
% 
% honeySpec = sw_egrid(honeySpec,'Evect',linspace(0,7,500),'component',{'0.5*Sxx+0.5*Syy' 'Szz'});
% figure
% sw_plotspec(honeySpec,'mode','color','axLim',[0 1],'dE',0.4)
% sw_plotspec(honeySpec,'mode','disp','axLim',[0 7],'colormap',[0 0 0],'colorbar',false,'lineStyle','-','legend',false)