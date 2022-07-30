AFMhoneycomb = spinw;
AFMhoneycomb.genlattice('lat_const', [3, 3,10],'angled',[90 90 120],'spgr','P -3');
AFMhoneycomb.addatom('r',[2/3 1/3 0], 'S', 1, 'label', 'MCu1','color','r');
fprintf('Magnetic atom positions:\n');
AFMhoneycomb.table('matom')

AFMhoneycomb.gencoupling('maxDistance',4);

AFMhoneycomb.addmatrix('value',-5,'label','J1','color','SteelBlue')
AFMhoneycomb.addmatrix('value',-0.2,'label','J2','color','SteelBlue')
AFMhoneycomb.addmatrix('value',38,'label','J3','color','y')
AFMhoneycomb.addmatrix('value',diag([0 0 0.1]),'label','D','color','r')

AFMhoneycomb.addcoupling('mat','J1','bond',1)
AFMhoneycomb.addcoupling('mat','J2','bond',2)
AFMhoneycomb.addcoupling('mat','J3','bond',3)
AFMhoneycomb.addaniso('D')

AFMhoneycomb.genmagstr('mode','direct','S',cat(2,[0; 1; 0],[0;1;0]),'k',[1/2, 0,  0],'n', [0 0 1]);
% AFMhoneycomb.genmagstr('mode','direct','S',cat(2,[0; 0; 1],[0;0;1]),'k',[0, 0, 0],'n', [0 0 1]);
%   plot(AFMhoneycomb,'range',[6 6 1],'magColor','purple','baseShift',[0;-1;0],'atomLegend',false)
%plot(AFMhoneycomb,'range',[3 3 1],'atomColor','gold')
% swplot.zoom(2);

momentum_slice_num = 100;
momentum_slice_line = cell(1, 2 * momentum_slice_num + 1);
for i = 1:momentum_slice_num
    momentum_slice_line{ 2 * i-1} = [0, 1/momentum_slice_num*(i-1),0];
    momentum_slice_line{ 2 * i} = [1, 1/momentum_slice_num*(i-1),0];
end
momentum_slice_line{end} = momentum_slice_num;
honeySpec = AFMhoneycomb.spinwave(momentum_slice_line);
omega = honeySpec.omega(honeySpec.omega>0);
histogram(omega)

% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% ylabel('$\omega$(meV)','Interpreter','latex');
% xlabel('Momentum $k$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24);
% set(get(gca,'YLabel'),'FontSize',24);
