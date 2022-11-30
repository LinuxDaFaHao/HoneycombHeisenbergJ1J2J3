hp = [0.01, 0.1,0.2];
Tc1 = [0.5, 0.5, 0.5];
Tc2 = [0.908, 0.908, 0.908];
plot(Tc1, hp,'-o');hold on;
plot(Tc2, hp,'-o');hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$T$','Interpreter','latex');
ylabel('$h_p$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
set(gca, 'Ylim',[0,inf]);
set(gca, 'Xlim',[0,inf]);