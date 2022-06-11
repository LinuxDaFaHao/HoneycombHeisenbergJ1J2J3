distance = 1:8;
sxsx = [0.121731583604, -0.278289110388, 0.10326408674, 0.525893138739, 0.0915916708728, -0.236979496978, 0.0761102644151, 0.366409994962];
szsz = [0.0607101844996, -0.132275875106, 0.0512804782053, 0.255269081079, 0.0453833443074, -0.111603088736, 0.037575718528, 0.175627404103 ];
h1 = plot(distance, abs(sxsx)/2,'-x');hold on;
h2 = plot(distance, abs(szsz),'-o');hold on;


l=legend([h1,h2],'$|\langle S^x(x_0) S^x(x_0+r)\rangle| = |\langle S^y(x_0) S^y(x_0+r)\rangle|$','$|\langle S^z(x_0) S^z(x_0 + r)\rangle |$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'Children'),'markersize',10); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('spin correlations','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
