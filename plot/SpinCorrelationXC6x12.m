distance = 1:12;
sxsx = [-0.650453372402, 0.65093148338, -0.675021557106, 0.635726947422, ...
        -0.623123844322, 0.609090324233,-0.605759210996, 0.595396588469, ...
        -0.594014741393, 0.585477046366, -0.586500093821, 0.578102028535];
szsz = [-0.12044848442, 0.119771224347, -0.138227003202, 0.109661328578, ...
        -0.101065964308, 0.0907179199567, -0.0883253253962, 0.0812490659094, ...
        -0.0795202851681, 0.0745501872496, -0.0735820609636, 0.0698473462981];
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
