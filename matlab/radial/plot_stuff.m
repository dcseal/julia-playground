%PLOT_STUFF - plotting routine(s)

figure(1);
clf
p = pcolor(XX',YY',q');
set(p,'EdgeColor','none');
%contourf(XX',YY',q');
t1 = title(['Advection equation curvilinear example t = ' num2str(0)]);
set(t1,'FontSize',12);
xlabel('x', 'FontSize', 16); ylabel('y', 'FontSize', 16);

figure(2);
clf
p = pcolor(XI',ET',qt');
set(p,'EdgeColor','none');
%contourf(XI',ET',qt');
t1 = title(['Advection equation curvilinear example - computational domain t = ' num2str(0)]);
set(t1,'FontSize',12);
xlabel('\xi','FontSize', 16); ylabel('\eta', 'FontSize', 16);
pause(0.1);

% TODO - save data to output files