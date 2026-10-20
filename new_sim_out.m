simout = out.simout;
t = simout.Time;  x = simout.Data(:,1);  th = simout.Data(:,2);  u = simout.Data(:,3);

blue = [0.00 0.45 0.74];  red = [0.85 0.33 0.10];  grey = [0.4 0.4 0.4];

fig = figure('Color','w','Position',[100 100 760 700]);
tl  = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
plot(t,x,'LineWidth',1.6,'Color',blue); hold on
yline(0.5,'--','x_{ref} = 0.5 m','Color',grey,'LabelHorizontalAlignment','left');
xline(1,':','Color',grey);
ylabel('x (m)'); title('Cart position'); grid on

ax2 = nexttile;
plot(t,th,'LineWidth',1.6,'Color',red); hold on
yline(0,'-','Color',grey);
xline(1,':','Color',grey);
ylabel('\theta (rad)'); title('Pendulum angle'); grid on

ax3 = nexttile;
m = t > 0.01;                              % hide the start-up spike
plot(t(m),u(m),'LineWidth',1.6,'Color',blue); hold on
yline(0,'-','Color',grey);
xline(1,':','Color',grey);
ylabel('u (N)'); title('Control force (start-up spike omitted)'); grid on
xlabel('Time (s)');

linkaxes([ax1 ax2 ax3],'x');  xlim([0 12]);
set(findall(fig,'Type','axes'),'FontSize',11,'Box','on');
exportgraphics(fig,'phase2_closed_loop.png','Resolution',300);