cart_pos  = out.cart_pos;
cart_ref  = out.cart_ref;
angle_q   = out.angle_q;
force_raw = out.force_raw;
force_out = out.force_out;

figure;

subplot(3,1,1);
plot(cart_pos.Time, cart_pos.Data, 'LineWidth',1.5); hold on;
plot(cart_ref.Time, cart_ref.Data, '--');
ylabel('Cart (m)'); legend('p','ref'); grid on;

subplot(3,1,2);
plot(angle_q.Time, rad2deg(angle_q.Data), 'LineWidth',1.5);
yline(rad2deg(0.1908),'--');
ylabel('Angle (deg)'); grid on;

subplot(3,1,3);
plot(force_raw.Time, force_raw.Data, '--'); hold on;
plot(force_out.Time, force_out.Data, 'LineWidth',1.5);
yline(20,':r'); yline(-20,':r');
ylabel('Force (N)'); xlabel('Time (s)');
legend('before sat','after sat'); grid on;

fprintf('Peak force: %.2f N\n', max(abs(force_out.Data)));
fprintf('Final cart position: %.4f m\n', cart_pos.Data(end));
fprintf('Final angle: %.4f rad\n', angle_q.Data(end));