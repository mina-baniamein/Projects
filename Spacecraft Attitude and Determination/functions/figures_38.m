%% PLOT

figure
plot(out.time, out.Am_1 , 'LineWidth', 2)
grid on
hold on
plot(out.time, out.Am_2 , 'LineWidth', 2)
plot(out.time, out.Am_3 , 'LineWidth', 2)
legend('A11', 'A12', 'A13','A21', 'A22', 'A23','A31', 'A32', 'A33')
xlabel('time')
ylabel('A_{B/N}^m')
title('Measured Attitude')

figure
plot(out.time, out.Ad_1 , 'LineWidth', 2)
grid on
hold on
plot(out.time, out.Ad_2 , 'LineWidth', 2)
plot(out.time, out.Ad_3 , 'LineWidth', 2)
legend('A11', 'A12', 'A13','A21', 'A22', 'A23','A31', 'A32', 'A33')
xlabel('time')
ylabel('A_{B/N}')
title('Desired Attitude')
ylim([-1.5 1.5])


figure
plot(out.time, out.Alpha_err , 'LineWidth', 1)
grid on
hold on
xlabel('time [s]')
ylabel('Alpha_{err} [°]')
title('Pointing error')
ylim([-1.5 4.5])

figure
plot(out.time, out.w_x , 'LineWidth', 1)
grid on
hold on
plot(out.time, out.w_y , 'LineWidth', 1)
plot(out.time, out.w_z , 'LineWidth', 1)
legend('w_x', 'w_y', 'w_z')
xlabel('time [s]')
ylabel('Angular velocity [rad/s]')
title('Angular velocity')


figure
plot(out.time, out.Ae_1 , 'LineWidth', 1)
grid on
hold on
plot(out.time, out.Ae_2 , 'LineWidth', 1)
plot(out.time, out.Ae_3 , 'LineWidth', 1)
legend('A11', 'A12', 'A13','A21', 'A22', 'A23','A31', 'A32', 'A33')
xlabel('time [s]')
ylabel('A_{err} [-]')
title('Error matrix')


%% 2Bp data
% Plot the results
Earth_plot
hold on
plot3( Y(:,1), Y(:,2), Y(:,3), '-' , 'LineWidth', 2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal;
grid on;

plot3(r_0(1),r_0(2),r_0(3), 'o')