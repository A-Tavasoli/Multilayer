% Parameters
m = 1;   % mass
l = 0.7;  % length
g = 9.81;  % acceleration due to gravity

% Define the energy function
energy_function = @(theta, theta_dot) 0.5 * m * l^2 * theta_dot.^2 + m * g * l * (1 - cos(theta));

% Generate a grid of theta and theta_dot values
theta_range = linspace(-2*pi, 2*pi, 100);
theta_dot_range = linspace(-10, 10, 100);
[Theta, ThetaDot] = meshgrid(theta_range, theta_dot_range);

% Compute energy values for each combination of theta and theta_dot
EnergyValues = energy_function(Theta, ThetaDot);

% Plot the smooth energy surface with transparency
figure;
surf(Theta, ThetaDot, EnergyValues, 'EdgeColor', 'none', 'FaceAlpha', 0.7); % Remove edges and adjust transparency
hold on;

% Add level sets using contour lines
contour(Theta, ThetaDot, EnergyValues, 3, 'LineColor', 'k', 'LineWidth', 2);  % Adjust the number of contours and line properties

xlabel('\theta (radians)');
ylabel('d\theta/dt');
zlabel('Energy');
%title('Smooth Energy Surface with Level Sets for Simple Pendulum');
% grid off;

set(gca, 'XTick', [-pi,0,pi], 'YTick', [], 'ZTick', []);
ax.XAxis.Axle.LineStyle = 'none';
ax.YAxis.Axle.LineStyle = 'none';
% Display a colorbar
%colorbar;
