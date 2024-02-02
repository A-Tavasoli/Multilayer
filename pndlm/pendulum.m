% Parameters
m = 1;   % mass
l = 10;  % length

% Equilibrium points
a = 1;
equilibrium_points = [0 ] ; 

% Plot the phase plane trajectories
figure;
xlabel('\theta (radians)');
ylabel('d\theta/dt Velocity');
% title('Phase Plane Trajectories of Simple Pendulum');
grid on;
hold on;

for i = 1:length(equilibrium_points)
    % Define the equilibrium point
    equilibrium_point = equilibrium_points(i);
    
    % Initial conditions around the equilibrium point
    initial_conditions = [equilibrium_point; 0.3];
    
    % Solve the differential equations
    [t, angles] = ode45(@(t, y) [y(2); -9.81 / l * sin(y(1))], [0 7], initial_conditions);
    
    % Plot the trajectory
    plot(angles(:,1), angles(:,2), 'LineWidth', 1.5);
    
    % Display the initial condition point
    % plot(initial_conditions(1), initial_conditions(2), 'ro', 'MarkerSize', 10);
    % text(initial_conditions(1), initial_conditions(2), sprintf('  Equilibrium Point %.2f', equilibrium_point));
end

% Show a legend
legend('Trajectories', 'Location', 'Northeast');
xlim([-2*pi, 2*pi])

set(gca, 'XTick', [], 'YTick', []);
ax.XAxis.Axle.LineStyle = 'none';
ax.YAxis.Axle.LineStyle = 'none';
