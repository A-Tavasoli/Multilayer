% Parameters
g = 9.81;  % acceleration due to gravity
l = 1.0;   % length of the pendulum
m = 0.1;   % mass at the end of the pendulum

% Initial conditions
theta0 = pi;    % initial angle (in radians)
omega0 = 0.1;       % initial angular velocity

% Time span for simulation
tspan = [0 20];    % simulation time span

% Define the system of differential equations
ode = @(t, y) [y(2); -(g/l)*sin(y(1))];

% Solve the differential equations
[t, y] = ode45(ode, tspan, [theta0; omega0]);

% Create a video file in MP4 format
videoFile = VideoWriter('pendulum_simulation_0.mp4', 'MPEG-4');
open(videoFile);

% Create a figure
figure;

% Loop through each time step and plot the pendulum
for i = 1:length(t)
    % Convert polar coordinates to Cartesian coordinates
    x = l*sin(y(i, 1));
    y_pendulum = -l*cos(y(i, 1));

    % Plot the pendulum line in blue
    plot([0, x], [0, y_pendulum], 'LineWidth', 7, 'Color', 'blue');
    hold on;

    % Plot the end mass three times larger in red
    mass_size = 30;  % increase the size of the mass
    plot(x, y_pendulum, 'o', 'MarkerSize', mass_size, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');

    % Plot the hinge in black
    plot(0, 0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black');

    axis equal;
    axis([-1.5*l, 1.5*l, -1.5*l, 1.5*l]);
    title(['Time: ' num2str(t(i)) ' seconds']);

    % Remove x- and y-ticks
    set(gca, 'XTick', [], 'YTick', []);

    % Capture the frame and add it to the video file
    frame = getframe(gcf);
    writeVideo(videoFile, frame);
    
    % Pause for a short duration to control the speed of the simulation
    pause(0.001);
    
    % Clear the figure for the next frame
    clf;
end

% Close the video file
close(videoFile);

disp('Simulation complete.');
