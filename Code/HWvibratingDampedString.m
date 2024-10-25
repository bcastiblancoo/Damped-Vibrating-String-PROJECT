clear all; %initialize the program clearing and closing all previous variables
close all;

%1. Define input data and discretization parameters of time

L=0.328; %Lenght of the string
rho=0.660E-3; %Density of the string
F=55.0; %Tension over the string
c=sqrt(F/rho); %Wave speed


gamma=400; % being gamma proportional to damping factor (multiply by c^2/2)

tmax=log(100)/gamma;

nst=13;
t = linspace(0, tmax, nst); % Discretize time over time in which atenues a factor of 10^-2


%2. Define array in space
npt=200;
step=L/npt;
x=linspace(0,L,npt);

%3. Computing initial conditions

%Let us define the two possible pinches

f1 = zeros(1, npt); % Initial displacement two pinches
f2 = zeros(1, npt); % Initial displacement two pinches

%Define parameters for the initial conditions
z1=0.15;
p1=L/4;

z2=-0.15;
p2=(3/4)*L;

%Initial deformations: two types

%One pinch f1(x)

for i=1:npt
    if (x(i)<=p1)
        f1(i)=z1*x(i)/p1;
    else 
        f1(i)=z1*(x(i)-L)/(p1-L);
    end    
end

%Two pinches f2(x)

for i=1:npt
    if (x(i)<=p1)
        f2(i)=z1*x(i)/p1;
    else 
        if (x(i)<=p2)
            f2(i)=((z1-z2)/(p1-p2))*(x(i)-p1)+z1;
        else 
            f2(i)=(z2*(x(i)-L))/(p2-L);
        end
    end    
end

%And the initial distribution of the velocity field g(x):

for i=1:npt
    g(i)=0;
end

%Plotting initial deformations:

figure(1);
plot(x,f1);
title('Initial condition: One pinch');
xlabel('x[m]');
ylabel('f1(x)[m]');
xlim([0,L]);
grid on;
saveas(gcf, 'OnePinchCondition.png');


figure(2);
plot(x,f2);
title('Initial condition: Two pinches');
xlabel('x[m]');
ylabel('f2(x)[m]');
xlim([0,L]);
grid on;
saveas(gcf, 'TwoPinchesCondition.png');

%4. Compute normal modes

%Now let us compute the normal modes
nmax=30;

for n=1:nmax
    for i=1:npt
        phi(n,i)=sqrt(2/L)*sin(n*pi*x(i)/L);
    end
end

%Let's plot the first seven
for n=1:7
 figure(3);
 plot(x/L,phi(n,:)); hold on;
 legendEntries{n} = sprintf('n = %d', n); % Store legend string for each plot
end

 xlim([0,1]); ylim([-sqrt(2/L),sqrt(2/L)]);
 legend(legendEntries, 'Location', 'northeastoutside');
 title('Normal eigenmodes');
 xlabel('x/L');
 ylabel('\phi(x)[m]');
 grid on;
 saveas(gcf, 'Eigenmodes.png');

%5. Display values in a table respecting significant digit rule

%Let us print some values

% print a table

for n=1:nmax
    k(n) = n*pi/L;
    omega(n)=n*pi*c/L;
    Omega(n) = sqrt((n*pi*c/L)^2-gamma^2);
    nu(n) = Omega(n)/(2*pi);  
    period(n) = 1/nu(n);  
end 
header{1}='n'; header{2}='k_{n}'; header{3}='\omega_{n}'; header{4}='\Omega_{n}'; header{5}='\nu_{n}'; header{6}='T_{n}=\nu_{n}^{-1}';   
printtable([(1:nmax)' k' omega' Omega' nu' period'], 'ColName', header, 'Precision',[1 3 3 3 3 3]);

%Let us make a graph to compare normal modes frequencies and damped ones
%with respect to \gamma


% Define the range of gamma
gammav = linspace(10, omega(1)/2, 100); % gamma from 10 to omega_1/2 with 100 points

%Let us consider the first eigenmode n=1
omega_1 = pi * c / L; % Constant omega_n for each n=1
Omega_1 = sqrt((omega_1)^2 - gammav.^2); % Omega_1 as a function of gamma

% Plot omega_1 vs gamma (dotted line)
figure;
plot(gammav, omega_1 * ones(size(gammav)), 'b:', 'LineWidth', 2); % Dotted line for omega_n
hold on;
% Plot Omega_n vs gamma (solid line)
plot(gammav, Omega_1, 'r-', 'LineWidth', 1); % Solid line for Omega_n
xlabel('\gamma [Hz]', 'FontSize', 12);   % Label for the x-axis
ylabel('\omega_1 and \Omega_1 [Hz]', 'FontSize', 12); % Label for the y-axis
title('Plot of \omega_1 and \Omega_1 vs. \gamma', 'FontSize', 14);
legend('\omega_1', '\Omega_1', 'Location', 'northeast');
grid on;
saveas(gcf, 'OmegasVsgamma.png');


%6. Test numerically the orthonormalisation with a table

%Let us verify orthonormalization


% Initialize the overlap matrix
overlapCol= zeros(n, 1);

% Compute the overlap integral for each pair of eigenfunctions
for i = 1:nmax
        % Compute the overlap integral using the trapezoidal rule
        overlapCol(i, :) = trapz(x, phi(i, :) .* phi(i, :));
end

% Display the overlap matrix
disp('Overlap:');
disp(overlapCol);


for i=1:nmax
    n(i) = i;      
    overlap(i) = overlapCol(i);  
end 
header{1}='n'; header{2}='<\phi_{n},\phi_{n}>';   
printtable([(1:nmax)' overlap'], 'ColName', header, 'Precision',[1 3]);

%7. Overlap with the initial conditions

% Initialize arrays to store overlap results
overlap_f = zeros(nmax, 1); % Overlap of f1 with each eigenfunction
overlap_g = zeros(nmax, 1); % Overlap of f2 with each eigenfunction

% Loop over each eigenfunction (mode number n)
for n = 1:nmax    
    % Compute the overlap of f1(x) with phi_n(x) using trapz
    overlap_f(n) = trapz(x, f2 .* phi(n,:));
    
    % Compute the overlap of f2(x) with phi_n(x) using trapz
    overlap_g(n) = trapz(x, g .* phi(n,:));
end

%Displaying the table

for i=1:nmax
    n(i) = i;      
    overlapwith_f(i) = overlap_f(i);
    overlapwith_g(i) = overlap_g(i);
end 
header{1}='n'; header{2}='<\phi_{n},f>'; header{3}='<\phi_{n},g>';   
printtable([(1:nmax)' overlapwith_f' overlapwith_g'], 'ColName', header, 'Precision',[1 3 3]);

% 8. Truncate series and plot

% Fourier series reconstruction for f2 and g
f_fourier = zeros(1, npt); % Initialize Fourier series for f2
g_fourier = zeros(1, npt); % Initialize Fourier series for g

for n = 1:nmax
    f_fourier = f_fourier + overlap_f(n) * phi(n, :); % Sum Fourier terms for f2
    g_fourier = g_fourier + overlap_g(n) * phi(n, :); % Sum Fourier terms for g
end

% Plotting original functions and their Fourier series expansions:
figure;
subplot(2, 1, 1);
plot(x, f2, 'b', 'LineWidth', 1); hold on;
plot(x, f_fourier, 'r--', 'LineWidth', 1);
title('Original two pinches f(x) and its Fourier series approx.');
legend('f(x)', 'Fourier series approx');
xlabel('x[m]');
ylabel('f(x)[m]');
xlim ([0 L]);
grid on;

subplot(2, 1, 2);
plot(x, g, 'b', 'LineWidth', 1); hold on;
plot(x, g_fourier, 'r--', 'LineWidth', 1);
title('Initial distribution of the velocity field g(x) and its Fourier series approx.');
legend('g(x)', 'Fourier Series approx');
xlabel('x[m]');
ylabel('g(x)[m/s]');
xlim ([0 L]);
grid on;

saveas(gcf, 'FourierApproximation.png');

% Compute the R^2 coefficient for f2
f2_mean = mean(f2); % Mean of the original f2
SS_res_f2 = sum((f2 - f_fourier).^2); % Residual sum of squares for f2
SS_tot_f = sum((f2 - f2_mean).^2); % Total sum of squares for f2
R2_f2 = 1 - SS_res_f2 / SS_tot_f; % R^2 for f2

% Compute the R^2 coefficient for g
g_mean = mean(g); % Mean of the original g
SS_res_g = sum((g - g_fourier).^2); % Residual sum of squares for g
SS_tot_g = sum((g - g_mean).^2); % Total sum of squares for g
R2_g = 1 - SS_res_g / SS_tot_g; % R^2 for g

% Display the R^2 values
fprintf('R^2 for f(x) Fourier series approximation: %.4f\n', R2_f2);
fprintf('R^2 for g(x) Fourier series approximation: %.4f\n', R2_g);

% Evaluate goodness of fit using the gof function for f and g
gof_f = gof(f_fourier, f2);
gof_g = gof(g_fourier, g);


% Display the results

disp('Goodness of fit for f(x) and its Fourier series expansion:');
disp(gof_f);

disp('Goodness of fit for g(x) and its Fourier series expansion:');
disp(gof_g);

%9. ploting succesive paterns of the amplitude

% Time evolution of the wavefunction Psi(x, t) for f1 pinch
Psi = zeros(nst, npt); % Initialize wavefunction Psi(x, t)

for n = 1:nmax
    omega_n = (n * pi * c) / L; % Angular frequency for mode n
    for j = 1:nst
        Psi(j, :) = Psi(j, :) + phi(n,:)*((overlap_f(n)) .* cos(omega_n*t(j)))*exp(-gamma*t(j));
    end
end

% Plot successive amplitude patterns for the wavefunction Psi(x, t)
figure;
hold on;

% Define a color map for better distinction between patterns
colors = jet(nst);

% Plot with labels for each pattern corresponding to a time step
for j = 1:nst
    plot(x, Psi(j, :), 'Color', colors(j, :), 'LineWidth', 1.5); % Plot the spatial wave pattern for each time step
    % Add label for each time step
    labels{j} = sprintf('t = %.3f s', t(j));
end

title('Successive amplitude patterns for \Psi(x,t) over one period');
xlabel('x[m]');
ylabel('Displacement[m]');
grid on;

% Add a legend with labels for the successive patterns
legend(labels, 'Location', 'northeastoutside');
hold off;
hold off;
saveas(gcf,'SuccesiveAmplitudePatterns.png');

% 10. Making a movie

% Set video parameters
fps = 60; % Frames per second for video
duration = tmax; % Total time duration (matching original request)
nframes = 283; % 100 frames for smooth video
time_values = linspace(0, duration, nframes); % Discretized time

% Prepare video writer
v = VideoWriter('dampedvibratingstring_two_pinches.mp4', 'MPEG-4');
v.FrameRate = fps; % Set frames per second
open(v); % Open the video file for writing

% Initialize figure for plotting
figure;

% Loop through 100 frames and plot the wavefunction evolution
for j = 1:nframes
    % Compute wavefunction for the current frame (cosine and sine evolution)
    Psi = zeros(1, npt);
    for n = 1:nmax
        omega_n = (n * pi * c) / L; % Angular frequency for mode n
        % Full wavefunction evolution with A_n (cosine term) and B_n (sine term)
         Psi = Psi + phi(n, :)*((overlap_f(n)).*cos(omega_n*time_values(j)))*exp(-gamma*time_values(j));
    end

    % Plot the wavefunction at the current time step
    plot(x, Psi, 'b', 'LineWidth', 2);
    title(sprintf('String Vibration at t = %.3f s, \\gamma = %.2f', time_values(j), gamma));
    xlabel('x[m]');
    ylabel('Displacement[m]');
    grid on;
    axis([0 L -0.4 0.4]); % Fixed axis limits for better visualization
    drawnow;

    % Capture the frame
    frame = getframe(gcf);

    % Write the frame to the video
    writeVideo(v, frame);
end

% Close the video writer
close(v); 

disp('The wavefunction evolution video with two pinches has been saved as dampedvibratingstring_two_pinches.mp4');