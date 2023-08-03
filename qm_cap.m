% Constants
e = 1.602e-19; % electron charge
k_B = 1.38064852e-23; % Boltzmann constant
T = 300; % temperature in Kelvin
mass_mos2 = 2.12488e-21;
mass_nc = 3.18161e-21;

% Density of states (DOS) calculated from DFT


% Define the thermal broadening function
thermal_broadening = @(E, Vex, T) sech((e*(E-Vex))/(2*k_B*T)).^2;

% Define the range of external potentials
Vext_range = -0.5:0.01:0.5;

% Initialize an array to store the capacitance values
Cq = zeros(size(Vext_range));

% Sweep the external potential
for i = 1:length(Vext_range)
    Vex = Vext_range(i);
    %Calculate the product of DOS and thermal broadening function
    dosF_T = nc_o_dos(:, 2) .* thermal_broadening(nc_o_dos(:, 1),Vex, T);
    % Define the integral function
    integrand = @(E) interp1(nc_o_dos(:, 1), dosF_T, E, 'spline');
    % Calculate the definite integral using integral
    Cq(i) = e^2*(4*k_B*T)^-1*integral(integrand, -inf, inf);
end

% Plot the results

plot(Vext_range, Cq/mass_nc);
xlabel('External Potential (V)');
ylabel('Quantum Capacitance (F/g)');
