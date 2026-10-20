
%% INVERTED PENDULUM ON A CART - PHASE 1
% Lagrangian modelling, nonlinear equations, linearization,
% state-space representation, transfer functions and stability.
%
% Verify assumptions 1-4 against your SolidWorks assembly.
% All calculations use SI units.

clear; clc; close all;

%% 1. PARAMETERS

g = 9.81;                       % Gravity (m/s^2)

% Moving cart / clevis mass
M = 3.085144;                   % kg

% Pendulum rod
m = 0.578458;                   % kg
rodLength = 0.610;              % m
l = rodLength/2;                % assumed pivot-to-COM distance (m)

% SolidWorks inertia at centre of mass:
% conversion: 1 g.mm^2 = 1e-9 kg.m^2
Ixx_G = 2.007e7 * 1e-9;         % kg.m^2
Iyy_G = 23365.534 * 1e-9;       % kg.m^2
Izz_G = 2.006e7 * 1e-9;         % kg.m^2

% Assumed rotation axis: SolidWorks Z axis
I_G = Izz_G;

% Damping coefficients (unknown; initially idealized)
b = 0;                          % N.s/m
c = 0;                          % N.m.s/rad

%% 2. INERTIA ABOUT THE PIVOT

% Parallel-axis theorem
I_p = I_G + m*l^2;

% Determinant of the coupled acceleration matrix
Delta = (M+m)*I_p - (m*l)^2;

assert(Delta > 0, ...
    'Invalid parameters: determinant must be positive.');

fprintf('\n===== MODEL PARAMETERS =====\n');
fprintf('Cart mass M              = %.6f kg\n', M);
fprintf('Pendulum mass m          = %.6f kg\n', m);
fprintf('Rod length               = %.6f m\n', rodLength);
fprintf('Pivot-to-COM distance l  = %.6f m (assumed)\n', l);
fprintf('COM inertia I_G          = %.8f kg.m^2\n', I_G);
fprintf('Pivot inertia I_p        = %.8f kg.m^2\n', I_p);
fprintf('Delta                    = %.8f\n', Delta);

%% 3. LAGRANGIAN DERIVATION
%
% Generalized coordinates:
% q = [x; theta]
%
% x     = cart position
% theta = pendulum angle from the upright vertical
%
% Pendulum COM coordinates:
% X_G = x + l*sin(theta)
% Y_G = l*cos(theta)
%
% Squared COM speed:
% v_G^2 = xdot^2 + 2*l*xdot*thetadot*cos(theta)
%         + l^2*thetadot^2
%
% Kinetic energy:
% T = 0.5*(M+m)*xdot^2
%     + m*l*xdot*thetadot*cos(theta)
%     + 0.5*I_p*thetadot^2
%
% Potential energy:
% V = m*g*l*cos(theta)
%
% Lagrangian:
% L = T - V
%
% Lagrange equations:
% d/dt(dL/dxdot) - dL/dx = u - b*xdot
% d/dt(dL/dthetadot) - dL/dtheta = -c*thetadot
%
% Nonlinear equations:
% (M+m)*xdd + b*xdot + m*l*cos(theta)*thetadd
%   - m*l*sin(theta)*thetadot^2 = u
%
% m*l*cos(theta)*xdd + I_p*thetadd
%   - m*g*l*sin(theta) + c*thetadot = 0

%% 4. NONLINEAR EQUATIONS IN ACCELERATION FORM
%
% Define:
% r1 = u - b*xdot + m*l*sin(theta)*thetadot^2
% r2 = m*g*l*sin(theta) - c*thetadot
%
% Then:
% [M+m,     m*l*cos(theta)] [xdd    ] = [r1]
% [m*l*cos, I_p           ] [thetadd]   [r2]
%
% The nonlinearODE function at the end solves this system.

%% 5. LINEARIZATION ABOUT THETA = 0
%
% Small-angle approximations:
% sin(theta) ~ theta
% cos(theta) ~ 1
% theta_dot^2 is neglected as a second-order term.
%
% Linearized equations:
% (M+m)*xdd + b*xdot + m*l*thetadd = u
% m*l*xdd + I_p*thetadd - m*g*l*theta
%   + c*thetadot = 0
%
% State vector:
% z = [x; xdot; theta; thetadot]
%
% State-space form:
% zdot = A*z + B*u
% y = C*z + D*u
%
% Outputs: y = [x; theta]

A = [ ...
    0, 1, 0, 0;
    0, -I_p*b/Delta, -(m^2*g*l^2)/Delta, m*l*c/Delta;
    0, 0, 0, 1;
    0, m*l*b/Delta, (M+m)*m*g*l/Delta, -(M+m)*c/Delta
    ];

B = [ ...
    0;
    I_p/Delta;
    0;
    -m*l/Delta
    ];

C = [ ...
    1, 0, 0, 0;
    0, 0, 1, 0
    ];

D = [0; 0];

sys = ss(A,B,C,D);

sys.StateName = {'x','xdot','theta','thetadot'};
sys.InputName = {'u'};
sys.OutputName = {'cart_position','pendulum_angle'};

fprintf('\n===== STATE-SPACE MODEL =====\n');
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);

%% 6. TRANSFER FUNCTIONS

% Transfer function from cart force u to cart position x
G_cart = minreal(tf(sys(1,1)));

% Transfer function from cart force u to pendulum angle theta
G_theta = minreal(tf(sys(2,1)));

fprintf('\n===== TRANSFER FUNCTIONS =====\n');
disp('Cart position X(s)/U(s):');
G_cart

disp('Pendulum angle Theta(s)/U(s):');
G_theta

%% 7. OPEN-LOOP EIGENVALUES AND STABILITY

lambda = eig(A);

fprintf('\n===== OPEN-LOOP EIGENVALUES =====\n');
disp(lambda);

if any(real(lambda) > 1e-8)
    fprintf(['The upright equilibrium is unstable: ', ...
        'at least one eigenvalue has positive real part.\n']);
else
    fprintf(['No positive-real-part eigenvalue detected. ', ...
        'Check the angle convention and parameters.\n']);
end

%% 8. PLOT LINEAR MODEL RESPONSE
% Unit-force input. This is not a closed-loop controller.

figure('Name','Linear State-Space Response');
step(sys,5);
grid on;
title('Open-Loop Linear Response to a Unit Force');
legend('Cart position x (m)', ...
       'Pendulum angle theta (rad)', ...
       'Location','best');

%% 9. OPTIONAL NONLINEAR OPEN-LOOP SIMULATION

% State initial conditions:
% x = 0 m
% xdot = 0 m/s
% theta = 2 degrees
% thetadot = 0 rad/s

uConstant = 0;                  % N
tspan = [0 3];                  % seconds
z0 = [0; 0; deg2rad(2); 0];

odefun = @(t,z) nonlinearODE( ...
    t,z,M,m,l,I_p,g,b,c,uConstant);

[tNL,zNL] = ode45(odefun,tspan,z0);

figure('Name','Nonlinear Pendulum Angle');
plot(tNL,rad2deg(zNL(:,3)),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Pendulum angle (degrees)');
title('Nonlinear Open-Loop Pendulum Response');

figure('Name','Nonlinear Cart Position');
plot(tNL,zNL(:,1),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Cart position (m)');
title('Nonlinear Open-Loop Cart Motion');

%% 10. SAVE THE MODEL

save('inverted_pendulum_phase1.mat', ...
    'M','m','l','rodLength','g','I_G','I_p', ...
    'Ixx_G','Iyy_G','Izz_G','b','c', ...
    'A','B','C','D','sys','G_cart','G_theta','lambda');

fprintf('\nModel saved as inverted_pendulum_phase1.mat\n');

%% LOCAL FUNCTION: NONLINEAR EQUATIONS OF MOTION
function dz = nonlinearODE(~,z,M,m,l,I_p,g,b,c,u)

    xdot = z(2);
    theta = z(3);
    thetadot = z(4);

    % Right-hand-side terms
    r1 = u - b*xdot + m*l*sin(theta)*thetadot^2;
    r2 = m*g*l*sin(theta) - c*thetadot;

    % Determinant of the nonlinear mass matrix
    DeltaNL = (M+m)*I_p - (m*l*cos(theta))^2;

    % Solve for cart and pendulum angular acceleration
    xdd = (I_p*r1 - m*l*cos(theta)*r2)/DeltaNL;

    thetadd = ((M+m)*r2 - m*l*cos(theta)*r1)/DeltaNL;

    dz = [xdot; xdd; thetadot; thetadd];
end
