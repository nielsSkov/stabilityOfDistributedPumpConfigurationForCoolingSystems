clear all, close all; clc; clear classes                                   %#ok<CLCLS>
                                                                           %#ok<*UNRCH>
%% Load Model Parameters

run('modelParameters')

%% Plot Options

plotPressure  = 1;
plotFlow      = 1;
plotPumpSpeed = 1;

%% Prep for Python (used in simSys_nEq4)

%add current folder to python search path
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

%option to turn python module reload on/off
reloadPyMod = 1;

%reload module (this will allow script changes to take effect)
if reloadPyMod
	mod = py.importlib.import_module('pyMinimize');
	py.importlib.reload(mod);
end

%% Simulation Setup

theta_c = 10 + 273.15; %[K]
T_a     = 30 + 273.15; %[K]

%sample time
Ts = .2; %[s]

%length of simulation
T_final = 15*60; %[s]

%initial states
theta_0 = T_a;% 30 + 273.15; %[K]
T_0     = T_a;% 30 + 273.15; %[K]
zeta_0  = 0;

%initialization for ode15s
tspan = 0:Ts:T_final;
init = [ theta_0  T_0 zeta_0 ...
         theta_0  T_0 zeta_0 ...
         theta_0  T_0 zeta_0 ...
         theta_0  T_0 zeta_0 ];

%option to lower relative tollerence (default 1e-3)
options = odeset('RelTol',1e-3);

%pump speeds
w_init = [ 0 0 0 0 ]';

%switch control on/off
conOn = 1;

%temperature refference
T_eq = ones(1,4)*(20 + 273.15);

%control gain vector
K = -[ -0.07538319 -0.03915181 -0.07705162 ];

%% Simulation

[ t, x ] = ode15s( @(t,x) simSys_nEq4( t, x, w_init, T_eq, r, R, R_c,...
                                       a, b, theta_c, T_a, K, conOn, ...
                                       B, C_w, C_a, V_w, V_a, Q, pth ), ...
                   tspan, init, options                                 );

%initializing for simulation re-run
q     = zeros(length(t),4);
p     = zeros(length(t),4);
w     = zeros(length(t),4);
w_err = zeros(length(t),4);

%re-running sim-function in loop to extract parameters at each time step
for i = 1:length(t)

  [ ~,  q(i,:), ...
        p(i,:), ...
        w(i,:), ...
	  w_err(i,:) ] = simSys_nEq4( t(i), x(i,:)', w_init, T_eq, r, R, R_c, ...
                                a, b, theta_c, T_a, K, conOn,           ...
                                B, C_w, C_a, V_w, V_a, Q, pth           );
end

%% Plot Simulation Results

%plot output temperatures (return water and exhaust air)
figure
tiledlayout(2,2); hAx = nexttile;
for i = 1:4
	if i > 1, nexttile; end
	plot(t,x(:,3*i-2)-273.15)
	hold on
	plot(t,x(:,3*i-1)-273.15)
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Temperature [\circC]')
	ylim([ 0 40 ])
end
legend( hAx, 'Return Water Temperature, \theta', ...
             'Exhaust Air Temperature, T',       ...
             'Location','NorthOutside'           );

figure
plot(t,x(:,3)), hold on, plot(t,x(:,6)), plot(t,x(:,9)), plot(t,x(:,12))
grid on, grid minor
xlabel('time [s]')
ylabel('Integral State')
legend( 'integral(T1-T*)', 'integral(T2-T*)', ...
        'integral(T3-T*)', 'integral(T4-T*)', 'Location', 'northeast')

if plotFlow
	figure
	plot(t,q(:,1)), hold on, plot(t,q(:,2)), plot(t,q(:,3)), plot(t,q(:,4))
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Flow [m^3/h]')
	legend('q1', 'q2','q3', 'q4', 'Location', 'northwest')
end

if plotPressure
	figure
	plot(t,p(:,1)), hold on, plot(t,p(:,2)), plot(t,p(:,3)), plot(t,p(:,4))
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Pump Pressure [bar]')
	legend('p1', 'p2','p3', 'p4', 'Location', 'northwest')
end

if plotPumpSpeed
	figure
	plot(t,w(:,1)), hold on, plot(t,w(:,2))
	plot(t,w(:,3)),          plot(t,w(:,4))
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Pump Speed, \omega')
	legend( '\omega_1', '\omega_2','\omega_3', '\omega_4', ...
	        'Location', 'northwest'                        )
end

if plotPumpSpeed
	figure
	plot(t,w_err(:,1)), hold on, plot(t,w_err(:,2))
	plot(t,w_err(:,3)),          plot(t,w_err(:,4))
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Error in Pump Speed from Optimization, \omega')
	legend( '\omega_1_{,err}', '\omega_2_{,err}',   ...
	        '\omega_3_{,err}', '\omega_4_{,err}',   ...
	        'Location', 'northwest'               )
end
