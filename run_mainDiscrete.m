clear all, close all; clc; clear classes                                   %#ok<CLCLS>
                                                                           %#ok<*UNRCH>
%% Load Model Parameters

run('modelParameters')

%% Prep for Python

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

%step size
h = .2; %[s]

%length of simulation
T_final = 8*60; %15*60; %[s]

%initial states
theta_0 = T_a;% 30 + 273.15; %[K]
T_0     = T_a;% 30 + 273.15; %[K]
zeta_0  = 0;

%time vector
t = 0:h:T_final;

%initialization
init = [ theta_0  T_0  ; % for subsystem n == 1
         theta_0  T_0  ; % for subsystem n == 2
         theta_0  T_0  ; % for subsystem n == 3
         theta_0  T_0 ]; % for subsystem n == 4

%switch control on/off
conOn = 1;

%temperature refference
T_eq = ones(1,4)*(20 + 273.15);

%control gain vector
K = -[ -0.07538319 -0.03915181 -0.07705162 ];

%initialize flows
qInit = q;

%init pump speeds - used as const pump speeds if control is off (con == 0)
w_init = [ .1 .1 .1 .1 ]';

%initialize vectors for discrete sim
x          = zeros(2,4,length(t));
q          = zeros(length(t),4);
p          = zeros(length(t),4);
w          = zeros(length(t),4);
err        = zeros(length(t),1);
T_integral = zeros(length(t),4);

%% First Step of Discrete Sim

%convert flow and pump speed to a format accepted by python
pth.q = py.list( qInit );
pth.w = py.list( w_init );

%solve using the Broyden–Fletcher–Goldfarb–Shanno algorithm
min_out = py.pyMinimize.BFGS( pth.q, pth.Q, pth.w, pth.R, ...
															pth.r, pth.a, pth.b, R_c    );

q(1,:) = [ min_out{1} min_out{2} min_out{3} min_out{4} ]';
err(1) = min_out{5};

for n = 1:4 %looping through each subsystem
	x(:,n,1) = init(n,:)';
end

%% Discrete Simulation - Euler's Method
for i = 2:length(t)
	
	if conOn
	
		%initialize integral anti windup
		antiWind = zeros(1,4);
		
		for n = 1:4 %looping through each subsystem
			
			if i > 2
				T_errThen = x(2,n,i-2) - T_eq(n);
			else
				T_errThen = 0;
			end
			
			T_errNow  = x(2,n,i-1) - T_eq(n);
			
			T_integral(i,n) = T_integral(i-1,n) + h*( T_errNow + T_errThen )/2;
			
			%calculate flow in equilibrium (T_eq is the desired temperature)
			theta_eq = ( -C_a*Q(n)*( T_a - T_eq(n) ) + B(n)*T_eq(n) )/B(n);

			%calculate return water temperature nn equilibrium
			q_eq = B(n)*(theta_eq - T_eq(n))/( C_w*(theta_c - theta_eq) );
			
			%states for calculating controlled punp speed
			y = [ x(1,n,i-1) - theta_eq  ;
						x(2,n,i-1) - T_eq(n)   ;
					  T_integral(i,n)       ];
			
			%pump speed (controlled)
			w(i,n) = K*y + q_eq;
			
			%set minimum allowed pump speed
			w_min  = 0.01;
			
			%option for hard upper limmit of pump speed
			w_max = inf;
			
			%limit pump speed
			if     w(i,n) < w_min, w(i,n) = w_min;
			elseif w(i,n) > w_max, w(i,n) = w_max;
			end
			
		end
	else
		w(i,:) = w_init;
	end
	
	%convert pump speed to a format accepted by python
	pth.w = py.list( w(i,:) );
	
	%initial guess for minimization
	pth.q = py.list( q(i-1,:) ); %use last flows
	
	%solve using the Broyden–Fletcher–Goldfarb–Shanno algorithm
	min_out = py.pyMinimize.BFGS( pth.q, pth.Q, pth.w, pth.R, ...
																pth.r, pth.a, pth.b, R_c    );
	
	q(i,:) = [ min_out{1} min_out{2} min_out{3} min_out{4} ]';
	err(i) = min_out{5};

	for n = 1:4 %looping through each subsystem
		A = [ ...
		-(q(i,n)/V_w(n) +B(n)/(C_w*V_w(n)))  B(n)/(C_w*V_w(n))                ;
		  B(n)/(C_a*V_a(n))                -(Q(n)/V_a(n) +B(n)/(C_a*V_a(n))) ];
		
		f = A*x(:,n,i-1) + [ (q(i,n)/V_w(n))*theta_c  ;
		                     (Q(n)/V_a(n))*T_a       ];
		
		x(1:2,n,i) = x(1:2,n,i-1) + h*f;
	end
end

%% Plot results

%plot output temperatures (return water and exhaust air)
figure
tiledlayout(2,2); hAx = nexttile;
for n = 1:4 %looping through each subsystem
	if n > 1, nexttile; end
	
	plot( t, squeeze( x(1,n,:) ) - 273.15 ), hold on
	plot( t, squeeze( x(2,n,:) ) - 273.15 )
	
	set(gca, 'XLimSpec', 'Tight');
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Temperature [\circC]')
	%ylim([ 0 40 ])
end
legend( hAx, 'Return Water Temperature, \theta', ...
             'Exhaust Air Temperature, T',       ...
             'Location','NorthOutside'           );

figure
tiledlayout(2,2); hAx = nexttile;
for n = 1:4 %looping through each subsystem
	if n > 1, nexttile; end
	
	plot( t, T_integral(:,n) )
	
	set(gca, 'XLimSpec', 'Tight');
	grid on, grid minor
	xlabel('time [s]')
	ylabel('Temperature [\circC]')
end
legend( hAx, 'Integral State',          ...
             'Location','NorthOutside'  );

figure
plot(t,err)
set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
ylabel('Optimization Error')






