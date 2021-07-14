clear all, close all; clc; clear classes                                   %#ok<CLCLS>
                                                                           %#ok<*UNRCH>
run('latexDefaults.m')

%matlab colors
matRed  = [ 0.85 0.325 0.098 ];
matBlue = [ 0    0.447 0.741 ];

figSavePath = 'figures/';

%% Load Model Parameters

run('modelParameters')

%% Simulation Setup

theta_c = 10 + 273.15; %[K]
T_a     = 30 + 273.15; %[K]

%step size
h = .2; %[s]

%length of simulation
T_final = 60*60;  %[s]

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
T_eq = ones(1,4)*(16 + 273.15);

%using LQR
K1 = -[ -0.07607165 -0.00605055 -0.00333333 ];

%initialize flows
qInit = q;

%init pump speeds - used as const pump speeds if control is off (con == 0)
w_init = [ .1 .1 .1 .1 ]';

%initialize vectors for discrete sim
x          = zeros(2,length(t));
q          = zeros(length(t),1);
p          = zeros(length(t),1);
w          = zeros(length(t),1);
q_eq       = zeros(length(t),1);
T_integral = zeros(length(t),1);
x_result   = zeros(2,length(t),1);

for j = 3

if j == 1 %case 1
	T_a  = 307.2395769621278;
	Q(1) = 8.9262;
elseif j == 2 %case 2
	T_a  = 307.24;
	Q(1) = 8.926;
elseif j == 3 %case 3
	T_a  = 307.2395;
	Q(1) = 8.9262;
end

%T_a  = 307.24-.0005;
%Q(1) = 8.9262;

%% First Step of Discrete Sim

q(1) = sqrt( b(1)*w_init(1)^2 / ( r(1) + R_c + 2*R(1) + a(1) ) );

x(:,1) = init(1,:)';


%% Discrete Simulation - Euler's Method
for i = 2:length(t)

	if conOn
		if i > 2
			T_errThen = x(2,i-2) - T_eq(1);
		else
			T_errThen = 0;
		end

		T_errNow  = x(2,i-1) - T_eq(1);

		T_integral(i) = T_integral(i-1) + h*( T_errNow + T_errThen )/2;
		
		%calculate return water temperature in equilibrium
		theta_eq = ( -C_a*Q(1)*( T_a - T_eq(1) ) + B(1)*T_eq(1) )/B(1);

		%calculate flow in equilibrium
		q_eq(i) = B(1)*(theta_eq - T_eq(1))/( C_w*(theta_c - theta_eq) );

		%states for calculating controlled pump speed
		y = [ x(1,i-1) - theta_eq  ;
		      x(2,i-1) - T_eq(1)   ;
		      T_integral(i)       ];

		%pump speed (controlled)
		w(i) = K1*y + q_eq(i);

		%set minimum allowed pump speed
		w_min  = 0.01;

		%option for hard upper limmit of pump speed
		w_max = inf;

		%limit pump speed
		if     w(i) < w_min, w(i) = w_min;
		elseif w(i) > w_max, w(i) = w_max;
		end
	else
		w(i,:) = w_init;
	end

	q(i) = sqrt( b(1)*w(i)^2 / ( r(1) + R_c + 2*R(1) + a(1) ) );

	A = [ ...
	 -(q(i)/V_w(1) +B(1)/(C_w*V_w(1)))   B(1)/(C_w*V_w(1))                ;
	   B(1)/(C_a*V_a(1))               -(Q(1)/V_a(1) +B(1)/(C_a*V_a(1))) ];

	f = A*x(:,i-1) + [ (q(i)/V_w(1))*theta_c  ;
	                   (Q(1)/V_a(1))*T_a     ];

	x(1:2,i) = x(1:2,i-1) + h*f;
end

%plot output temperatures (return water and exhaust air)
figure, tiledlayout(2,1), nexttile
plot( t, x(2,:) - 273.15, 'color', matRed  ), hold on
set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)
%
nexttile
plot( t, x(1,:) - 273.15, 'color', matBlue ), hold on
set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)
%leg = { sprintf('$\\theta_%i$', 1), sprintf('T$_%i$', 1) };
%legend( leg, 'Location','southeast' );

if j == 1 %case 1
	saveCroppedPdf( gcf, [figSavePath 'case1' '.pdf'] )
elseif j == 2 %case 2
	saveCroppedPdf( gcf, [figSavePath 'case2' '.pdf'] )
elseif j == 3 %case 3
	saveCroppedPdf( gcf, [figSavePath 'case3' '.pdf'] )
end

figure
plot(t,w), hold on
grid on, grid minor
xlabel('time [s]')
ylabel('Pump Speed, $\omega$')

figure
plot(t(2:end),q_eq(2:end)), hold on
grid on, grid minor
xlabel('time [s]')
ylabel('$q_{eq}$')

end