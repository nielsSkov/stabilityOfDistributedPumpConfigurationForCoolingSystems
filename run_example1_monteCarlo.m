clear all, close all; clc; clear classes                                   %#ok<CLCLS>
<<<<<<< HEAD
<<<<<<< HEAD
oldBar = findall(0,'type','figure','tag','TMWWaitbar'); delete(oldBar);    %#ok<*UNRCH>

=======
                                                                           %#ok<*UNRCH>
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
=======
                                                                           %#ok<*UNRCH>
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
run('latexDefaults.m')

%matlab colors
matRed  = [ 0.85 0.325 0.098 ];
matBlue = [ 0    0.447 0.741 ];

figSavePath = 'figures/';

%simulation options (select ONE)
simNoSave  = 0;
simAndSave = 1;
loadData   = 0;

%% Load Model Parameters

run('modelParameters')

%% Simulation Setup

theta_c = 10 + 273.15; %[K]
T_a     = 30 + 273.15; %[K]

%step size
h = .2; %[s]

%length of simulation
<<<<<<< HEAD
<<<<<<< HEAD
T_final = 30*60; %15*60; %[s]
=======
T_final = 20*60; %15*60; %[s]
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
=======
T_final = 20*60; %15*60; %[s]
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6

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

%control gain vector
%K = -[ -0.07538319  -0.03915181 -0.07705162 ]; %using pole placement
%K = -[ -0.07791562 -0.01685086 -0.00333333 ];  %using LQR

%using LQR
%K1 = -[ -0.07750724 -0.01797215 -0.00333333 ];
%K2 = -[ -0.07281985 -0.02400036 -0.00333333 ];
%K3 = -[ -0.07131688 -0.02596342 -0.00333333 ];
%K4 = -[ -0.0759679  -0.01994393 -0.00333333 ];

<<<<<<< HEAD
<<<<<<< HEAD
K1 = -[ -0.07607165 -0.00605055 -0.00333333 ];
%K2 = -[ -0.07791738 -0.00879464 -0.00333333 ];
%K3 = -[ -0.0779497  -0.00986821 -0.00333333 ];
%K4 = -[ -0.07698114 -0.00687199 -0.00333333 ];
=======
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
K1 = [ -0.07607165 -0.00605055 -0.00333333 ];
%K2 = [ -0.07791738 -0.00879464 -0.00333333 ];
%K3 = [ -0.0779497  -0.00986821 -0.00333333 ];
%K4 = [ -0.07698114 -0.00687199 -0.00333333 ];
<<<<<<< HEAD
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6

%K_all = [ K1; K2; K3; K4 ];

%initialize flows
qInit = q;

%init pump speeds - used as const pump speeds if control is off (con == 0)
w_init = [ .1 .1 .1 .1 ]';

%monte carlo iterations
simIterations = 1000;

%initialize vectors for discrete sim
x          = zeros(2,length(t));
q          = zeros(length(t),1);
p          = zeros(length(t),1);
w          = zeros(length(t),1);
T_integral = zeros(length(t),1);
x_result   = zeros(2,length(t),simIterations);


%% Parameter variation

% Points to simulate
% (Ta min, Q1 min)
% (Ta min, Q1 max)
% (Ta midt, Q1 midt)
% (Ta max, Q1 min)
% (Ta max, Q1 max)

%variation in ambient temperature
Ta_min  = 25 + 273.15; %[K]
Ta_midt = 30 + 273.15; %[K]
Ta_max  = 35 + 273.15; %[K]

%variation in air flow
Q_min  = 1/3*Q;
Q_midt = 2/3*Q;
Q_max  = 3/3*Q;

%random signals for parameter variations bounded by parameter min/max
Q_rand  = Q_min(1) + (Q_max(1)-Q_min(1))*rand(simIterations,1);
Ta_rand = Ta_min(1) + (Ta_max(1)-Ta_min(1))*rand(simIterations,1);

%decide weather to simulate and save (or load results from mat-file)
if simNoSave || simAndSave

<<<<<<< HEAD
<<<<<<< HEAD
for j = 1:simIterations+3

progressBar( simIterations+3, 'Running simulations: ' )

%setting randomized variables
if j <= simIterations
	T_a  = Ta_rand(j);
	Q(1) = Q_rand(j);
elseif j == simIterations+1
	T_a  = Ta_min;
	Q(1) = Q_min(1);
elseif j == simIterations+2 
	T_a  = Ta_midt;
	Q(1) = Q_midt(1);
elseif j == simIterations+3 
	T_a  = Ta_max;
	Q(1) = Q_max(1);
end

%T_a  = 307.24-.0005;
%Q(1) = 8.9262;
%T_a  = 3.072395769621278e+02; Q(1) = 8.9262;
=======
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
for j = 1:simIterations

progressBar( simIterations, 'Running simulations: ' )

%setting randomized variables
T_a  = Ta_rand(j);
Q(1) = Q_rand(j);

<<<<<<< HEAD
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
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

<<<<<<< HEAD
<<<<<<< HEAD
		%calculate return water temperature in equilibrium
		theta_eq = ( -C_a*Q(1)*( T_a - T_eq(1) ) + B(1)*T_eq(1) )/B(1);

		%calculate flow in equilibrium
		q_eq = B(1)*(theta_eq - T_eq(1))/( C_w*(theta_c - theta_eq) );

		%states for calculating controlled pump speed
=======
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
		%calculate flow in equilibrium (T_eq is the desired temperature)
		theta_eq = ( -C_a*Q(1)*( T_a - T_eq(1) ) + B(1)*T_eq(1) )/B(1);

		%calculate return water temperature nn equilibrium
		q_eq = B(1)*(theta_eq - T_eq(1))/( C_w*(theta_c - theta_eq) );

		%states for calculating controlled punp speed
<<<<<<< HEAD
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
=======
>>>>>>> c38b8706fa1f50ba613546da32b91c988cb8f6b6
		y = [ x(1,i-1) - theta_eq  ;
		      x(2,i-1) - T_eq(1)   ;
		      T_integral(i)       ];

		%pump speed (controlled)
		w(i) = K1*y + q_eq;

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

x_result(1:2,:,j) = x(1:2,:);

end

%load data from mat-file
elseif loadData
	load('dataExample1');
end

%save data to mat-file
if simAndSave
	save('dataExample1.mat','x_result');
end

%plot output temperatures (return water and exhaust air)
figure, tiledlayout(2,1), nexttile
for i = 1:simIterations
plot( t, x_result(2,:,i) - 273.15, 'color', matRed  ), hold on
end
ylim([0 30])
set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)
%
nexttile
for i = 1:simIterations
plot( t, x_result(1,:,i) - 273.15, 'color', matBlue ), hold on
end
ylim([0 30])
set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)
%leg = { sprintf('$\\theta_%i$', 1), sprintf('T$_%i$', 1) };
%legend( leg, 'Location','southeast' );

saveCroppedPdf( gcf, [figSavePath 'singleAHU_nonlinManySim_1' '.pdf'] )

sys1 = squeeze(x_result(2,:,:) - 273.15);

min1 = min( sys1, [], 2 );
max1 = max( sys1, [], 2 );

X = [ t'          ;
			flipud(t') ];

Y_sys1 = [ min1          ;
				   flipud(max1) ];

figure
h_fill = fill(X,Y_sys1,[0 .7 0],'edgecolor','none', 'facealpha', '.2');
hold on

min_idx     = simIterations+1;
nominal_idx = simIterations+2;
max_idx     = simIterations+3;

plot( t, sys1(:,min_idx),     'color', '[0 .65 0]', 'linewidth', 1.4 )
plot( t, sys1(:,nominal_idx), 'color', matBlue,     'linewidth', 1.4 )
plot( t, sys1(:,max_idx),     'color', matRed,      'linewidth', 1.4 )

set(gca, 'XLimSpec', 'Tight');
grid on, grid minor
xlabel('time [s]')
labelY = sprintf('Temp. %s', '[$^\circ$C]');
ylabel(labelY)

legend( ['$T_{a,min} \leq T_a \leq T_{a,max}\ ,\ \ '   ...
          'Q_{min}   \leq Q   \leq Q_{max}$'        ], ...
        '$T_{a,min}, Q_{min}$',                            ...
        '$T_{a,nom}, Q_{nom}$',                            ...
        '$T_{a,max}, Q_{max}$'                             )


saveCroppedPdf( gcf, [figSavePath 'singleAHU_nonlinManySim_2' '.pdf'] )


oldBar = findall(0,'type','figure','tag','TMWWaitbar'); delete(oldBar);


