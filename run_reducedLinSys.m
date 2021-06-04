clear all, close all; clc

run('latexDefaults.m')

figSavePath = 'figures/';

%% Load Model Parameters

run('modelParameters')

%setting cooling water and ambient air temperatures
theta_c = 10 + 273.15; %[K]
T_a     = 30 + 273.15; %[K]

%sample time
Ts = 0.1; %[s]

%length of simulation
t_final = 2*60; %[s]

%initial states
theta_0 = T_a;%0 + 273.15; %[K]
T_0     = T_a;%0 + 273.15; %[K]
zeta_0  = 0;

%initialization for ode15s
tspan = 0:Ts:t_final;
init = [ theta_0  T_0  zeta_0 ];

%option to lower relative tollerence (default 1e-3)
options = odeset('RelTol',1e-3);

%simulate and plot the iTh WAHE reduced linear system
for iTh = 1:4
	
	%pole placement
	%K = -[ -0.07538319 -0.03915181 -0.07705162 ]; %poles: [ -1 -2 -3 ]
	K = -[ -0.00242309 -0.06580131 -0.03082065 ]; %poles: [ -1, -1.5, -1.6 ]
	%K =  -[ 0.07821702 -0.11822958 -0.00269681 ]; %poles: [ -0.3, -0.7, -1 ]
	%LQR
	%K = -[ -0.07791562 -0.01685086 -0.00333333 ];
	
	[ t, x ] = ode15s( @(t,x) reducedLinSys( t, x, K, r, R, R_c, a, b, ...
	                                         iTh, theta_c, T_a, Ts,    ...
	                                         B, C_w, C_a, V_w, V_a, Q  ), ...
	                     tspan, init, options                             );
	
	%initializing flow and motor speed vectors
	q = zeros(length(t),1);
	w = zeros(length(t),1);

	%re-running sim-function in loop to extract flows at each time step
	for j = 1:length(t)
  [ ~, q(j), ...
       w(j)  ] = reducedLinSys( t(j), x(j,:)', K, r, R, R_c, a, b, ...
	                              iTh, theta_c, T_a, Ts,             ...
	                              B, C_w, C_a, V_w, V_a, Q           );
	end
	
	%saving sim result for each subsystem
	x_out.(['sub' num2str(iTh)]) = x;
	t_out.(['sub' num2str(iTh)]) = t;
	q_out.(['sub' num2str(iTh)]) = q;
	w_out.(['sub' num2str(iTh)]) = w;
		
end

%plot output temperatures (return water and exhaust air)
figure
margY = 3;
tiledlayout(2,2)
hAx = nexttile;
plot( t_out.sub1, x_out.sub1(:,1)-273.15 ), hold on
plot( t_out.sub1, x_out.sub1(:,2)-273.15 )
ylim([0 30])
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
legend( '$\theta_1$', 'T$_1$', 'Location','southeast' );
xlabel('time [s]')
ylabel('Temp. [$^\circ$C]')
nexttile
plot( t_out.sub2, x_out.sub2(:,1)-273.15 ), hold on
plot( t_out.sub2, x_out.sub2(:,2)-273.15 )
ylim([0 30])
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
xlabel('time [s]')
ylabel('Temp. [$^\circ$C]')
legend( '$\theta_2$', 'T$_2$', 'Location','southeast' );
nexttile
plot( t_out.sub3, x_out.sub3(:,1)-273.15 ), hold on
plot( t_out.sub3, x_out.sub3(:,2)-273.15 )
ylim([0 30])
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
xlabel('time [s]')
ylabel('Temp. [$^\circ$C]')
legend( '$\theta_3$', 'T$_3$', 'Location','southeast' );
nexttile
plot( t_out.sub4, x_out.sub4(:,1)-273.15 ), hold on
plot( t_out.sub4, x_out.sub4(:,2)-273.15 )
ylim([0 30])
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
xlabel('time [s]')
ylabel('Temp. [$^\circ$C]')
legend( '$\theta_4$', 'T$_4$', 'Location','southeast' );
saveCroppedPdf( gcf, [figSavePath 'polePlacementIsolatedSystems' '.pdf'] )

figure
plot( t_out.sub4, x_out.sub4(:,1)-273.15 ), hold on
plot( t_out.sub4, x_out.sub4(:,2)-273.15 )
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
xlabel('time [s]')
ylabel('Temp. [$^\circ$C]')
xlim([ 0 200 ])
limY = ylim; ylim([ 0 limY(2) ])
legend( 'Return Water Temperature, $\theta$', ...
        'Exhaust Air Temperature, T',         ...
        'Location','northeast   '             );

saveCroppedPdf( gcf, [figSavePath 'tempSingleUnitPolePlacement' '.pdf'] )

%plot water flow
figure
plot( t_out.sub1, q_out.sub1(:,1) ), hold on
plot( t_out.sub2, q_out.sub2(:,1) )
plot( t_out.sub3, q_out.sub3(:,1) )
plot( t_out.sub4, q_out.sub4(:,1) )
grid on, grid minor
xlabel('time [s]')
ylabel('Flow [$m^3/h$]')

%plot pump speeds
figure
plot( t_out.sub1, w_out.sub1(:,1) ), hold on
plot( t_out.sub2, w_out.sub2(:,1) )
plot( t_out.sub3, w_out.sub3(:,1) )
plot( t_out.sub4, w_out.sub4(:,1) )
grid on, grid minor
set(gca, 'XLimSpec', 'Tight');
xlabel('time [s]')
ylabel('Pump Speed, $\omega$')
legend( '$\omega_1$', '$\omega_2$','$\omega_3$', '$\omega_4$', ...
				'Location', 'northeast'                                )


%% Contorl Design

theta_c = 10 + 273.15; %[K]
T_a     = 30 + 273.15; %[K]

for i = 1:4
	
	T_eq = 16 + 273.15; %[K]

	%calculate flow in equilibrium (T_eq is the desired temperature)
	theta_eq = ( -C_a*Q(i)*( T_a - T_eq ) + B(i)*T_eq )/B(i);

	%calculate return water temperature in equilibrium
	q_eq = B(i)*(theta_eq - T_eq)/( C_w*(theta_c - theta_eq) );

	%AWHE System Matrix (with two integral states)
	A_lin = [ ...
	-(q_eq/V_w(i) + B(i)/(C_w*V_w(i)))  B(i)/(C_w*V_w(i))                0  ;
		B(i)/(C_a*V_a(i))               -(Q(i)/V_a(i) + B(i)/(C_a*V_a(i))) 0  ;
		0                                 1                                0 ];

	%AWHE Input Matrix
	B_lin = [ (theta_c - theta_eq)/V_w(i)  ;
						 0                           ;
						 0                          ];
	
	%calculate controllability matrix
	calC = [ B_lin  A_lin*B_lin  A_lin^2*B_lin  ];

	
	fprintf('\nRank of Controllability Matrix is %i\n', rank(calC))
	disp('A = '),    disp(A_lin)
	disp('B = '),    disp(B_lin)
	disp('calC = '), disp(calC)
	
	%rref(calC)

	%K = place(A,B,p)
	%K = [ 0.0813625  -0.20673981 -0.00867699] ];
	
end