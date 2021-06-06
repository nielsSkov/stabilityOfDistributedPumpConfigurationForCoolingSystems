clear all, close all; clc; clear classes                                   %#ok<CLCLS>
                                                                           %#ok<*UNRCH>
run('latexDefaults.m')

%simulation options (select ONE)
simNoSave  = 0;
simAndSave = 1;
loadData   = 0;

%matlab colors
matRed  = [ 0.85 0.325 0.098 ];
matBlue = [ 0    0.447 0.741 ];

figSavePath = 'figures/';

%% Load Model Parameters

run('modelParameters')

%variation in air flow
Q_min  = 1/3*Q;
Q_midt = 2/3*Q;
Q_max  = 3/3*Q;

%variation in ambient temperature
Ta_min  = 25 + 273.15; %[K]
Ta_midt = 30 + 273.15; %[K]
Ta_max  = 35 + 273.15; %[K]

%structure for combination of
% 1: min, 2: midt, 3: max
% between each of the four Q's
%
%create grid of vectors to combine
[v1, v2, v3, v4] = ndgrid( 1:3, 1:3, 1:3, 1:3 );
%stack columns of each sub-grid 
v1 = reshape(v1,[],1);
v2 = reshape(v2,[],1);
v3 = reshape(v3,[],1);
v4 = reshape(v4,[],1);
%combine vectors to arrive at all combinations of 3 values into 4 places
combi = [ v1, v2, v3, v4 ];

%substitute values for Q
v1 = (combi(:,1)==1)*Q_min(1)  + ...
     (combi(:,1)==2)*Q_midt(1) + ...
     (combi(:,1)==3)*Q_max(1);

v2 = (combi(:,2)==1)*Q_min(2)  + ...
     (combi(:,2)==2)*Q_midt(2) + ...
     (combi(:,2)==3)*Q_max(2);

v3 = (combi(:,3)==1)*Q_min(3)  + ...
     (combi(:,3)==2)*Q_midt(3) + ...
     (combi(:,3)==3)*Q_max(3);

v4 = (combi(:,4)==1)*Q_min(4)  + ...
     (combi(:,4)==2)*Q_midt(4) + ...
     (combi(:,4)==3)*Q_max(4);

%create combination matrix of variation over airflows Q
Q_var = [ v1, v2, v3, v4 ];

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
T_final = 30*60; %15*60; %[s]

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
%K = -[ -0.07791562 -0.01685086 -0.00333333 ];         %using LQR

%using LQR
K1 = -[ -0.07607165 -0.00605055 -0.00333333 ];
K2 = -[ -0.07791738 -0.00879464 -0.00333333 ];
K3 = -[ -0.0779497  -0.00986821 -0.00333333 ];
K4 = -[ -0.07698114 -0.00687199 -0.00333333 ];

K_all = [ K1; K2; K3; K4 ];

%initialize flows
qInit = q;

%init pump speeds - used as const pump speeds if control is off (con == 0)
w_init = [ .1 .1 .1 .1 ]';

%initialize vectors for discrete sim
x          = zeros(2,4,length(t));
x1         = zeros(length(Q_var)*3,4,length(t));
x2         = zeros(length(Q_var)*3,4,length(t));
q          = zeros(length(t),4);
p          = zeros(length(t),4);
w          = zeros(length(t),4);
err        = zeros(length(t),1);
T_integral = zeros(length(t),4);

%decide weather to simulate and save (or load results from mat-file)
if simNoSave || simAndSave

%% loops through all combinations for each T_a
T_a = 0;
lenQ = length(Q_var);
jj = 0;
for j = [ 1:lenQ 1:lenQ 1:lenQ ]

progressBar( length(Q_var)*3, 'Running simulations: ' )

jj = jj+1;

if j == 1 && T_a == 0
	T_a = Ta_min;
elseif j == 1 && T_a == Ta_min
	T_a = Ta_midt;
elseif j == 1 && T_a == Ta_midt
	T_a = Ta_max;
end

Q = Q_var(j,:);

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
	
	x1(jj,n,1)  = x(1,n,1);
	x2(jj,n,1)  = x(2,n,1);
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
			
			%select subsystem specific controller
			K = K_all(n,:);
			
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
		
		x1(jj,n,i)  = x(1,n,i);
		x2(jj,n,i)  = x(2,n,i);
	end
end

end

elseif loadData
	load('dataExample2');
end

%save data to mat-file
if simAndSave
	save('dataExample2.mat','x1','x2');
end


%% Plot Results and Save

run('plotExample2.m')






