%% Air to Water Heat Exchanger Parameters

%number of AHU
n = 4;

%ambient temperature
T_a = 30 + 273.15; % [K]
T_aMax = 35 + 273.15;
T_aMin = 25 + 273.15;

%heat transfer constant between water and air
B = [ 24000.000  14000.000  12000.000  20000.000 ]; %[J/K]

%specific heat coff for water (Cw) and air (Ca)
C_w = 4170451.0; %[J/m3]
C_a = 891.8;     %[J/m3]

%water volumes
V_w = [ 0.230  0.134  0.115  0.192 ]; %[m^3]

%air volumes
V_a = [ 9.0  5.2  4.5  7.5 ]; %[m^3]

%water flows
q = [ 20.7 12.1 10.4 17.3 ]*1/3600; %[m^3/h]*1/3600 = [m^3/s]

%air flows
Q = [ 32294.2  18838.3  16147.1  26911.9 ]*1/3600;%[m^3/h]*1/3600 = [m^3/s]
Qmin = Q*(1/3);

%  singular pertubated systems (S 430)

%delay estimates (Duck volume 6.8 [m^3])
%air flow: 8097 [m3/h]  ->  delay 3.0 [sec]
%air flow: 12146 [m3/h]  ->  delay 2.0 [sec]


%% Hydraulic Circuit Model Parameters

%calculating guess of R (hydraulic resistance)
%assuming r1=R1=R_c and pressure change deltaP1=0.15
%R = 0.15/( q(1)^2 + (q(1)+q(2)+q(3)+q(4))^2 + 2*(q(1)+q(2)+q(3)+q(4))^2 )
%R = 1.3147e-05
% 
% R_guess = 1.3147e-05;
% 
% R = ones(1,4);
% r = ones(1,4);
% d_p = zeros(1,4); % [bar]
% 
% R(1) = R_guess+.1*R_guess;  R(2) = R_guess+.2*R_guess;
% R(3) = R_guess-.2*R_guess;  R(4) = R_guess-.1*R_guess;
% 
% r(1) = R_guess-.2*R_guess;  r(2) = R_guess-.24*R_guess;
% r(3) = R_guess-.28*R_guess; r(4) = R_guess-.3*R_guess;
% 
% R_c = R_guess+.3*R_guess;
% 
% %pump parameters (guessed)
% a = [ 1 1 1 1 ];
% b = [ 1 1 1 1 ];

%Dymola parameters
Pzones = [120e3, 70e3, 60e3, 100e3]; %[J] Room power 
Dp_transport = 5;                    %[m] pressure drop over pipes
Dp_heatexchanger = 5;                %[m] pressure drop over heat exchanger
n_num = 1;                           %[-] nominal pump speed

% Pump and pipe parameters

R = zeros(1,n);
a = zeros(1,n);
b = zeros(1,n);
for i=1:length(Pzones)
  R(i) = Dp_transport/(2*sum(q(i:end))^2);
  b(i) = 2*(2*i*Dp_transport+Dp_heatexchanger)/(n_num^2);
  a(i) = (2*i*Dp_transport+Dp_heatexchanger)/((q(i))^2);
	r(i) = Dp_heatexchanger/(q(i)^2);
end

clear i

%cooler hydraulic resistance (guessed)
R_c = Dp_heatexchanger/(sum(q(:))^2);

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

%% create python model parameter variables
pth.q = py.list(q);
pth.Q = py.list(Q);
pth.R = py.list(R);
pth.r = py.list(r);
pth.a = py.list(a);
pth.b = py.list(b);



