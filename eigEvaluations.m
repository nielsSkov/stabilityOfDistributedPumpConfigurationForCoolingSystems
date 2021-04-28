clear; 
%close all; 
clc;

% Load model parameters
run('hydraulicCircuit')

% Setup
T_a = 30 + 273.15+5;                  % Ambient temperature [K]
T_r = ones(1,4)*(20 + 273.15);
% T_r = [293 294 295 296];    % AHU reference air temperatures [K]
th_c = 10 + 273.15;                 % AHU inflow water temperature [K] 
Ki = -[ -0.07538319 -0.03915181 -0.07705162 ]*1;

Q = Q*0.5;
% Compute \theta^*_i and q^*_i eq. (7) and (8)
th_r = zeros(n,1);
q_r = zeros(n,1);
As = zeros(3,3,n);
Bs = zeros(3,1,n);
F_H_q = zeros(3*n,3*n);
G_H = zeros(3*n,n);
alpha_bar = zeros(n,1);
rho = zeros(n,1);
Df_q = zeros(n,n);

for i=1:n
    % (7) and (8)
    th_r(i) = (-C_a*Q(i)*(T_a-T_r(i)) + B(i)*T_r(i))/B(i);
    q_r(i) = (B(i)*(th_r(i)-T_r(i)))/(C_w*(th_c-th_r(i)));
    
    % (10) and (11)
    As(:,:,i) = [-(q_r(i)/V_w(i)+B(i)/(C_w*V_w(i)))     B(i)/(C_w*V_w(i))                      0;
                 B(i)/(C_a*V_a(i))                         -(Q(i)/V_a(i)+B(i)/(C_a*V_a(i)))    0;
                 0                                         1                                   0];
    
    Bs(:,:,i) = [(th_c-th_r(i))/V_w(i); 0; 0];
    % Check controllability
    rank(ctrb(As(:,:,i),Bs(:,:,i)));
    
    % Compute (27) and (28)
    F_H_q((i-1)*3+1:i*3,(i-1)*3+1:i*3) = As(:,:,i);
    G_H((i-1)*3+1:i*3,i) = Bs(:,:,i);
    
    % Compute alpha_bar
    tt = 0;
    for jj=1:i
        tt = tt + 2*R(jj)/b(i);
    end
    alpha_bar(i) = sqrt((r(i)+a(i))/b(i) + R_c/b(i) + tt);
end

Lambda_bar = diag(alpha_bar);

for i=1:n
    % Compute rho
    rho(i) = 1/(2*sqrt(q_r'*S(:,:,i)*q_r));
    
    % Compute Df_q
    Df_q(:,i) = rho(i)*S(:,:,i)*q_r;
end
Dg_w = inv(Df_q');


%% Decoupled system eq. (26b)
K = blkdiag(Ki,Ki,Ki,Ki);
x0i = [1; -1; 0];
x0 = repmat(x0i,4,1);
sysCL_de = ss(F_H_q+G_H*K,zeros(3*n),eye(3*n),0);

eig(sysCL_de)


%% Coupled system eq. (37b)
sysCL_co = ss(F_H_q+G_H*Dg_w*Lambda_bar*K, zeros(3*n),eye(3*n),0);

eig(sysCL_co)

figure
initial(sysCL_de,x0)
hold on
initial(sysCL_co,x0,'r')
legend('Decoupled','Coupled')

figure
pzmap(sysCL_de)
hold on
pzmap(sysCL_co, 'r')
legend('Decoupled','Coupled')
