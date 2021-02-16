function [ x_dot, q, w ] = reducedLinSys_n3( ~, x, K, r, R, R_c, a, b, ...
                                             i, theta_c, T_a, Ts,      ...
                                             B, C_w, C_a, V_w, V_a, Q  )

T_eq = 20 + 273.15;

%calculate flow in equilibrium (T_eq is the desired temperature)
theta_eq = ( -C_a*Q(i)*( T_a - T_eq ) + B(i)*T_eq )/B(i);

%calculate return water temperature in equilibrium
q_eq = B(i)*(theta_eq - T_eq)/( C_w*(theta_c - theta_eq) );

%states
y = [ x(1) - theta_eq  ;
      x(2) - T_eq      ;
      x(3)            ];

%calculate alpha
alpha = sqrt( ( r(i) + R_c + 2*sum(R(1:i)) + a(i) )/b(i) );

%pump speed (controlled)
w = alpha*K*y + alpha*q_eq;

w_min = 0.01;
w_max = inf;  %1;
if     w < w_min, w = w_min;
elseif w > w_max, w = w_max;
end

%calculate flow, q, given pump speed, w, for iTh AWHE's hydraulic circuit
q = (1/alpha)*w;

%AWHE System Matrix (with two integral states)
lin_A = [ ...
  -(q_eq/V_w(i) + B(i)/(C_w*V_w(i)))  B(i)/(C_w*V_w(i))                0  ;
    B(i)/(C_a*V_a(i))               -(Q(i)/V_a(i) + B(i)/(C_a*V_a(i))) 0  ;
    0                                 1                                0 ];

lin_B = [ (theta_c - theta_eq)/V_w(i)  ;
           0                           ;
           0                          ];

u = (q - q_eq);

x_dot = lin_A*y + lin_B*u;

end

