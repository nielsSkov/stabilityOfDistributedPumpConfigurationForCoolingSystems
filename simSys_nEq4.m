function [ x_dot,   ...
           q, p, w, ...
           w_err    ] = simSys_nEq4( ~, x, w_init, T_eq, r, R, R_c, ...
                                     a, b,theta_c, T_a, K, conOn,   ...
                                     B, C_w, C_a, V_w, V_a, Q, pth  )        
persistent qPrev                                                           %#ok<*UNRCH>

if isempty(qPrev)
	qPrev = pth.q;
end

%initialize w to w_init, the values set here are used only if conOn == 0
w = w_init;

%initialize integral anti windup
antiWind = zeros(1,4);

if conOn
	for i = 1:4
		%calculate flow in equilibrium (T_eq is the desired temperature)
		theta_eq = ( -C_a*Q(i)*( T_a - T_eq(i) ) + B(i)*T_eq(i) )/B(i);
		
		%calculate return water temperature in equilibrium
		q_eq = B(i)*(theta_eq - T_eq(i))/( C_w*(theta_c - theta_eq) );
		
		%set minimum allowed pump speed
		w_min  = 0.01;
		
		%require integral state min-limit that keeps pump speed above w_min
		x1 = x(3*i-2) - theta_eq;
		x2 = x(3*i-1) - T_eq(i);
		x3_min = (w_min - K(1)*x1 - K(2)*x2 + q_eq)/K(3);
		
		%option for hard max-limit of integral state
		x3_max = inf;
		
		%calculate integral saturation
		if     x(3*i) < x3_min,  x3_sat = x3_min;
		elseif x(3*i) > x3_max,  x3_sat = x3_max;
		else,                    x3_sat = x(3*i);
		end
		
		%integral anti windup saturation block feedback gain
		k_wind = 2;
		
		%integral anti windup (implemented in final states)
		antiWind(i) = k_wind*( x(3*i) - x3_sat );
		
		%states for calculating controlled punp speed
		y = [ x(3*i-2) - theta_eq  ;
		      x(3*i-1) - T_eq(i)   ;
		      x(3*i)              ];
		
		%pump speed (controlled)
		w(i) = K*y + q_eq;
		
		%option for hard upper limmit of pump speed
		w_max = inf;
		
		%limit pump speed
		if     w(i) < w_min, w(i) = w_min;
		elseif w(i) > w_max, w(i) = w_max;
		end
	end
end

%convert pump speed to a format accepted by python
pth.w = py.list(w);

%solve using the Broyden–Fletcher–Goldfarb–Shanno algorithm
q_out = py.pyMinimize.BFGS( qPrev, pth.Q, pth.w, pth.R, ...
                            pth.r, pth.a, pth.b, R_c    );

q = [ q_out{1} q_out{2} q_out{3} q_out{4} ]';

%save flows as persistant variable for initial guess in optimization
qPrev = py.list(q);

%calculate optimization error in therms of pump speed, w
w_postOpt(1) = ...
	sqrt(((r(1) + a(1))*q(1)^2 + R_c*sum(q)^2 + 2*R(1)*sum(q)^2)/b(1) );

w_postOpt(2) = ...
	sqrt(((r(2) + a(2))*q(2)^2 + R_c*sum(q)^2 + 2*R(1)*sum(q)^2   + ...
                                              2*R(2)*sum(q(2:4))^2)/b(2));
w_postOpt(3) = ...
	sqrt(((r(3) + a(3))*q(3)^2 + R_c*sum(q)^2 + 2*R(1)*sum(q)^2   + ...
                                              2*R(2)*sum(q(2:4))^2   + ...
                                              2*R(3)*sum(q(3:4))^2)/b(3));
w_postOpt(4) = ...
	sqrt(((r(4) + a(4))*q(4)^2 + R_c*sum(q)^2 + 2*R(1)*sum(q)^2 + ...
                                              2*R(2)*sum(q(2:4))^2 + ...
                                              2*R(3)*sum(q(3:4))^2 + ...
                                              2*R(4)*q(4)^2)/b(4)  );
w_err     = w-w_postOpt';

%calculate pressure
p = a.*q'.^2 + b.*w'.^2;

%% Air to Water Heat Exchangers
for i = 1:4
	A = [ ...
		-(q(i)/V_w(i) +B(i)/(C_w*V_w(i)))  B(i)/(C_w*V_w(i))               0  ;
		  B(i)/(C_a*V_a(i))              -(Q(i)/V_a(i) +B(i)/(C_a*V_a(i))) 0  ;
		  0                                0                               1 ];
	
	%building state vector including anti windup integral state
	xx = [ x(3*i-2)                          ;
		     x(3*i-1)                          ;
		     x(3*i-1) - T_eq(i) - antiWind(i) ];
	
	x_dot(3*i-2:3*i) = A*xx + [ (q(i)/V_w(i))*theta_c  ;
	                            (Q(i)/V_a(i))*T_a      ;
	                             0                     ];
end

x_dot = x_dot';

end

