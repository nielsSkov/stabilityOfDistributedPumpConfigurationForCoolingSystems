clear; 
close all; 
clc;

% Load model parameters
run('hydraulicCircuit')

% Setup
T_r = ones(1,4)*(20 + 273.15);
% T_r = [293 294 295 296];    % AHU reference air temperatures [K]
th_c = 10 + 273.15;                 % AHU inflow water temperature [K] 
% Ki = -[ -0.07538319 -0.03915181 -0.07705162 ]*1;
Ki = -[ -0.07791562 -0.01685086 -0.00333333 ];

nn = 3; % number of paramter values

% Allocate memory
th_r = zeros(n,1);
q_r = zeros(n,1);
As = zeros(3,3,n);
Bs = zeros(3,1,n);
F_H_q = zeros(3*n,3*n);
G_H = zeros(3*n,n);
alpha_bar = zeros(n,1);
rho = zeros(n,1);
Df_q = zeros(n,n);



% Set flows
Q = Q;

gran = 2;
% Change parameters loop
Ta = linspace(T_aMin,T_aMax,gran);
Qq = [linspace(Q(1)/3,Q(1),gran); 
      linspace(Q(2)/3,Q(2),gran);
      linspace(Q(3)/3,Q(3),gran);
      linspace(Q(4)/3,Q(4),gran)];

Ta_nom = T_aMin+(T_aMax-T_aMin)/2;
Qq_nom = Qq(:,1)+(Qq(:,2)-Qq(:,1))./2;

nn = length(Ta);
nnn = n^nn;

eigDecou = zeros(3*n,nn);
eigCou = zeros(3*n,nn);

%% Simulate nominal system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n % number of systems
    % (7) and (8)
    th_r(i) = (-C_a*Qq_nom(i)*(Ta_nom-T_r(i)) + B(i)*T_r(i))/B(i);
    q_r(i) = (B(i)*(th_r(i)-T_r(i)))/(C_w*(th_c-th_r(i)));

    % (10) and (11)
    As(:,:,i) = [-(q_r(i)/V_w(i)+B(i)/(C_w*V_w(i)))     B(i)/(C_w*V_w(i))                      0;
                 B(i)/(C_a*V_a(i))                      -(Qq_nom(i)/V_a(i)+B(i)/(C_a*V_a(i)))    0;
                 0                                      1                                      0];

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

for ii=1:n
    % Compute rho
    rho(ii) = 1/(2*sqrt(q_r'*S(:,:,ii)*q_r));

    % Compute Df_q
    Df_q(:,ii) = rho(ii)*S(:,:,ii)*q_r;
end
Dg_w = inv(Df_q');

% Decoupled system eq. (26b)
K = blkdiag(Ki,Ki,Ki,Ki);
x0i = [1; -1; 0];
x0 = repmat(x0i,4,1);
sysCL_de_nom = ss(F_H_q+G_H*K,zeros(3*n),eye(3*n),0);

eigDecou_nom = eig(sysCL_de_nom);

% Coupled system eq. (37b)
sysCL_co_nom = ss(F_H_q+G_H*Dg_w*Lambda_bar*K, zeros(3*n),eye(3*n),0);

eigCou_nom = eig(sysCL_co_nom);

%% Eigenvalue numerical analysis
% Form the Qq object that has 16 different configurations
Qqq = fliplr(combvec(Qq(4,:),Qq(3,:),Qq(2,:),Qq(1,:))')';
for vv = 1:nn % Ta
    for j=1:nnn % Qq 
        Qq_t = Qqq(:,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:n % number of systems
            % (7) and (8)
            th_r(i) = (-C_a*Qq_t(i)*(Ta(vv)-T_r(i)) + B(i)*T_r(i))/B(i);
            q_r(i) = (B(i)*(th_r(i)-T_r(i)))/(C_w*(th_c-th_r(i)));

            % (10) and (11)
            As(:,:,i) = [-(q_r(i)/V_w(i)+B(i)/(C_w*V_w(i)))     B(i)/(C_w*V_w(i))                      0;
                         B(i)/(C_a*V_a(i))                      -(Qq_t(i)/V_a(i)+B(i)/(C_a*V_a(i)))    0;
                         0                                      1                                      0];

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

        for ii=1:n
            % Compute rho
            rho(ii) = 1/(2*sqrt(q_r'*S(:,:,ii)*q_r));

            % Compute Df_q
            Df_q(:,ii) = rho(ii)*S(:,:,ii)*q_r;
        end
        Dg_w = inv(Df_q');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %% Decoupled system eq. (26b)
        K = blkdiag(Ki,Ki,Ki,Ki);
        x0i = [1; -1; 0];
        x0 = repmat(x0i,4,1);
        sysCL_de(:,:,vv,j) = ss(F_H_q+G_H*K,zeros(3*n),eye(3*n),0);
        % Only first subsystem
        sysCL_de_1(:,:,vv,j) = ss(sysCL_de(:,:,vv,j).a(1:3,1:3),sysCL_de(:,:,vv,j).b(1:3,1:3),sysCL_de(:,:,vv,j).c(1:3,1:3),sysCL_de(:,:,vv,j).d(1:3,1:3));

        eigDecou(:,vv,j) = eig(sysCL_de(:,:,vv,j));
        


        %% Coupled system eq. (37b)
        sysCL_co(:,:,vv,j) = ss(F_H_q+G_H*Dg_w*Lambda_bar*K, zeros(3*n),eye(3*n),0);

        eigCou(:,vv,j) = eig(sysCL_co(:,:,vv,j));
    end
end

%% Find maximum real value of all eigenvalues
maxRealEigVal = max(real(eigCou(:)));

%% Plotting
% Plot nominal system
figure(1)
jj=1;
for vv=1:nn
    for j=1:nnn
%         if vv==i || (vv==1 && i==3) || (vv==3 && i==1)
%             if vv==1
%                 textTa = '\overline{T}_a';
%             elseif vv==2
%                 textTa = 'T_a^*';
%             elseif vv==3
%                 textTa = '\underline{T}_a';
%             end
%             if i==1
%                 textQ1 = '\overline{Q}_1$';
%             elseif i==2
%                 textQ1 = 'Q_1^*$';
%             elseif i==3
%                 textQ1 = '\underline{Q}_1$';
%             end
%             legendT{jj} = ['$T_a = ' textTa ',\; Q_1 = ' textQ1];
%             jj = jj+1;
            plot(real(eigCou(:,vv,j)),imag(eigCou(:,vv,j)),'xb');
            hold on
%             legend(legendT, 'interpreter', 'latex','Location','NorthWest')
%         end
    end
end
xlabel('Real Axis')
ylabel('Imaginary Axis')

figure(2)
plot(real(eigCou(end,:)), imag(eigCou(end,:)),'-x')

figure(3)
pzmap(sysCL_co, 'r')

figure(4)
pzmap(sysCL_de_nom)
hold on
pzmap(sysCL_co_nom, 'r')
legend('Decoupled System','Coupled System','Location','northwest')
% figure(4)
% plot(real(eigDecou_nom),imag(eigDecou_nom),'xb',real(eigCou_nom),imag(eigCou_nom),'xr')
% xlabel('Real Axis')
% ylabel('Imaginary Axis')
% ylim([-2 2])
% grid
% legend('Decoupled','Coupled')
% title('Pole Placement of Nominal System')
% ,'Location','NorthOutSide','Orientation','Horizontal')
% legend boxoff

%% Plotting only first subsystem
% figure(5)
% jj = 1;
% for vv=1:nn
%     hold all
%     for i=1:nn
%         if vv==i || (vv==1 && i==3) || (vv==3 && i==1)
%             [y,t,~] = initial(sysCL_de_1(:,:,vv,i),ones(3,1),150);
%             subplot(3,1,1)
%             plot(t,y(:,1))
%             ylabel('$x_{1,1}$', 'interpreter', 'latex')
%             if vv==1
%                 textTa = '\overline{T}_a';
%             elseif vv==2
%                 textTa = 'T_a^*';
%             elseif vv==3
%                 textTa = '\underline{T}_a';
%             end
%             if i==1
%                 textQ1 = '\overline{Q}_1$';
%             elseif i==2
%                 textQ1 = 'Q_1^*$';
%             elseif i==3
%                 textQ1 = '\underline{Q}_1$';
%             end
%             legendText{jj} = ['$T_a = ' textTa ',\; Q_1 = ' textQ1];
%             jj = jj+1;
%             legend(legendText, 'interpreter', 'latex')
%             xlim([0 5])
%             hold all
%             subplot(3,1,2)
%             plot(t,y(:,2))
%             ylabel('$x_{1,2}$', 'interpreter', 'latex')
%             xlim([0 5])
%             hold all
%             subplot(3,1,3)
%             plot(t,y(:,3))
%             ylabel('$x_{1,3}$', 'interpreter', 'latex')
%             hold all
%             xlabel('Time [sec]')
%             
%         end
%     end
% end
