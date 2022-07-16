function [x, P] = Giorgi_Kalman(timestamps, x1, P1, var_C, var_d, var_theta, var_gamma)

%%  Giorgi's Kalman Filtering
%   This simple source code implements the Giorgi's Kalman Filtering proposed in
%   G. Giorgi and G. Narduzzi, 
%   "Performance Analysis of Kalman-Filter-Based Clock Synchronization in IEEE 1588 Networks," 
%   IEEE Transactions on Instrumentation and Measurement, vol. 60, no. 8, pp. 2902–2909, Aug. 2011.
%
%%  Assumption of this implementation  
%   - Local clock is running without modification: u(n) = [0 0]'.
%
%%  State Space Model defined by Eq.(23) and (26)
%   z(n) = Hx(n) + v(n), v(n) ~ N(0, R(n))
%   x(n) = A(n)x(n-1) + B(n)u(n-1) + w(n-1), w(n) ~ N(0, Q(n)) by Eq.(23)
%
%   z: observation,  v: measurement noise     
%   x: state,        u: control input,     w: process noise 
% 
%   Eq.(24) and Eq.(25) define A(n) and B(n)
%   : A(n) = [1 dt(n)   B(n) = [-1 -dt(n)
%             0    1 ]           0    -1 ]
%
%   Eq.(26) assumes the H is identity matrix
%   : H = [1 0   
%          0 1]
%
%   Eq.(33) defines R(n)
%   : R(n) = [var_theta_M         var_theta_M/dt(n)
%             var_theta_M/dt(n) 2(var_theta_M/dt(n))^2]
%
%   Line 4 of right column in page 2906 defines Q(n)
%   : Q(n) = [var_theta*dt(n)             0
%                          0 var_gamma*dt(n)]
%
%   Eq.(25) assumes that 
%    dt corresponds to the time interavl between two IEEE 1588 sync messages.
%   dt(n) = T1(n) - T1(n-1)
%
%   Eq.(18) defined var_theta_M
%   : var_theta_M = 0.5*(var_C + var_T + var_d)
%     - var_T = 0 : the master node is assumed to 
%                   have negligible timestamping uncertainty
%     - var_C : timestamping uncertainty of the slave node
%     - var_d : asymmetry uncertainty 
%
%%  Kalman Filtering
%   
%   By Eq.(27) and Eq.(28)
%   : x(n|n-1) ~ N(x_hat(n|n-1), P(n|n-1))
%     - xhat(n|n-1) = A(n)xhat(t-1) + Bu(n-1)
%     -    P(n|n-1) = A(n)P(n-1)A(n)' + Q(n)
%
%   By Eq.(29) and Eq.(30)
%   : K(t) = P(n|n-1)inv(P(n|n-1)+R(n))
%
%   By Eq.(31) and Eq.(32)
%   : xhat(n) = xhat(n|n-1)+K(n)(z_M(n)-xhat(n|n-1))
%        P(n) = [I-K(n)]P(n|n-1)
%
%   By Eq.(16) and Eq.(17) 
%   : z_M(n) = [theta_M(n) gamma_M(n)]'
%
%% Implementation
%
%  Input: timestamps, x1, P1, var_T, var_C, var_d, var_theta, var_gamma
%   - timestamps: IEEE 1588 PTP timestamps [t1 C(t2) C(t3) t4]
%   - x1: initial value of x, x(1)
%   - P1: initial matrix of P, P{1}
%   - var_C: timestamping uncertainty of the slave
%   - var_d: uncertainty of the asymmetry
%   - var_theta: process noise of the offset  
%   - var_gamma: process noise of the skew
%  
%  Output: x, P 
%   - x: set of [offset skew]' vectors
%   - P: set of cov(x(n))
%
%%  1. Initialize

N = length(timestamps); % total length of the observations

[dim_x, ~] = size(x1);  % compute the dimension of x1
[dim_P, ~] = size(P1);  % compute the dimension of P1

% check the equality of the dimensions
if dim_x ~= dim_P
    quit;
end

x = zeros(dim_x, N);
for n = 1:N
    P{n} = zeros(dim_x,dim_x);
end

x(:,1) = x1;
P{1} = P1;

var_T = 0;  %   the paper assumes var_T = 0 
var_theta_M = 0.5*(var_C + var_T + var_d);

u = [0,0]'; %   this implementation assumes free running of the slave clock

%%  2. Kalman Filtering  
theta_M_old = 0;
t1_old = 0;
for n = 1:N
  
    t1 = timestamps(n,1);
    t2 = timestamps(n,2);
    t3 = timestamps(n,3);
    t4 = timestamps(n,4);

    if n == 1

        packet_delay = 0.5*( (t4-t1) - (t3-t2) );
        theta_M_old = (t2 - t1) - packet_delay;
        t1_old = t1;

    else

        % 2.1. compute dt
        dt = t1 - t1_old;

        % 2.2. IEEE 1588 PTP offset estimation
        packet_delay = 0.5*( (t4-t1) - (t3-t2) );
        theta_M = (t2 - t1) - packet_delay;         %   by Eq.(16)
        gamma_M = (theta_M - theta_M_old) / dt;     %   by Eq.(17)
        
        % 2.3. set measurement
        z_M = [theta_M; gamma_M];

        % 2.4. set state process matrices
        A = [1, dt;
             0,  1];        %   by Eq.(24)

        B = [-1, -dt;
              0,  -1];      %   by Eq.(25)

        % 2.5. set noise covariance matrices
        Q = diag([var_theta*dt, var_gamma*dt]);
        R = [var_theta_M,         var_theta_M/dt^2;
             var_theta_M/dt^2, 2*(var_theta_M/dt)^2];   %   by Eq.(33)
        
        % 2.6. do the prediction
        x_n_n_1 = A*x(:,n-1) + B*u;     %   by Eq.(27)
        P_n_n_1 = A*P{n-1}*A' + Q;      %   by Eq.(28)

        % 2.7.compute the gain K
        K_n = P_n_n_1*inv( P_n_n_1 + R );            %   by Eq.(30)

        % 2.8. do the correction
        x_n = x_n_n_1 + K_n*(z_M - x_n_n_1);         %   by Eq.(31)
        P_n = ( eye(length(z_M)) - K_n )*P_n_n_1;    %   by Eq.(32)

        % 2.9. finish the filtering for n-step
        x(:,n) = x_n;
          P{n} = P_n; 

        % 2.10. prepare the next step
        theta_M_old = theta_M;
        t1_old = t1;
    end
end

end
