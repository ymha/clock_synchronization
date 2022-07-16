function [x, X] = ha_lkf(timestamps, gamma_bar, x1, X1, zeta1, var_eta, var_theta_bar, var_gamma_bar)
%%  Ha's LKF
%   This simple source code implements the Ha's Kalman Filtering proposed in
%   Y. Ha, et al., 
%   "Clock Offset Estimation for Systems with Asymmetric Packet Delays,"
%   TechRxiv, 2022. 
%
%%  Assumptions of this example
%   1. Local clock is running without a modification.
%      The assumption implies that u_{x,k} = [0, 0]'.
%   2. For simplicity, assume that observation errors are mutually independent
%      where the observation errors include var_theta_bar and var_gamma_bar.
%   3. Hyper-parameters are given by the INLA or other optimizations.
%      Hyper-parameters consist of {var_theta_bar, var_gamma_bar, var_eta}
%     
%%  State Space Model
%   y_k = C_k*x_k   + D_k*u_{y,k} + v_k,    v_k ~ N(0,R)   by Eq.(18a)
%   x_k = A_k*x_{k-1} + B_k*u_{x_k} + w_k,  w_k ~ N(0,Q)   by Eq.(18b)
%
%   Eq.(18) defines u_{y,k} as [-zeta_k/2; 0].
%   
%   zeta_k = zeta_{k-1} + dzeta_k
%
%   dzeta_k = 2*dt1_k*gamma_k + [1, 1]*alpha_k + [-1, -1]*beta_k   by Eq.(18d)
%   where dti_k = ti_k - ti_{k,1} for i \in {1,2,3,4},
%         alpha_k = [dt1_k; dt4_k], and
%         beta_k  = [dt2_k; dt3_k].
%
%   var_dzeta_k = 4*( (dt1_k)^2*var_gamma_k + (var_eta*(dt1_k)^3)/3 )
%   because var_{w_zeta,k} = (4/3)*var_eta*(dt1_k)^3 by Eq.(22c) for given var_eta.
%
%   Let denote diag([var_theta_bar, var_gamma_bar]) as V.
%   Eq.(19b) and Eq.(23c) define R_k as follow:
%     R_k = diag([var_dzeta_k/4, 0]) + V
%
%   Eq.(22) also defines var_{theta,k} and var_{gamma,k} as:
%     var_{w_theta,k} = ( var_eta * (dt1_k)^3 ) / 3
%     var_{w_gamma,k} = ( var_eta *  dt1_k )
%   
%   Eq.(23b) defines Q_k as follows:
%     Q_k = diag([var_{w_theta,k}, var_{w_gamma,k}])
%     Q_k(1,2) = ( eta_var * dt1_k.^2 ) / 2
%     Q_k(2,1) = Q(1,2)
%
%     S_k(1,1) = Q(1,1)
%     S_k(1,2) = 0
%     S_k(2,1) = Q(2,1)
%     S_k(2,2) = 0
%
%   Eq.(18) defines C_k and D_k as follows:
%     C_k = [1, 0        D_k = [1, 0 
%            0, 1]              0, 0]
%
%   Eq.(21) defines A_k and B_k as follows:
%     A_k = [1, dt1_k    B_k = [-1 -dt1_k
%            0,     1]           0     -1]
% 
%%  Kalman Filtering
%   x_k|y_bar_{k-1} ~ N(x_hat_k, X_hat_k) 
%   x_hat_k = A_k*x_{k-1} + B_k*u_{x,k}
%   X_hat_k = A_k*X_{k-1}*A_k' + Q_k
%
%   y_k|y_bar_{k-1} ~ N(y_hat_k, Y_hat_k)
%   y_hat_k = C_k*x_hat_k + D_k*u_{y,k}
%   Y_hat_k = C_k*X_hat_k*C_k' + R_k
%
%   K_k = X_hat_k*C_k'*inv(Y_hat_k)
%
%   x_k|y_bar_k     ~ N(x_est_k, X_est_k)
%   x_est_k = x_hat_k - K_k*( y_hat_k - y_bar_k )
%   X_est_k = X_hat_k - K_k*C_k*X_hat(t)
%   
%   y_bar_k = [theta_bar; gamma_bar]
%   where theta_bar is given by IEEE 1588 PTP
%
%   x_k = x_est_k
%   X_k = X_est_k
%
%%  If S is not equal to 0 matrix
%   Eq (2.4) and (2.7) in SIEW WAH CHAN et. al. 
%   Convergence Properties of the RDE Optimal Filtering of Nonstabilizable Systems.
%   IEEE Transactions on Automatic Control. 1984
%
%   x_k|y_bar_{k-1}  ~ N(x_hat_k, X_hat_k) 
%   x_hat_k = A_k*x_{k-1} + B_k*u_{x,k} + S_{k-1}*inv(R_{k-1})*err_{k-1} Eq.(2.4)
%   X_hat_k = A_k*X_{t-1}*A_k' + Q_k - S_{k-1}*inv(R_{k-1})*S_{k-1}'     Eq.(2.7)
%
%   err_{k-1} = y_hat_{k-1} - y_bar{k-1}
%
%% Implementation
%
%  Input: timestamps, gamma_bar, x1, X1, var_eta, var_theta_bar, var_gamma_bar
%   - timestamps: IEEE 1588 PTP timestamps [t1 C(t2) C(t3) t4]
%   - gamma_bar: observation of clock skew
%   - x1, X1, zeta1: initial condition of x, X, and zeta
%   - var_eta: the step variance of random walk skew
%   - var_theta_bar: uncertainty of clock offset observation
%   - var_gamma_bar: uncertainty of clock skew observation  
%  
%  Output: x, X 
%   - x: set of [offset skew]' vectors
%   - X: set of cov(x(k))
%
%% 1. Initialization
[N, ~] = size(timestamps);   %   N is the number of observations

x    = zeros(length(x1), N); %   save the history for debugginge slave clock
zeta = zeros(1,N);           %   save the history for debugging

for k = 1:N
    X{k} = zeros(size(X1));
end

V = diag([var_theta_bar, var_gamma_bar]);   %   see assumption 2

%% 2. Linear Kalman Filtering
t1_old = 0;
t2_old = 0;
t3_old = 0;
t4_old = 0;

S_old = zeros(2,2);
R_old = zeros(2,2);

err_old = 0;

for k = 1:N
    
    t1_k = timestamps(k,1);
    t2_k = timestamps(k,2);
    t3_k = timestamps(k,3);
    t4_k = timestamps(k,4);

    if k == 1
        x(:,k) = x1;
        X{k}   = X1;

        zeta(:,k) = zeta1;
        
        S_k = zeros(2,2);
        Q_k = zeros(size(X));
        R_k = zeros(2,2);

        err_k = zeros(size(x1));
    else        
        % 2.1. explore observation 
        packet_delay = 0.5 * ( (t4_k - t1_k) - (t3_k - t2_k) );
        ptp_offset = (t2_k - t1_k) - packet_delay;

        theta_bar = ptp_offset;
        y_bar_k = [theta_bar; gamma_bar(k)];

        % 2.2. predict zeta_k|y_bar_{k-1}
        dt1_k = t1_k - t1_old;
        dt2_k = t2_k - t2_old;
        dt3_k = t3_k - t3_old;
        dt4_k = t4_k - t4_old;

        alpha_k = [dt1_k; dt4_k];
        beta_k  = [dt2_k; dt3_k];

        gamma_k = x(2,k-1);
        var_gamma_k = X{i-1}(2,2);
        
        % by Eq.(18d)
        dzeta_k_k_1 = 2*dt1_k*gamma_k + [1, 1]*alpha_k + [-1, -1]*beta_k;
        var_dzeta_k = 4*( (dt1_k)^2*var_gamma_k + (var_eta*(dt1_k)^3)/3 );
        
        % predict zeta_k
        zeta_k_pred = zeta(k-1) + dzeta_k_k_1;        

        % 2.3. compute u_{y,k}
        u_y_k = [-zeta_k_pred/2; 0];
        
        u_x_k = [0; 0];     % by assumption 1
        V_k = V;            %   we do not update V_k in this example
                            %   role of INLA or sequential INLA

        % 2.4. compute R,Q,S
        % by Eq.(19b) and Eq.(23c)
        R_k = diag([var_dzeta_k/4, 0]) + V_k;

        % by Eq.(22)
        var_w_theta_k = ( var_eta * (dt1_k)^3 ) / 3;
        var_w_gamma_k = ( var_eta *  dt1_k );

        % by Eq.(23b)
        Q_k = diag( [var_w_theta_k, var_w_gamma_k] );
        Q_k(1,2) = ( eta_var * dt1_k.^2 ) / 2;
        Q_k(2,1) = Q(1,2);

        S_k(1,1) = Q(1,1);
        S_k(1,2) = 0;
        S_k(2,1) = Q(2,1);
        S_k(2,2) = 0;

        % 2.5. compute A,B,C,D
        % by Eq.(18)
        C_k = [1, 0; 0, 1];
        D_k = [1, 0; 0, 0];
        
        % by Eq.(21)
        A_k = [1,  dt1_k; 0, 1];
        B_k = [-1 -dt1_k; 0, -1];  

        % 2.6. do predict
        x_hat_k = A_k*x(:,k-1) + B_k*u_x_k + S_old*inv(R_old)*err_old;
        X_hat_k = A_k*X{k-1}*A_k' + (Q_k-S_old*inv(R_old)*S_old');

        % 2.7. do innovate
        y_hat_k = C_k*x_hat_k + D_k*u_y_k;
        Y_hat_k = C_k*X_hat_k*C_k' + R_k;

        % 2.8. do Kalman filtering
        K_k = X_hat_k * C_k' * inv(Y_hat_k);

        err_k = y_hat_k - y_bar_k;

        x_est_k = x_hat_k - K_k * err_k;
        X_est_k = X_hat_k - K_k*C_k*X_hat_k;

        x(:,k) = x_est_k;
        X{k} = X_est_k;

        % 2.9. do apply the filtering to zeta_k
        gamma_k = x(2,k);
        dzeta_k_k = 2*dt1_k*gamma_k + [1, 1]*alpha_k + [-1, -1]*beta_k;
        
        zeta(1,k) = zeta(1,k-1) + dzeta_k_k;
    end
    
    % prepare the next step
    t1_old = t1_k;
    t2_old = t2_k;
    t3_old = t3_k;
    t4_old = t4_k;

    S_old = S_k;
    R_old = R_k;

    err_old = err_k;
end

end
