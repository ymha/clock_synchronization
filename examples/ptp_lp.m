function [gamma, theta, x_ub, x_lb] = ptp_lp(timestamps)

%%  PTP-LP
%   This simple source code implements the PTP-LP proposed in
%   H. Puttnies, P. Danielis, and D. Timmermann, 
%   “PTP-LP: Using Linear Programming to Increase the Delay Robustness of IEEE 1588 PTP,” 
%   in 2018 IEEE Global Communications Conference, IEEE, 2018.

%% Implementation
%
%  Input: timestamps
%   - timestamps: IEEE 1588 PTP timestamps [t1 C(t2) C(t3) t4]
%  
%  Output: gamma, theta, x_ub, x_lb
%   - gamma: clock skew
%   - theta: clock offset
%   - x_ub : information related to the upper bound of the slave clock
%   - x_lb : information related to the lower bound of the slave clock
%

% set the number of timestamps
[N,~] = size(timestamps);

% get set of t1, t2, t3, t4
t1 = timestamps(:,1);
t2 = timestamps(:,2);	%	t2 == C(t2)
t3 = timestamps(:,3);	%	t3 == C(t3)
t4 = timestamps(:,4);

%% upper bound 
% Eq.(11) alpha1.*T1_n + beta1_n <= C(T2_n)
A_ub(1:N,1) = t1;
A_ub(1:N,2) = 1;

b_ub(1,1:N) = t2';

% Eq.(12) minimization of f(alpha1,beta1) = Sum_n^N {C(T2)_n - alpha1.*T1_n - beta1}
% Because Sum_n^N C(T2)_n is deterministic for given C(T2)_n
% we try to minimize the terms 
% Sum_n^N {-alpha1.*T1_n - beta1} = -(Sum_n^N T1_n).*alpha1 - N.*beta1 
f_ub = [-sum(t1), -N];

% lp-solver
x_ub = linprog(f_ub, A_ub, b_ub);

%% lower bound
% Eq.(15) alpha2.*T4_n + beta2_n >= C(T3_n)
A_lb(1:N,1) = -t4;
A_lb(1:N,2) = -1;

b_lb(1,1:N) = -t3';

% Eq.(12) minimization of f(alpha2,beta2) = Sum_n^N {alpha2.*T4_n + beta2 - C(T3)_n}
% Because Sum_n^N C(T3)_n is deterministic
% we try to minimize the terms 
% Sum_n^N {alpha2.*T4_n + beta2} = (Sum_n^N T4_n).*alpha2 + N.*beta2
f_lb = [sum(t4), N];

% lp-solver
x_lb = linprog(f_lb, A_lb, b_lb);

%% Eq.(17)
gamma = (x_ub(1,1) + x_lb(1,1)) * 0.5;
theta = (x_ub(2,1) + x_lb(2,1)) * 0.5;

end
