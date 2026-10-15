function [C, D] = simulate_sis_cum(p, t_query, N, b, mu)
%SIMULATE_SIS_CUM  Cumulative cases and deaths from the demographic SIS model.
%
%   Model:
%       dS/dt = b*N - beta*S*I + gamma*I - mu*S
%       dI/dt = beta*S*I - gamma*I - mu_tilde*I
%       dC/dt = beta*S*I                (C = cumulative new infections = cases)
%       dD/dt = mu_tilde*I              (D = cumulative disease-induced deaths)
%
%   p        = [betaN, gamma, mu_tilde, I0]
%                betaN    = beta*N, transmission rate (per day), so that
%                           beta = betaN/N (per person per day)
%                gamma    = recovery rate (per day)
%                mu_tilde = disease-induced death rate (per day)
%                I0       = infectious individuals at t = 0
%   t_query  = times in days since the start of the fitting window
%              (column vector, first value 0, values increasing)
%   N        = population at the start of the window
%   b, mu    = birth and natural death rates (per person per day)
%
%   Returns C(t) and D(t) at the requested times, with C(0) = D(0) = 0.
%   Cases and deaths in any reporting interval are given by diff(C) and
%   diff(D). Callers that only need the cases can use  C = simulate_sis_cum(...).

beta     = p(1) / N;
gamma    = p(2);
mu_tilde = p(3);
I0       = p(4);

rhs = @(~, x) [ ...
    b*N - beta*x(1)*x(2) + gamma*x(2) - mu*x(1); ...   % dS/dt
    beta*x(1)*x(2) - gamma*x(2) - mu_tilde*x(2); ...   % dI/dt
    beta*x(1)*x(2); ...                                % dC/dt
    mu_tilde*x(2) ];                                   % dD/dt

x0   = [N - I0; I0; 0; 0];                             % [S; I; C; D]
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:4);
sol  = ode45(rhs, [0, max(t_query)], x0, opts);
X    = deval(sol, t_query(:)')';      % one row per requested time

C = X(:, 3);
D = X(:, 4);
end