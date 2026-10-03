function C = simulate_cum(theta, t_query, N, mu, sigma, gamma)
%SIMULATE_CUM  Cumulative new cases predicted by an SEIR model.
%
%   theta   = [beta, E0, I0]   parameters being fitted
%   t_query = times in days since the first report (column vector,
%             first value must be 0, values increasing)
%   N       = population at risk
%   mu      = birth/death rate per day
%               mu = 0    -> basic SEIR
%               mu > 0    -> demographic SEIR (births = deaths = mu)
%   sigma   = 1 / incubation period (per day)
%   gamma   = 1 / infectious period (per day)
%
%   Returns C(t): cumulative new cases (the E -> I flow) since t = 0.
%   Cases in any interval are then given by diff(C).

    beta = theta(1);
    E0   = theta(2);
    I0   = theta(3);

    x0 = [N - E0 - I0; E0; I0; 0; 0];     % [S; E; I; R; C]

    rhs = @(~, x) [ ...
        mu*N - beta*x(1)*x(3)/N - mu*x(1); ...       % dS/dt
        beta*x(1)*x(3)/N - (sigma + mu)*x(2); ...    % dE/dt
        sigma*x(2) - (gamma + mu)*x(3); ...          % dI/dt
        gamma*x(3) - mu*x(4); ...                    % dR/dt
        sigma*x(2) ];                                % dC/dt (cumulative cases)

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:5);
    sol  = ode45(rhs, [0, max(t_query)], x0, opts);
    X    = deval(sol, t_query(:)')';      % one row per requested time
    C    = X(:, 5);
end
