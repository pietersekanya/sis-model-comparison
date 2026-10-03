function [theta_best, ssr_best] = fit_model(t, cum_obs, N, mu, sigma, gamma, n_starts)
%FIT_MODEL  Fit [beta, E0, I0] by nonlinear least squares (lsqnonlin).
%
%   Compares the model's predicted cases in each reporting interval with
%   the observed cases in that interval (no interpolation needed).
%   Runs lsqnonlin from n_starts random starting points and keeps the
%   best fit, to avoid getting stuck in a poor local minimum.
%
%   t        = days since first report (column, t(1) = 0)
%   cum_obs  = observed cumulative cases at those report dates (column)
%   mu       = 0 for basic SEIR, >0 for demographic SEIR

    lb = [0.01,   0,   0];     % lower bounds for [beta, E0, I0]
    ub = [2.00, 1000, 1000];     % upper bounds

    opts = optimoptions('lsqnonlin', 'Display', 'off');

    theta_best = [];
    ssr_best   = inf;

    for k = 1:n_starts
        theta0 = lb + (ub - lb) .* rand(1, 3);
        try
            [th, ssr] = lsqnonlin( ...
                @(th) interval_residuals(th, t, cum_obs, N, mu, sigma, gamma), ...
                theta0, lb, ub, opts);
        catch
            continue            % skip starting points where the solver fails
        end
        if ssr < ssr_best
            ssr_best   = ssr;
            theta_best = th;
        end
    end

    if isempty(theta_best)
        error('fit_model:noFit', 'All starting points failed.');
    end
end


function r = interval_residuals(theta, t, cum_obs, N, mu, sigma, gamma)
% Residuals on the square-root scale (stabilises variance of count data).
    C       = simulate_cum(theta, t, N, mu, sigma, gamma);
    new_mod = diff(C);          % model cases in each reporting interval
    new_obs = diff(cum_obs);    % observed cases in each reporting interval
    r = sqrt(max(new_mod, 0)) - sqrt(new_obs);
    r(~isfinite(r)) = 1e3;      % penalise failed simulations
end
