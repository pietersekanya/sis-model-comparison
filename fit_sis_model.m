function [p_best, ssr_best, res_best, J_best] = fit_sis_model( ...
    t, cum_obs, deaths_obs, N, b, mu, p_fixed, free, lb, ub, n_starts)
%FIT_SIS_MODEL  Fit the demographic SIS model to cases and deaths.
%
%   Parameter vector p = [betaN, gamma, mu_tilde, I0] (see simulate_sis_cum).
%
%   t           = days since the start of the fitting window (column, t(1) = 0)
%   cum_obs     = observed cumulative cases on those report dates (column)
%   deaths_obs  = observed cumulative deaths on the same dates (column),
%                 or [] to fit the cases only
%   free        = 1x4 logical: which parameters are fitted
%   p_fixed     = 1x4 values used for the parameters that are NOT fitted
%                 (entries for fitted parameters are ignored)
%   lb, ub      = 1x4 bounds (only the fitted parameters are used)
%
%   The model is fitted to the observed cases AND deaths in each actual
%   reporting interval (no interpolation), on the square-root scale. The
%   case residuals and the death residuals are stacked into one vector, so
%   each observed interval counts equally. lsqnonlin is run from n_starts
%   random starting points and the best fit is kept.
%
%   Outputs:
%     p_best   = 1x4 full parameter vector
%     ssr_best = sum of squared residuals
%     res_best = residual vector at the optimum (cases first, then deaths)
%     J_best   = Jacobian with respect to the FITTED parameters only

idx  = find(free);
lb_f = lb(idx);
ub_f = ub(idx);

opts = optimoptions('lsqnonlin', 'Display', 'off');

x_best   = [];
ssr_best = inf;
res_best = [];
J_best   = [];

for k = 1:n_starts
    x0 = lb_f + (ub_f - lb_f) .* rand(1, numel(idx));
    try
        [x, ssr, res, ~, ~, ~, J] = lsqnonlin( ...
            @(x) interval_residuals(x, idx, p_fixed, t, cum_obs, deaths_obs, N, b, mu), ...
            x0, lb_f, ub_f, opts);
    catch
        continue            % skip starting points where the solver fails
    end
    if ssr < ssr_best
        ssr_best = ssr;
        x_best   = x;
        res_best = res;
        J_best   = J;
    end
end

if isempty(x_best)
    error('fit_sis_model:noFit', 'All starting points failed.');
end

p_best = p_fixed;
p_best(idx) = x_best;


% Warn if a fitted parameter sits on one of its bounds
tol = 1e-3 * (ub_f - lb_f);
near_lower = abs(x_best - lb_f) < tol;
near_upper = abs(x_best - ub_f) < tol;
at_bound = near_lower | near_upper;

if any(at_bound)
    names4 = {'betaN', 'gamma', 'mu_tilde', 'I0'};
    names_fitted = names4(idx);

    for j = find(at_bound)
        if near_lower(j)
            bound_type = 'LOWER';
            bound_value = lb_f(j);
        else
            bound_type = 'UPPER';
            bound_value = ub_f(j);
        end

        warning('fit_sis_model:atBound', ...
            ['Parameter %s is near its %s bound. ' ...
            'Estimate = %.6g; bound = %.6g.'], ...
            names_fitted{j}, bound_type, ...
            x_best(j), bound_value);
    end
end
end


function r = interval_residuals(x, idx, p_fixed, t, cum_obs, deaths_obs, N, b, mu)
% Residuals on the square-root scale, one per reporting interval:
% the case residuals first, then (if deaths are supplied) the death residuals.
p = p_fixed;
p(idx) = x;
[C, D] = simulate_sis_cum(p, t, N, b, mu);

r = sqrt(max(diff(C), 0)) - sqrt(max(diff(cum_obs), 0));      % cases

if ~isempty(deaths_obs)
    r_d = sqrt(max(diff(D), 0)) - sqrt(max(diff(deaths_obs), 0));   % deaths
    r = [r; r_d];
end

r(~isfinite(r)) = 1e3;                   % penalise failed simulations
end