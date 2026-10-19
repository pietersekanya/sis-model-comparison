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
%                (used for approximate standard errors)

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
