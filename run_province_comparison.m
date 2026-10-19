%% run_province_comparison.m
% Fit a demographic SIS model separately to Ituri and North Kivu.
% Includes parameter estimation with approximate standard errors, a
% comparison of the parameters between the provinces (approximate z-tests),
% identifiability diagnostics, sensitivity analysis and
% leave-future-out cross-validation.
%
% Required files:
% load_province_cases.m
% load_province_deaths.m
% province_settings.m
% fit_sis_model.m
% simulate_sis_cum.m
%
% Run from the top-level project folder.

clear; clc; close all;
rng(1);

%% 1. MODEL SETTINGS

% Parameter vector: p = [betaN, gamma, mu_tilde, I0]
% betaN    = beta*N, units 1/day
% gamma    = recovery rate, units 1/day
% mu_tilde = disease-related mortality rate, units 1/day
% I0       = initial infected individuals, units people

names4 = {'betaN', 'gamma', 'mu_tilde', 'I0'};

free = [true, false, true, true];
p_fixed = [NaN, 1/10, NaN, NaN];

% Parameter bounds
lb = [0.01, 1/30, 0.0, 0.1];
ub = [5.00, 1/2,  0.5, 2000];

n_starts = 10;

% Sensitivity analysis: assumed mean recovery periods
infectious_periods = [7, 10, 14];
gamma_values = 1 ./ infectious_periods;

% Leave-future-out cross-validation settings
first_origin = 40;
origin_step = 5;
h = 5;
min_train = 15;

%% 2. LOAD DATA

[P, ~] = load_province_cases();
[Dth, ~] = load_province_deaths();
S = province_settings();

if ~isequal(P.Date, Dth.Date)
    error('Case and death dates do not match.');
end

report_date = P.Date;
t = days(report_date - report_date(1));
n = numel(t);
np = numel(S);

if np ~= 2
    error('The comparison between provinces assumes exactly two provinces.');
end

%% 3. FIT THE MAIN MODEL TO EACH PROVINCE

% Results include parameter estimates, approximate standard errors,
% derived quantities, and in-sample MAE, RMSE and WAPE for cases and deaths.

results = table();

% Derived quantities as functions of the full parameter vector
Rfun   = @(pp) pp(1) / (pp(2) + pp(3));      % initial reproduction measure
CFRfun = @(pp) pp(3) / (pp(2) + pp(3));      % model-implied case fatality proportion

% Estimates and approximate standard errors.
% Columns: betaN, gamma, mu_tilde, I0, R_start, CFR
est = nan(np, 6);
se  = zeros(np, 6);
cond_main = zeros(np, 1);
corr_main = cell(np, 1);

% Save each province's fitting window for sensitivity analysis and CV.
start_idx = zeros(np, 1);
tw_all = cell(np, 1);
cw_all = cell(np, 1);
dw_all = cell(np, 1);

% Store in-sample observations and predictions for plotting.
fit_dates = cell(np, 1);
fit_obs_c = cell(np, 1);
fit_pred_c = cell(np, 1);
fit_obs_d = cell(np, 1);
fit_pred_d = cell(np, 1);

% Largest reporting intervals, for the discussion of unusually
% large counts (one table per province and series).
largest = cell(np, 2);

fig_fit = figure('Color', 'w', ...
    'Position', [100 100 1000 650]);

for q = 1:np

    cum_cases = P.(S(q).column);
    cum_deaths = Dth.(S(q).column);

    % Begin at the report before the first positive cumulative case,
    % where possible.
    first_pos = find(cum_cases > 0, 1);

    if isempty(first_pos)
        error('%s has no confirmed cases.', S(q).name);
    end

    s = max(first_pos - 1, 1);
    start_idx(q) = s;

    % Use actual reporting dates; do not interpolate.
    tw = t(s:end) - t(s);
    cases = cum_cases(s:end);
    deaths = cum_deaths(s:end);

    if any(diff(cases) < 0) || any(diff(deaths) < 0)
        error('Cumulative counts decrease for %s.', S(q).name);
    end

    tw_all{q} = tw;
    cw_all{q} = cases;
    dw_all{q} = deaths;

    fprintf('\n=== MAIN FIT: %s ===\n', S(q).name);

    % Estimate betaN, mu_tilde and I0; gamma is fixed at 1/10.
    [p_hat, ~, res, J] = fit_sis_model( ...
        tw, cases, deaths, S(q).N, S(q).b, S(q).mu, ...
        p_fixed, free, lb, ub, n_starts);

    betaN = p_hat(1);
    gamma = p_hat(2);
    mu_tilde = p_hat(3);
    I0 = p_hat(4);

    beta = betaN / S(q).N;
    R_start = Rfun(p_hat);
    CFR = CFRfun(p_hat);

    % Approximate standard errors from the Jacobian at the optimum
    % (Gauss-Newton approximation), and by the delta method for the
    % derived quantities.
    [C_cov, R_corr, cond_JTJ] = parameter_covariance(res, J);

    se4 = zeros(1, 4);
    se4(free) = sqrt(max(diag(C_cov), 0))';

    se_R = delta_se(Rfun, p_hat, free, C_cov);
    se_CFR = delta_se(CFRfun, p_hat, free, C_cov);

    est(q, :) = [p_hat, R_start, CFR];
    se(q, :)  = [se4, se_R, se_CFR];
    cond_main(q) = cond_JTJ;
    corr_main{q} = R_corr;

    % Simulate cumulative cases and deaths.
    [C_fit, D_fit] = simulate_sis_cum( ...
        p_hat, tw, S(q).N, S(q).b, S(q).mu);

    % Convert cumulative counts into reporting-interval counts.
    obs_cases = diff(cases);
    pred_cases = diff(C_fit);
    obs_deaths = diff(deaths);
    pred_deaths = diff(D_fit);

    % In-sample error metrics.
    MAE_cases = mean(abs(pred_cases - obs_cases));
    RMSE_cases = sqrt(mean((pred_cases - obs_cases).^2));

    MAE_deaths = mean(abs(pred_deaths - obs_deaths));
    RMSE_deaths = sqrt(mean((pred_deaths - obs_deaths).^2));

    % Weighted absolute percentage error (scale-free):
    % sum of absolute errors divided by the sum of observed counts.
    WAPE_cases = sum(abs(pred_cases - obs_cases)) / sum(obs_cases);
    WAPE_deaths = sum(abs(pred_deaths - obs_deaths)) / max(sum(obs_deaths), eps);

    % Print estimates and fit metrics.
    fprintf('betaN       = %.5g per day (SE %.3g)\n', betaN, se4(1));
    fprintf('beta        = %.5g per person per day\n', beta);
    fprintf('gamma       = %.5g per day (fixed)\n', gamma);
    fprintf('mu_tilde    = %.5g per day (SE %.3g)\n', mu_tilde, se4(3));
    fprintf('I0          = %.5g people (SE %.3g)\n', I0, se4(4));
    fprintf('R_start     = %.5g (SE %.3g)\n', R_start, se_R);
    fprintf('Model CFR   = %.2f%% (SE %.2f percentage points)\n', 100*CFR, 100*se_CFR);

    fprintf('Cases:  MAE = %.3f, RMSE = %.3f, WAPE = %.3f\n', ...
        MAE_cases, RMSE_cases, WAPE_cases);
    fprintf('Deaths: MAE = %.3f, RMSE = %.3f, WAPE = %.3f\n', ...
        MAE_deaths, RMSE_deaths, WAPE_deaths);

    % Identifiability check for the main fit.
    fprintf('Identifiability: cond(J''J) = %.2e\n', cond_JTJ);
    if numel(R_corr) > 1
        disp(array2table(R_corr, 'VariableNames', names4(free), ...
            'RowNames', names4(free)));
    end

    % Store results.
    new_row = table(string(S(q).name), betaN, beta, gamma, ...
        mu_tilde, I0, R_start, CFR, ...
        MAE_cases, RMSE_cases, WAPE_cases, ...
        MAE_deaths, RMSE_deaths, WAPE_deaths, ...
        'VariableNames', {'Province', 'betaN', 'beta', ...
        'gamma', 'mu_tilde', 'I0', 'R_start', 'CFR', ...
        'MAE_cases', 'RMSE_cases', 'WAPE_cases', ...
        'MAE_deaths', 'RMSE_deaths', 'WAPE_deaths'});

    results = [results; new_row];

    % Save plotting data.
    fit_dates{q} = report_date(s+1:end);
    fit_obs_c{q} = obs_cases;
    fit_pred_c{q} = pred_cases;
    fit_obs_d{q} = obs_deaths;
    fit_pred_d{q} = pred_deaths;

    % Largest reporting intervals (cases and deaths per day,
    % together with the number of days each interval covers).
    dt_int = diff(tw);                      % days in each reporting interval
    rate_c = obs_cases ./ dt_int;           % cases per day in each interval
    rate_d = obs_deaths ./ dt_int;          % deaths per day in each interval
    z_c = (rate_c - mean(rate_c)) ./ std(rate_c);
    z_d = (rate_d - mean(rate_d)) ./ std(rate_d);

    largest{q, 1} = top_intervals(S(q).name, 'Cases', ...
        fit_dates{q}, obs_cases, dt_int, z_c);
    largest{q, 2} = top_intervals(S(q).name, 'Deaths', ...
        fit_dates{q}, obs_deaths, dt_int, z_d);

    fprintf(['Reporting intervals (days): min = %g, median = %g, max = %g. ' ...
        'Intervals more than 2 SD above the mean: %d (cases per day), %d (deaths per day)\n'], ...
        min(dt_int), median(dt_int), max(dt_int), sum(z_c > 2), sum(z_d > 2));

    % Plot cases.
    subplot(2, np, q);
    plot(fit_dates{q}, obs_cases, 'ko', ...
        'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on;
    plot(fit_dates{q}, pred_cases, '-', 'LineWidth', 1.5);
    xlabel('Date');
    ylabel('Cases per reporting interval');
    title([S(q).name ': cases']);
    legend({'Observed', 'Fitted SIS'}, 'Location', 'best');
    grid on;

    % Plot deaths.
    subplot(2, np, np+q);
    plot(fit_dates{q}, obs_deaths, 'ko', ...
        'MarkerFaceColor', [0.7 0.7 0.7]);
    hold on;
    plot(fit_dates{q}, pred_deaths, '-', 'LineWidth', 1.5);
    xlabel('Date');
    ylabel('Deaths per reporting interval');
    title([S(q).name ': deaths']);
    legend({'Observed', 'Fitted SIS'}, 'Location', 'best');
    grid on;

end

%% 3b. LARGEST REPORTING INTERVALS

% The five largest case counts and the five largest death counts in a
% single reporting interval, for each province, with the length of the
% interval behind each count. Long intervals suggest that delayed reports
% were batched together; 1-day intervals suggest a genuine surge or a
% correction to earlier counts.

largest_all = vertcat(largest{:});

fprintf('\n=== LARGEST REPORTING INTERVALS ===\n');
disp(largest_all);

%% 3c. COMPARISON OF THE PROVINCES (APPROXIMATE Z-TESTS)

% Difference between the provinces (first minus second) for each parameter,
% using the approximate standard errors above and assuming that the two
% estimates are independent. These tests are approximate: the residuals
% from cumulative data are autocorrelated and the model does not capture
% every feature of the data, so the standard errors are probably too small.

col_names = {'betaN', 'gamma', 'mu_tilde', 'I0', 'R_start', 'CFR'};

d_est = est(1, :) - est(2, :);
d_se  = sqrt(se(1, :).^2 + se(2, :).^2);
d_lo  = d_est - 1.96 * d_se;
d_hi  = d_est + 1.96 * d_se;
z     = d_est ./ d_se;
z(d_se == 0) = NaN;                       % parameters that are fixed
pval  = erfc(abs(z) / sqrt(2));           % two-sided normal p-value

comparison = table(col_names', est(1, :)', est(2, :)', ...
    d_est', d_se', d_lo', d_hi', z', pval', ...
    'VariableNames', {'Parameter', ...
    [S(1).column '_Estimate'], [S(2).column '_Estimate'], ...
    'Difference', 'SE_Difference', 'CI_Lower', 'CI_Upper', 'z', 'p_value'});

fprintf('\n=== COMPARISON OF PARAMETERS (%s minus %s) ===\n', ...
    S(1).name, S(2).name);
disp(comparison);

% Long-format table of estimates with approximate 95% intervals.
Prov = {};  Par = {};  Est = [];  SEv = [];
for q = 1:np
    for j = 1:numel(col_names)
        Prov{end+1, 1} = S(q).name;             %#ok<SAGROW>
        Par{end+1, 1}  = col_names{j};          %#ok<SAGROW>
        Est(end+1, 1)  = est(q, j);             %#ok<SAGROW>
        SEv(end+1, 1)  = se(q, j);              %#ok<SAGROW>
    end
end
params_ci = table(Prov, Par, Est, SEv, Est - 1.96*SEv, Est + 1.96*SEv, ...
    'VariableNames', {'Province', 'Parameter', 'Estimate', 'SE', ...
    'CI_Lower', 'CI_Upper'});

% Plot the estimates with approximate 95% intervals.
fig_par = figure('Color', 'w', 'Position', [100 100 900 650]);
plot_cols = [1 3 5 6];                 % betaN, mu_tilde, R_start, CFR
for j = 1:4
    subplot(2, 2, j);
    col = plot_cols(j);
    errorbar(1:np, est(:, col), 1.96 * se(:, col), 'o', 'LineWidth', 1.5);
    xlim([0.5, np + 0.5]);
    xticks(1:np);
    xticklabels({S.name});
    title(col_names{col});
    grid on;
end

%% 4. SENSITIVITY ANALYSIS FOR GAMMA

% Refit both provinces using assumed mean recovery periods
% of 7, 10 and 14 days.

sens_rows = [];

fprintf('\n=== SENSITIVITY ANALYSIS ===\n');

for g = 1:numel(gamma_values)

    for q = 1:np

        pf = p_fixed;
        pf(2) = gamma_values(g);

        ph = fit_sis_model( ...
            tw_all{q}, cw_all{q}, dw_all{q}, ...
            S(q).N, S(q).b, S(q).mu, ...
            pf, free, lb, ub, n_starts);

        betaN_s = ph(1);
        gamma_s = ph(2);
        mu_s = ph(3);
        I0_s = ph(4);

        R_s = betaN_s / (gamma_s + mu_s);
        CFR_s = mu_s / (gamma_s + mu_s);

        sens_rows(end+1, :) = [ ...
            gamma_s, 1/gamma_s, q, betaN_s, ...
            mu_s, I0_s, R_s, CFR_s]; %#ok<SAGROW>

    end
end

sens_table = table( ...
    sens_rows(:,1), sens_rows(:,2), ...
    string({S(sens_rows(:,3)).name})', ...
    sens_rows(:,4), sens_rows(:,5), sens_rows(:,6), ...
    sens_rows(:,7), sens_rows(:,8), ...
    'VariableNames', {'Gamma', 'MeanRecoveryDays', ...
    'Province', 'betaN', 'mu_tilde', 'I0', ...
    'R_start', 'CFR'});

disp(sens_table);

%% 4b. IDENTIFIABILITY DIAGNOSTICS

% Why gamma is fixed: refit the CASES ONLY with betaN, gamma, mu_tilde and
% I0 all estimated, and compare with the main fit (cases and deaths,
% gamma fixed). A very large condition number and correlations close to +1
% or -1 mean that the parameters cannot be estimated separately.

free_all = true(1, 4);
p_fixed_all = [NaN, NaN, NaN, NaN];

diag_rows = {};
ws = warning('off', 'all');            % bound warnings are expected here

fprintf('\n=== IDENTIFIABILITY DIAGNOSTICS ===\n');

for q = 1:np
    [p_c, ~, res_c, J_c] = fit_sis_model( ...
        tw_all{q}, cw_all{q}, [], S(q).N, S(q).b, S(q).mu, ...
        p_fixed_all, free_all, lb, ub, n_starts);
    [~, R_c, cond_c] = parameter_covariance(res_c, J_c);
    max_corr_c = max(abs(R_c(~eye(4))));

    at_bound = abs(p_c - lb) < 1e-3 * (ub - lb) | abs(p_c - ub) < 1e-3 * (ub - lb);

    R_m = corr_main{q};
    max_corr_m = max(abs(R_m(~eye(size(R_m)))));

    diag_rows(end+1, :) = {S(q).name, ...
        'Cases only; betaN, gamma, mu_tilde, I0 estimated', ...
        cond_c, max_corr_c, strjoin(names4(at_bound), ', ')}; %#ok<SAGROW>
    diag_rows(end+1, :) = {S(q).name, ...
        'Cases and deaths; gamma fixed; betaN, mu_tilde, I0 estimated', ...
        cond_main(q), max_corr_m, 'none'}; %#ok<SAGROW>
end

warning(ws);

diag_table = cell2table(diag_rows, 'VariableNames', ...
    {'Province', 'Fit', 'ConditionNumber', 'LargestAbsCorrelation', 'ParametersAtBounds'});
disp(diag_table);

%% 5. LEAVE-FUTURE-OUT CROSS-VALIDATION

% At each forecast origin, fit only to data available up to that date.
% Forecast the next h reporting intervals.
% Forecast windows use actual report dates.

origins = first_origin:origin_step:(n-h);

if isempty(origins)
    error('Not enough data for the selected CV settings.');
end

n_origins = numel(origins);

mae_cv_c = nan(n_origins, np);
rmse_cv_c = nan(n_origins, np);
wape_cv_c = nan(n_origins, np);
mae_cv_d = nan(n_origins, np);
rmse_cv_d = nan(n_origins, np);
wape_cv_d = nan(n_origins, np);

fprintf('\n=== LEAVE-FUTURE-OUT CROSS-VALIDATION ===\n');

for q = 1:np

    cum_cases = P.(S(q).column);
    cum_deaths = Dth.(S(q).column);
    s = start_idx(q);

    for k = 1:n_origins

        c = origins(k);

        % Ensure sufficient training data for this province.
        if c - s + 1 < min_train
            continue;
        end

        % Fit using data up to and including report c.
        p_cv = fit_sis_model( ...
            t(s:c) - t(s), ...
            cum_cases(s:c), cum_deaths(s:c), ...
            S(q).N, S(q).b, S(q).mu, ...
            p_fixed, free, lb, ub, n_starts);

        % Simulate through the next h reporting intervals.
        [C_cv, D_cv] = simulate_sis_cum( ...
            p_cv, t(s:c+h) - t(s), ...
            S(q).N, S(q).b, S(q).mu);

        % Predicted and observed interval counts for c+1 ... c+h.
        pred_c = diff(C_cv(end-h:end));
        obs_c = diff(cum_cases(c:c+h));

        pred_d = diff(D_cv(end-h:end));
        obs_d = diff(cum_deaths(c:c+h));

        err_c = pred_c - obs_c;
        err_d = pred_d - obs_d;

        mae_cv_c(k,q) = mean(abs(err_c));
        rmse_cv_c(k,q) = sqrt(mean(err_c.^2));
        wape_cv_c(k,q) = sum(abs(err_c)) / max(sum(obs_c), eps);

        mae_cv_d(k,q) = mean(abs(err_d));
        rmse_cv_d(k,q) = sqrt(mean(err_d.^2));
        if sum(obs_d) > 0
            wape_cv_d(k,q) = sum(abs(err_d)) / sum(obs_d);
        end

    end

    fprintf('Finished CV for %s.\n', S(q).name);
end

% Average performance across valid forecast origins.
cv_summary = table( ...
    string({S.name})', ...
    mean(mae_cv_c, 1, 'omitnan')', ...
    mean(rmse_cv_c, 1, 'omitnan')', ...
    mean(wape_cv_c, 1, 'omitnan')', ...
    mean(mae_cv_d, 1, 'omitnan')', ...
    mean(rmse_cv_d, 1, 'omitnan')', ...
    mean(wape_cv_d, 1, 'omitnan')', ...
    sum(~isnan(mae_cv_c), 1)', ...
    'VariableNames', {'Province', ...
    'CV_MAE_cases', 'CV_RMSE_cases', 'CV_WAPE_cases', ...
    'CV_MAE_deaths', 'CV_RMSE_deaths', 'CV_WAPE_deaths', ...
    'Origins_used'});

disp(cv_summary);

% Save errors for each forecast origin.
cv_by_origin = table( ...
    report_date(origins), ...
    mae_cv_c(:,1), mae_cv_c(:,2), ...
    rmse_cv_c(:,1), rmse_cv_c(:,2), ...
    wape_cv_c(:,1), wape_cv_c(:,2), ...
    mae_cv_d(:,1), mae_cv_d(:,2), ...
    rmse_cv_d(:,1), rmse_cv_d(:,2), ...
    wape_cv_d(:,1), wape_cv_d(:,2), ...
    'VariableNames', {'ForecastOrigin', ...
    'MAE_cases_Ituri', 'MAE_cases_NorthKivu', ...
    'RMSE_cases_Ituri', 'RMSE_cases_NorthKivu', ...
    'WAPE_cases_Ituri', 'WAPE_cases_NorthKivu', ...
    'MAE_deaths_Ituri', 'MAE_deaths_NorthKivu', ...
    'RMSE_deaths_Ituri', 'RMSE_deaths_NorthKivu', ...
    'WAPE_deaths_Ituri', 'WAPE_deaths_NorthKivu'});

%% 6. PLOT CROSS-VALIDATION PERFORMANCE

fig_cv = figure('Color', 'w', ...
    'Position', [100 100 1000 400]);

subplot(1,2,1);
plot(report_date(origins), rmse_cv_c(:,1), '-o', ...
    'LineWidth', 1.2);
hold on;
plot(report_date(origins), rmse_cv_c(:,2), '-s', ...
    'LineWidth', 1.2);
xlabel('Forecast origin');
ylabel('RMSE (cases)');
title('Future case prediction error');
legend({S.name}, 'Location', 'best');
grid on;

subplot(1,2,2);
plot(report_date(origins), rmse_cv_d(:,1), '-o', ...
    'LineWidth', 1.2);
hold on;
plot(report_date(origins), rmse_cv_d(:,2), '-s', ...
    'LineWidth', 1.2);
xlabel('Forecast origin');
ylabel('RMSE (deaths)');
title('Future death prediction error');
legend({S.name}, 'Location', 'best');
grid on;

%% 7. SAVE RESULTS

if ~isfolder('results')
    mkdir('results');
end

writetable(results, fullfile('results', ...
    'sis_province_parameters.csv'));

writetable(params_ci, fullfile('results', ...
    'sis_province_parameters_with_se.csv'));

writetable(comparison, fullfile('results', ...
    'sis_province_comparison.csv'));

writetable(sens_table, fullfile('results', ...
    'sis_province_gamma_sensitivity.csv'));

writetable(diag_table, fullfile('results', ...
    'sis_province_identifiability.csv'));

writetable(cv_summary, fullfile('results', ...
    'sis_province_cv_summary.csv'));

writetable(cv_by_origin, fullfile('results', ...
    'sis_province_cv_by_origin.csv'));

writetable(largest_all, fullfile('results', ...
    'sis_province_largest_intervals.csv'));

saveas(fig_fit, fullfile('results', ...
    'sis_province_fits.png'));

saveas(fig_par, fullfile('results', ...
    'sis_province_parameters.png'));

saveas(fig_cv, fullfile('results', ...
    'sis_province_cv_error.png'));

fprintf('\nAll results saved in the results/ folder.\n');

%% LOCAL FUNCTIONS (must stay at the end of the script)

function [C, R, kappa] = parameter_covariance(res, J)
% Approximate covariance, correlation and conditioning of the fitted
% parameters from the Jacobian at the optimum (Gauss-Newton approximation).
    J   = full(J);
    n   = numel(res);
    k   = size(J, 2);
    s2  = sum(res.^2) / max(n - k, 1);
    JTJ = J' * J;
    C   = s2 * pinv(JTJ);
    d   = sqrt(max(diag(C), eps));
    R   = C ./ (d * d');
    kappa = cond(JTJ);
end

function se = delta_se(fun, p_hat, free, C_cov)
% Standard error of a function of the parameters (delta method, numerical
% gradient with respect to the fitted parameters).
    f0 = fun(p_hat);
    g  = zeros(1, numel(p_hat));
    for i = find(free)
        step = 1e-6 * max(abs(p_hat(i)), 1e-3);
        pp = p_hat;
        pp(i) = pp(i) + step;
        g(i) = (fun(pp) - f0) / step;
    end
    n_par = numel(p_hat);
    cov_full = zeros(n_par, n_par);
    cov_full(free, free) = C_cov;
    se = sqrt(max(g * cov_full * g', 0));
end

function T = top_intervals(province, series, dates, obs, dt, z)
% The five reporting intervals with the largest counts, with their length.
    [~, ord] = sort(obs, 'descend');
    k = ord(1:5);
    T = table(repmat(string(province), 5, 1), repmat(string(series), 5, 1), ...
        dates(k), obs(k), dt(k), obs(k) ./ dt(k), z(k), ...
        'VariableNames', {'Province', 'Series', 'Date', 'Count', ...
        'DaysInInterval', 'PerDay', 'ZScore'});
end