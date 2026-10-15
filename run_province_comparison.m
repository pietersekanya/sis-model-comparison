
%% run_province_comparison.m
% Fit a demographic SIS model separately to Ituri and North Kivu.
% Include parameter estimation, sensitivity analysis and
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

%% 3. FIT THE MAIN MODEL TO EACH PROVINCE

% Results include parameter estimates, derived quantities,
% and in-sample MAE and RMSE for cases and deaths.

results = table();

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
    p_hat = fit_sis_model( ...
        tw, cases, deaths, S(q).N, S(q).b, S(q).mu, ...
        p_fixed, free, lb, ub, n_starts);

    betaN = p_hat(1);
    gamma = p_hat(2);
    mu_tilde = p_hat(3);
    I0 = p_hat(4);

    beta = betaN / S(q).N;
    R_start = betaN / (gamma + mu_tilde);
    CFR = mu_tilde / (gamma + mu_tilde);

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

    % Print estimates and fit metrics.
    fprintf('betaN       = %.5g per day\n', betaN);
    fprintf('beta        = %.5g per person per day\n', beta);
    fprintf('gamma       = %.5g per day (fixed)\n', gamma);
    fprintf('mu_tilde    = %.5g per day\n', mu_tilde);
    fprintf('I0          = %.5g people\n', I0);
    fprintf('R_start     = %.5g\n', R_start);
    fprintf('Model CFR   = %.2f%%\n', 100*CFR);

    fprintf('Cases:  MAE = %.3f, RMSE = %.3f\n', ...
        MAE_cases, RMSE_cases);
    fprintf('Deaths: MAE = %.3f, RMSE = %.3f\n', ...
        MAE_deaths, RMSE_deaths);

    % Store results.
    new_row = table(string(S(q).name), betaN, beta, gamma, ...
        mu_tilde, I0, R_start, CFR, ...
        MAE_cases, RMSE_cases, MAE_deaths, RMSE_deaths, ...
        'VariableNames', {'Province', 'betaN', 'beta', ...
        'gamma', 'mu_tilde', 'I0', 'R_start', 'CFR', ...
        'MAE_cases', 'RMSE_cases', 'MAE_deaths', 'RMSE_deaths'});

    results = [results; new_row];

    % Save plotting data.
    fit_dates{q} = report_date(s+1:end);
    fit_obs_c{q} = obs_cases;
    fit_pred_c{q} = pred_cases;
    fit_obs_d{q} = obs_deaths;
    fit_pred_d{q} = pred_deaths;

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
mae_cv_d = nan(n_origins, np);
rmse_cv_d = nan(n_origins, np);

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

        mae_cv_d(k,q) = mean(abs(err_d));
        rmse_cv_d(k,q) = sqrt(mean(err_d.^2));

    end

    fprintf('Finished CV for %s.\n', S(q).name);
end

% Average performance across valid forecast origins.
cv_summary = table( ...
    string({S.name})', ...
    mean(mae_cv_c, 1, 'omitnan')', ...
    mean(rmse_cv_c, 1, 'omitnan')', ...
    mean(mae_cv_d, 1, 'omitnan')', ...
    mean(rmse_cv_d, 1, 'omitnan')', ...
    sum(~isnan(mae_cv_c), 1)', ...
    'VariableNames', {'Province', ...
    'CV_MAE_cases', 'CV_RMSE_cases', ...
    'CV_MAE_deaths', 'CV_RMSE_deaths', ...
    'Origins_used'});

disp(cv_summary);

% Save errors for each forecast origin.
cv_by_origin = table( ...
    report_date(origins), ...
    mae_cv_c(:,1), mae_cv_c(:,2), ...
    rmse_cv_c(:,1), rmse_cv_c(:,2), ...
    mae_cv_d(:,1), mae_cv_d(:,2), ...
    rmse_cv_d(:,1), rmse_cv_d(:,2), ...
    'VariableNames', {'ForecastOrigin', ...
    'MAE_cases_Ituri', 'MAE_cases_NorthKivu', ...
    'RMSE_cases_Ituri', 'RMSE_cases_NorthKivu', ...
    'MAE_deaths_Ituri', 'MAE_deaths_NorthKivu', ...
    'RMSE_deaths_Ituri', 'RMSE_deaths_NorthKivu'});

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

writetable(sens_table, fullfile('results', ...
    'sis_province_gamma_sensitivity.csv'));

writetable(cv_summary, fullfile('results', ...
    'sis_province_cv_summary.csv'));

writetable(cv_by_origin, fullfile('results', ...
    'sis_province_cv_by_origin.csv'));

saveas(fig_fit, fullfile('results', ...
    'sis_province_fits.png'));

saveas(fig_cv, fullfile('results', ...
    'sis_province_cv_error.png'));

fprintf('\nAll results saved in the results/ folder.\n');
