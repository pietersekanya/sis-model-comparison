clear; clc; close all;

%%\ 1. SETTINGS

N_values = [5e5, 1e6, 2e6, 4e6, 7e6, 1e7, 1.3e7];

sigma = 1/6.3;
gamma = 1/10;

mu_dem = 1/(62.07*365);

models(1).name = 'Basic SEIR';
models(1).mu = 0;

models(2).name = 'Demographic SEIR';
models(2).mu = mu_dem;

n_starts = 10;

rng(1);

%%\ 2. LOAD DATA

[report_date, cum_obs] = load_cases();

t = days(report_date - report_date(1));
dt = diff(t);

obs_per_day = diff(cum_obs) ./ dt;

nm = numel(models);
nN = numel(N_values);

%%\ 3. STORAGE

beta_results = zeros(nN, nm);
R0_results   = zeros(nN, nm);
rmse_results = zeros(nN, nm);

%%\ 4. REFIT EACH MODEL FOR EACH N

for j = 1:nN

    N = N_values(j);

    fprintf('\nTesting N = %g\n', N);

    for m = 1:nm

        theta = fit_model( ...
            t, cum_obs, N, models(m).mu, ...
            sigma, gamma, n_starts);

        beta_results(j,m) = theta(1);

        R0_results(j,m) = theta(1) / ...
            (gamma + models(m).mu);

        C = simulate_cum( ...
            theta, t, N, models(m).mu, ...
            sigma, gamma);

        pred_per_day = diff(C) ./ dt;

        rmse_results(j,m) = sqrt(mean( ...
            (pred_per_day - obs_per_day).^2));

        fprintf('  %s: beta = %.5f, R0 = %.5f, RMSE = %.3f\n', ...
            models(m).name, ...
            beta_results(j,m), ...
            R0_results(j,m), ...
            rmse_results(j,m));

    end
end

%%\ 5. RESULTS TABLE

sensitivity_results = table( ...
    N_values', ...
    beta_results(:,1), ...
    beta_results(:,2), ...
    R0_results(:,1), ...
    R0_results(:,2), ...
    rmse_results(:,1), ...
    rmse_results(:,2), ...
    'VariableNames', { ...
    'N', ...
    'BasicSEIR_beta', ...
    'DemographicSEIR_beta', ...
    'BasicSEIR_R0', ...
    'DemographicSEIR_R0', ...
    'BasicSEIR_RMSE', ...
    'DemographicSEIR_RMSE'});

disp(sensitivity_results);

%%\ 6. SAVE

if ~isfolder('results')
    mkdir('results');
end

writetable(sensitivity_results, ...
    fullfile('results', 'N_sensitivity.csv'));

%%\ 7. PLOT R0

fig1 = figure('Color','w');

plot(N_values, R0_results(:,1), '-o');
hold on;
plot(N_values, R0_results(:,2), '-s');

xlabel('Assumed at-risk population N');
ylabel('R_0');
legend({models.name}, 'Location', 'best');
title('Sensitivity of R_0 to assumed population size');
grid on;

saveas(fig1, ...
    fullfile('results', 'N_sensitivity_R0.png'));

%%\ 8. PLOT BETA

fig2 = figure('Color','w');

plot(N_values, beta_results(:,1), '-o');
hold on;
plot(N_values, beta_results(:,2), '-s');

xlabel('Assumed at-risk population N');
ylabel('\beta');
legend({models.name}, 'Location', 'best');
title('Sensitivity of fitted \beta to assumed population size');
grid on;

saveas(fig2, ...
    fullfile('results', 'N_sensitivity_beta.png'));

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
