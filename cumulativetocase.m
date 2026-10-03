%% cumulativetocase.m
% DRC 2026 Ebola data: cumulative -> daily cases
%
% Version 1: observed reporting intervals (no interpolation)
% Version 2: complete daily series from linear interpolation of the
%            cumulative cases (IsObserved marks the real report dates)
%
% The data and the cutoff date come from load_cases.m.
% Run from the top-level repo folder.

clear; clc; close all;

%% ============================================================
% 1. LOAD THE CUMULATIVE DATA
% =============================================================
[report_date, cumulative_cases] = load_cases();


%% ============================================================
% 2. VERSION 1: OBSERVED REPORTING INTERVALS
% =============================================================
% new_cases = C(t) - C(previous report), WITHOUT filling in missing dates.
% For example, if 93 cases are reported across a 2-day gap, the script
% records 93 cases over that 2-day reporting interval.

new_cases_interval    = [NaN; diff(cumulative_cases)];
days_between          = [NaN; days(diff(report_date))];
average_cases_per_day = new_cases_interval ./ days_between;

observed_data = table( ...
    report_date, cumulative_cases, new_cases_interval, ...
    days_between, average_cases_per_day, ...
    'VariableNames', {'Date', 'CumulativeCases', ...
    'NewCasesSincePreviousReport', 'DaysSincePreviousReport', ...
    'AverageCasesPerDay'});

disp('============================================================');
disp('VERSION 1: OBSERVED REPORTING INTERVALS');
disp('============================================================');
disp(observed_data);

if ~isfolder('data'), mkdir('data'); end
writetable(observed_data, ...
    fullfile('data', 'DRC_2026_observed_reporting_intervals.csv'));


%% ============================================================
% 3. VERSION 2: COMPLETE DAILY DATASET (INTERPOLATED)
% =============================================================
% The CUMULATIVE cases are interpolated linearly between reports.
% Example: 689 cases on June 11 and 782 on June 13 -> 93 cases over two
% days -> 46.5 cases per day.
%
% IMPORTANT: this is an ASSUMPTION. The data do not say how the cases
% were spread between the report dates. Use IsObserved to keep
% interpolated days out of any fitting or scoring.

full_date = (report_date(1):caldays(1):report_date(end))';

x_obs  = days(report_date - report_date(1));
x_full = days(full_date   - full_date(1));

full_cumulative = interp1(x_obs, cumulative_cases, x_full, 'linear');
full_new_cases  = [NaN; diff(full_cumulative)];
is_observed     = ismember(full_date, report_date);

daily_data = table( ...
    full_date, full_cumulative, full_new_cases, is_observed, ...
    'VariableNames', {'Date', 'EstimatedCumulativeCases', ...
    'EstimatedNewCases', 'IsObserved'});

disp('============================================================');
disp('VERSION 2: COMPLETE DAILY DATASET');
disp('============================================================');
disp(daily_data);

writetable(daily_data, ...
    fullfile('data', 'DRC_2026_complete_daily_interpolated.csv'));

n_interp = sum(~is_observed);
fprintf('\n%d of %d calendar days (%.0f%%) are interpolated.\n', ...
    n_interp, numel(full_date), 100*n_interp/numel(full_date));


%% ============================================================
% 4. PLOT THE CUMULATIVE DATA
% =============================================================
figure;
plot(report_date, cumulative_cases, 'o');
hold on;
plot(full_date, full_cumulative, '-');
xlabel('Date');
ylabel('Cumulative confirmed cases');
title('DRC 2026 cumulative cases');
legend('Observed data', 'Linear interpolation', 'Location', 'northwest');
grid on;


%% ============================================================
% 5. PLOT DAILY CASES
% =============================================================
figure;
plot(report_date, new_cases_interval, 'o-');
xlabel('Date');
ylabel('New cases since previous report');
title('DRC 2026 observed reporting intervals');
grid on;

figure;
plot(full_date, full_new_cases, '-');
xlabel('Date');
ylabel('Estimated new cases per day');
title('DRC 2026 estimated daily cases (interpolated)');
grid on;


%% ============================================================
% 6. CHECK THAT THE INTERPOLATED DATA REACH THE SAME TOTAL
% =============================================================
total_interpolated = sum(full_new_cases, 'omitnan');

fprintf('\n');
fprintf('Initial cumulative cases: %d\n', cumulative_cases(1));
fprintf('Final cumulative cases:   %d\n', cumulative_cases(end));
fprintf('Total increase:           %d\n', ...
    cumulative_cases(end) - cumulative_cases(1));
fprintf('Sum of interpolated daily increases: %.2f\n', total_interpolated);
