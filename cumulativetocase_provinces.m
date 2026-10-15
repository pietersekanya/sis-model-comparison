
%% cumulativetocase_provinces.m
%
% Ituri and North Kivu:
% Convert cumulative reported cases into observed reporting intervals.
%
% IMPORTANT:
% No missing calendar days are filled in and no zero cases are invented.
% New cases are calculated only between consecutive observed reports.
%
% For example:
%   Report date 1 -> Report date 2 = observed cases during that interval.
%
% The time intervals may therefore be uneven.
%
% The data, zone lists and cutoff date come from
% load_province_cases.m.
%
% Run from the top-level repository folder.

clear; clc; close all;

%% ============================================================
% 1. LOAD THE PROVINCE DATA
% =============================================================

[P, ~] = load_province_cases();

report_date = P.Date;
cum         = [P.Ituri, P.NorthKivu];

names = {'Ituri', 'North Kivu'};

%% ============================================================
% 2. OBSERVED REPORTING INTERVALS
% =============================================================

% IMPORTANT:
% We use ONLY the actual reporting dates in the data.
%
% If there is no report on a calendar day, that day is NOT added
% to the dataset and is NOT assigned zero cases.
%
% New cases are calculated between consecutive reports.

new_interval = [NaN(1,2); diff(cum)];

% Number of calendar days between consecutive reports.
days_between = [NaN; days(diff(report_date))];

% Average cases per day over each observed reporting interval.
%
% This is descriptive only and is NOT used as the model-fitting target.
avg_per_day = new_interval ./ days_between;

observed_data = table( ...
    report_date, ...
    days_between, ...
    cum(:,1), ...
    cum(:,2), ...
    new_interval(:,1), ...
    new_interval(:,2), ...
    avg_per_day(:,1), ...
    avg_per_day(:,2), ...
    'VariableNames', { ...
    'Date', ...
    'DaysSincePreviousReport', ...
    'Ituri_Cumulative', ...
    'NorthKivu_Cumulative', ...
    'Ituri_NewCases', ...
    'NorthKivu_NewCases', ...
    'Ituri_AverageCasesPerDay', ...
    'NorthKivu_AverageCasesPerDay'});

%% ============================================================
% 3. DISPLAY OBSERVED DATA
% =============================================================

disp('============================================================');
disp('OBSERVED REPORTING INTERVALS');
disp('============================================================');

disp(observed_data);

%% ============================================================
% 4. SAVE OBSERVED DATA
% =============================================================

if ~isfolder('data')
    mkdir('data');
end

if ~isfolder('results')
    mkdir('results');
end

writetable(observed_data, ...
    fullfile('data', 'province_observed_reporting_intervals.csv'));

%% ============================================================
% 5. REPORTING-INTERVAL DIAGNOSTICS
% =============================================================

valid_intervals = days_between(2:end);

fprintf('\n============================================================\n');
fprintf('REPORTING INTERVAL SUMMARY\n');
fprintf('============================================================\n');

fprintf('Number of report dates: %d\n', numel(report_date));
fprintf('Number of observed intervals: %d\n', numel(valid_intervals));
fprintf('Minimum interval length: %.0f day(s)\n', min(valid_intervals));
fprintf('Maximum interval length: %.0f day(s)\n', max(valid_intervals));
fprintf('Mean interval length: %.2f days\n', mean(valid_intervals));
fprintf('Median interval length: %.2f days\n', median(valid_intervals));

%% ============================================================
% 6. CHECK FOR NEGATIVE INTERVAL COUNTS
% =============================================================

for p = 1:2

    negative_intervals = new_interval(:,p) < 0;

    fprintf('\n%s:\n', names{p});
    fprintf('  Negative reporting intervals: %d\n', ...
        sum(negative_intervals, 'omitnan'));

    fprintf('  Total increase over analysis period: %.0f\n', ...
        cum(end,p) - cum(1,p));

end

%% ============================================================
% 7. PLOT CUMULATIVE CASES AT OBSERVED REPORT DATES
% =============================================================

fig1 = figure('Color','w');

plot(report_date, cum(:,1), 'o-');
hold on;
plot(report_date, cum(:,2), 's-');

xlabel('Report date');
ylabel('Cumulative confirmed cases');

legend(names, 'Location','northwest');

title('Cumulative confirmed cases at observed reporting dates');

grid on;

saveas(fig1, ...
    fullfile('results','province_cumulative_observed.png'));

%% ============================================================
% 8. PLOT OBSERVED INTERVAL CASES
% =============================================================

fig2 = figure('Color','w');

plot(report_date, new_interval(:,1), 'o-');
hold on;
plot(report_date, new_interval(:,2), 's-');

xlabel('Report date');
ylabel('New cases since previous report');

legend(names, 'Location','northwest');

title('New cases between consecutive reporting dates');

grid on;

saveas(fig2, ...
    fullfile('results','province_new_cases_by_reporting_interval.png'));

%% ============================================================
% IMPORTANT NOTE
% ============================================================
%
% The average cases-per-day values above are provided only to describe
% the uneven reporting intervals.
%
% The SEIR model is NOT fitted to these average values.
%
% Model fitting should use:
%
%   observed interval cases = diff(cumulative observed cases)
%
% and compare them with:
%
%   model interval cases = diff(model cumulative cases)
%
% at the actual reporting times.

