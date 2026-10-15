%% check_province_deaths.m
% Builds the Ituri and North Kivu cumulative confirmed deaths from the
% health-zone file, lines them up with the confirmed cases on the same
% report dates, prints diagnostics, plots the result and saves it as CSV
% files for later scripts. Run from the top-level repo folder.

clear; clc; close all;

[P, ~]    = load_province_cases();     % cumulative cases: Date, Ituri, NorthKivu, Unassigned
[D, info] = load_province_deaths();    % cumulative deaths on the same dates

%% Diagnostics
disp('Last rows of the cumulative deaths:');
disp(D(end-4:end, :));

disp('Names in the deaths file that were NOT counted (other provinces, or unassigned "NA"):');
disp(info.excluded_names);

disp('Ituri zones in your list that never appear in the deaths file (fine if they have no deaths yet):');
disp(info.ituri_never_seen);

disp('North Kivu zones in your list that never appear in the deaths file:');
disp(info.nkivu_never_seen);

%% New cases and deaths in each reporting interval
days_between = [NaN; days(diff(P.Date))];
new_cases    = [NaN(1, 2); diff([P.Ituri, P.NorthKivu])];
new_deaths   = [NaN(1, 2); diff([D.Ituri, D.NorthKivu])];

combined = table(P.Date, days_between, ...
    P.Ituri, P.NorthKivu, D.Ituri, D.NorthKivu, ...
    new_cases(:, 1), new_cases(:, 2), new_deaths(:, 1), new_deaths(:, 2), ...
    'VariableNames', {'Date', 'DaysSincePreviousReport', ...
    'Ituri_CumulativeCases', 'NorthKivu_CumulativeCases', ...
    'Ituri_CumulativeDeaths', 'NorthKivu_CumulativeDeaths', ...
    'Ituri_NewCases', 'NorthKivu_NewCases', ...
    'Ituri_NewDeaths', 'NorthKivu_NewDeaths'});

%% Crude deaths per confirmed case on the last date
% Descriptive only: deaths lag cases and not every case is confirmed, so this
% is not a true case fatality ratio.
names = {'Ituri', 'North Kivu'};
x = [D.Ituri(end), D.NorthKivu(end)];
n = [P.Ituri(end), P.NorthKivu(end)];

fprintf('\nCrude deaths per confirmed case on %s:\n', char(P.Date(end), 'yyyy-MM-dd'));
for q = 1:2
    ci = wilson_ci(x(q), n(q));
    fprintf('  %-10s %d deaths / %d cases = %.3f (95%% CI %.3f to %.3f)\n', ...
        names{q}, x(q), n(q), x(q) / n(q), ci(1), ci(2));
end

%% Plots
fig1 = figure('Color', 'w');
plot(D.Date, D.Ituri, '-o'); hold on;
plot(D.Date, D.NorthKivu, '-s');
plot(D.Date, D.Unassigned, '--');
xlabel('Date');
ylabel('Cumulative confirmed deaths');
legend({'Ituri', 'North Kivu', 'Unassigned (no health zone)'}, 'Location', 'northwest');
title('Cumulative confirmed deaths by province');
grid on;

ratio = [D.Ituri ./ P.Ituri, D.NorthKivu ./ P.NorthKivu];
ratio([P.Ituri, P.NorthKivu] == 0) = NaN;      % no cases yet

fig2 = figure('Color', 'w');
plot(P.Date, ratio(:, 1), '-o'); hold on;
plot(P.Date, ratio(:, 2), '-s');
xlabel('Date');
ylabel('Cumulative deaths per confirmed case');
legend(names, 'Location', 'best');
title('Crude deaths per confirmed case over time');
grid on;

%% Save
if ~isfolder('data'),    mkdir('data');    end
if ~isfolder('results'), mkdir('results'); end
writetable(D, fullfile('data', 'province_cumulative_deaths.csv'));
writetable(combined, fullfile('data', 'province_cases_and_deaths.csv'));
saveas(fig1, fullfile('results', 'province_cumulative_deaths.png'));
saveas(fig2, fullfile('results', 'province_deaths_per_case.png'));

fprintf('\nSaved data/province_cumulative_deaths.csv and data/province_cases_and_deaths.csv\n');


%% LOCAL FUNCTIONS -------------------------------------------------
function ci = wilson_ci(x, n)
% 95% Wilson score interval for a proportion x/n.
z = 1.96;
p = x / n;
centre = (p + z^2 / (2*n)) / (1 + z^2 / n);
half   = z * sqrt(p*(1 - p)/n + z^2 / (4*n^2)) / (1 + z^2 / n);
ci = [centre - half, centre + half];
end