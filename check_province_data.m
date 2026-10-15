%% check_province_data.m
% Builds the Ituri and North Kivu cumulative series from the health-zone
% file, prints diagnostics, plots the result and saves it for later scripts.
% Run from the top-level repo folder.

clear; clc; close all;

[P, info] = load_province_cases();

disp('Last rows of the province series:');
disp(P(end-4:end, :));

disp('Names in the file that were NOT counted (other provinces, or unassigned "NA"):');
disp(info.excluded_names);

disp('Ituri zones in your list that never appear in the file (fine if they have no cases yet; otherwise check the spelling):');
disp(info.ituri_never_seen);

disp('North Kivu zones in your list that never appear in the file:');
disp(info.nkivu_never_seen);

%% Plot
fig = figure('Color', 'w');
plot(P.Date, P.Ituri, '-o'); hold on;
plot(P.Date, P.NorthKivu, '-s');
plot(P.Date, P.Unassigned, '--');
xlabel('Date');
ylabel('Cumulative confirmed cases');
legend({'Ituri', 'North Kivu', 'Unassigned (no health zone yet)'}, ...
    'Location', 'northwest');
title('Cumulative confirmed cases by province');
grid on;

%% Save
if ~isfolder('data'),    mkdir('data');    end
if ~isfolder('results'), mkdir('results'); end
writetable(P, fullfile('data', 'province_cumulative_cases.csv'));
saveas(fig, fullfile('results', 'province_cumulative_cases.png'));
