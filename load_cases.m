function [report_date, cum_cases] = load_cases()
%LOAD_CASES  Load, check and trim the DRC cumulative case data.
%
%   Source: INRB-UMIE/BDBV2026-Data (GitHub), tag build-2026-10-01-b077a6c
%   data/insp_sitrep/processed/
%   insp_sitrep__national_cumulative_confirmed_cases__daily.csv
%   Accessed 3 October 2026.
%
%   This is the ONLY place where the data file and the analysis cutoff
%   date are set, so every script uses exactly the same data.
%   Keep this file in the top-level repo folder, next to the data/ folder.
%
%   Returns column vectors:
%     report_date : datetime, the dates on which cumulative cases were reported
%     cum_cases   : cumulative confirmed cases on those dates

    %% SETTINGS (set once, then do not change)
    data_file = fullfile(fileparts(mfilename('fullpath')), 'data', ...
        'insp_sitrep__national_cumulative_confirmed_cases__daily.csv');

    cutoff_date = datetime(2026, 9, 21);   % TODO: choose once; last report date used

    %% READ
    T = readtable(data_file);

    % The file should hold a single national series
    if numel(unique(T.nom)) > 1
        error('load_cases:multipleSeries', ...
            '"nom" has more than one value - filter to the national series first.');
    end

    if ~isdatetime(T.date)
        T.date = datetime(T.date, 'InputFormat', 'yyyy-MM-dd');
    end
    T = sortrows(T, 'date');

    %% APPLY THE CUTOFF (freeze the data used in the project)
    last_in_file = T.date(end);
    T = T(T.date <= cutoff_date, :);
    if height(T) < 2
        error('load_cases:tooFewRows', ...
            'Fewer than 2 reports on or before the cutoff date.');
    end

    %% CHECKS
    [~, keep] = unique(T.date);              % drop any duplicated dates
    T = T(keep, :);
    if any(diff(T.national_cumulative_confirmed_cases) < 0)
        error('load_cases:decreasing', ...
            'Cumulative cases decrease somewhere - check the raw data.');
    end

    report_date = T.date;
    cum_cases   = T.national_cumulative_confirmed_cases;

    fprintf('Data: %d reports from %s to %s (cutoff %s; last date in file %s)\n', ...
        numel(report_date), char(report_date(1), 'yyyy-MM-dd'), ...
        char(report_date(end), 'yyyy-MM-dd'), char(cutoff_date, 'yyyy-MM-dd'), ...
        char(last_in_file, 'yyyy-MM-dd'));
end
