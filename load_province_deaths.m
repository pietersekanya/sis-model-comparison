function [D, info] = load_province_deaths()
%LOAD_PROVINCE_DEATHS  Cumulative confirmed deaths for Ituri and North Kivu.
%
%   Reads the health-zone-level deaths file from INRB-UMIE/BDBV2026-Data
%   (insp_sitrep__cumulative_confirmed_deaths__daily.csv, columns:
%   nom, date, cumulative_confirmed_deaths), adds up the health zones in
%   each province, and returns one cumulative series per province.
%
%   The deaths are returned on the SAME report dates, and with the SAME
%   cutoff date, as the case data from load_province_cases.m, so that cases
%   and deaths line up row by row.
%
%   D    : table with columns Date, Ituri, NorthKivu, Unassigned
%   info : diagnostics (names that were excluded, data fixes applied, ...)
%
%   Keep this file in the top-level repo folder, next to load_province_cases.m
%   and the data/ folder.

    %% SETTINGS --------------------------------------------------------
    data_file = fullfile(fileparts(mfilename('fullpath')), 'data', ...
        'insp_sitrep__cumulative_confirmed_deaths__daily.csv');

    % The case data decide the report dates and the cutoff date
    [P_cases, info_cases] = load_province_cases();
    report_dates = P_cases.Date;
    cutoff_date  = info_cases.cutoff_date;

    % IMPORTANT: these lists and spelling variants must be identical to the
    % ones in load_province_cases.m.
    ituri_zones = [ ...
        "Boga", "Bunia", "Gethy", "Komanda", "Nyakunde", "Rwampara", ...
        "Lolwa", "Mambasa", "Mandima", "Nia-Nia", "Bambu", "Damas", ...
        "Drodro", "Fataki", "Jiba", "Kilo", "Linga", "Lita", "Mangala", ...
        "Mongbalu", "Nizi", "Rethy", "Tchomia", "Angumu", "Aungba", ...
        "Kambala", "Logo", "Mahagi", "Nyarambe", "Rimba", "Adi", "Adja", ...
        "Ariwara", "Aru", "Biringi", "Laybo"];

    nkivu_zones = [ ...
        "Goma", "Karisimbi", "Nyiragongo", "Katoyi", "Kirotshe", "Masisi", ...
        "Mweso", "Itebero", "Kibua", "Pinga", "Walikale", "Alimbongo", ...
        "Biena", "Kayna", "Lubero", "Manguredjipa", "Masereka", ...
        "Musienene", "Kalunguta", "Kamango", "Kyondo", "Mabalako", ...
        "Mutwanga", "Oicha", "Vuhovi", "Beni", "Butembo", "Katwa", ...
        "Bambo", "Binza", "Birambizo", "Kibirizi", "Rutshuru", "Rwanguba"];

    alias_from = ["gety", "nyankunde", "mongbwalu", "rumba"];
    alias_to   = ["gethy", "nyakunde", "mongbalu", "rimba"];

    %% READ EVERYTHING AS TEXT (the file has some malformed entries) ---
    T = readtable(data_file, 'Format', '%s%s%s', 'TextType', 'string');
    if width(T) ~= 3
        error('load_province_deaths:format', ...
            'Expected 3 columns: nom, date, cumulative_confirmed_deaths.');
    end
    T.Properties.VariableNames = {'nom', 'date', 'cum'};

    %% CLEAN ------------------------------------------------------------
    key = clean_key(T.nom, alias_from, alias_to);
    key(ismissing(T.nom)) = "na";          % a missing name = unassigned deaths

    % Date: pull out yyyy-mm-dd
    d = datetime(regexp(T.date, '\d{4}-\d{2}-\d{2}', 'match', 'once'), ...
        'InputFormat', 'yyyy-MM-dd');

    % Count: keep the digits only. "ND" (not determined) and "NA" give NaN,
    % and a stray character such as "1]" becomes 1.
    c = str2double(regexp(T.cum, '\d+', 'match', 'once'));

    %% DEATHS REPORT DATES UP TO THE CUTOFF -----------------------------
    deaths_dates = unique(d(~isnat(d) & d <= cutoff_date));
    if isempty(deaths_dates)
        error('load_province_deaths:noDates', ...
            'No deaths reports on or before the cutoff date.');
    end

    ituri_key = clean_key(ituri_zones(:), alias_from, alias_to);
    nkivu_key = clean_key(nkivu_zones(:), alias_from, alias_to);

    in_ituri = ismember(key, ituri_key);
    in_nkivu = ismember(key, nkivu_key);

    %% ADD UP THE HEALTH ZONES IN EACH PROVINCE -------------------------
    % A zone missing on a given date keeps its last reported value, and a
    % zone that has not reported yet (or is "ND") counts as 0.
    ituri = sum_zones(key, d, c, deaths_dates, ituri_key, true);
    nkivu = sum_zones(key, d, c, deaths_dates, nkivu_key, true);

    % Deaths not yet linked to a health zone ("NA" as the name in the file)
    unassigned = sum_zones(key, d, c, deaths_dates, "na", false);

    % Cumulative counts should never fall. Where a province total dips
    % (a correction downwards), keep the highest value reached so far.
    info.decreases_corrected_ituri = sum(diff(ituri) < 0);
    info.decreases_corrected_nkivu = sum(diff(nkivu) < 0);
    ituri = cummax(ituri);
    nkivu = cummax(nkivu);

    %% PUT THE DEATHS ON THE CASE REPORT DATES -------------------------
    % For each case report date, use the latest deaths report on or before
    % it (0 if there is none yet).
    vals = last_on_or_before(deaths_dates, [ituri, nkivu, unassigned], report_dates);

    D = table(report_dates, vals(:, 1), vals(:, 2), vals(:, 3), ...
        'VariableNames', {'Date', 'Ituri', 'NorthKivu', 'Unassigned'});

    %% DIAGNOSTICS ------------------------------------------------------
    info.cutoff_date       = cutoff_date;
    info.last_date_in_file = max(d);
    info.rows_bad_date     = sum(isnat(d));
    info.rows_no_count     = sum(isnan(c));        % "ND", "NA" or no number
    info.excluded_names    = unique(T.nom(~in_ituri & ~in_nkivu));
    info.ituri_never_seen  = ituri_zones(~ismember(ituri_key, key));
    info.nkivu_never_seen  = nkivu_zones(~ismember(nkivu_key, key));
    info.case_dates_without_deaths_report = sum(~ismember(report_dates, deaths_dates));

    fprintf('Deaths: %d deaths report dates up to %s (last date in file %s)\n', ...
        numel(deaths_dates), char(deaths_dates(end), 'yyyy-MM-dd'), ...
        char(info.last_date_in_file, 'yyyy-MM-dd'));
    fprintf('Last case date: deaths Ituri = %d, North Kivu = %d, unassigned = %d\n', ...
        D.Ituri(end), D.NorthKivu(end), D.Unassigned(end));
    fprintf('Decreases corrected: Ituri %d, North Kivu %d\n', ...
        info.decreases_corrected_ituri, info.decreases_corrected_nkivu);
    fprintf('Rows with an unreadable date: %d; rows with no numeric count (ND/NA): %d\n', ...
        info.rows_bad_date, info.rows_no_count);
    fprintf('Case dates with no deaths report on that exact date (value carried forward): %d of %d\n', ...
        info.case_dates_without_deaths_report, numel(report_dates));
end


function key = clean_key(nom, alias_from, alias_to)
% Lower-case letters only, so that "Nia Nia", "Nia-Nia" and "NIA-NIA" match.
    key = lower(strtrim(string(nom)));
    key = regexprep(key, '\(.*?\)', '');    % drop bracketed text
    key = regexprep(key, '[^a-z]', '');     % keep letters only
    [tf, loc] = ismember(key, alias_from);
    key(tf) = alias_to(loc(tf));
end


function total = sum_zones(key, d, c, all_dates, zone_keys, carry_forward)
% Total cumulative count over a set of zones, for every date in all_dates.
    rows = ismember(key, zone_keys) & ismember(d, all_dates);
    k  = key(rows);
    dd = d(rows);
    cc = c(rows);
    zones = unique(k);

    M = nan(numel(all_dates), numel(zones));     % dates x zones
    for z = 1:numel(zones)
        r = (k == zones(z));
        [ud, ~, g] = unique(dd(r));
        v = accumarray(g, cc(r), [], @max);      % one value per date
        [~, loc] = ismember(ud, all_dates);
        M(loc, z) = v;
    end

    if carry_forward && ~isempty(M)
        M = fillmissing(M, 'previous');          % keep the zone's last value
    end
    M(isnan(M)) = 0;                              % not yet reported = 0
    total = sum(M, 2);
end


function v = last_on_or_before(dates_src, vals_src, dates_query)
% For each query date, the row of vals_src at the latest source date that is
% on or before it (zeros if there is no earlier source date).
    v = zeros(numel(dates_query), size(vals_src, 2));
    for i = 1:numel(dates_query)
        k = find(dates_src <= dates_query(i), 1, 'last');
        if ~isempty(k)
            v(i, :) = vals_src(k, :);
        end
    end
end
