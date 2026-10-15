function [P, info] = load_province_cases()
%LOAD_PROVINCE_CASES  Cumulative confirmed cases for Ituri and North Kivu.
%
%   Reads the health-zone-level file from INRB-UMIE/BDBV2026-Data
%   (insp_sitrep__cumulative_confirmed_cases__daily.csv, columns:
%   nom, date, cumulative_confirmed_cases), adds up the health zones in
%   each province, and returns one cumulative series per province.
%
%   P    : table with columns Date, Ituri, NorthKivu, Unassigned
%   info : diagnostics (names that were excluded, data fixes applied, ...)
%
%   Keep this file in the top-level repo folder, next to the data/ folder.
%   This is the ONLY place where the data file, the cutoff date and the
%   health-zone lists are set.

    %% SETTINGS (set once, then do not change) ------------------------
    data_file = fullfile(fileparts(mfilename('fullpath')), 'data', ...
        'insp_sitrep__cumulative_confirmed_cases__daily.csv');

    cutoff_date = datetime(2026, 10, 12);   % TODO: choose once; last report date used

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

    % Spelling variants found in the data file -> spelling used in the lists
    % above. Names are compared after removing capitals, spaces and hyphens.
    %   Gety      = Gethy        Nyankunde = Nyakunde
    %   Mongbwalu = Mongbalu     Rumba     = Rimba (typo on 2026-07-11)
    alias_from = ["gety", "nyankunde", "mongbwalu", "rumba"];
    alias_to   = ["gethy", "nyakunde", "mongbalu", "rimba"];

    %% READ EVERYTHING AS TEXT (the file has a few malformed entries) --
    T = readtable(data_file, 'Format', '%s%s%s', 'TextType', 'string');
    if width(T) ~= 3
        error('load_province_cases:format', ...
            'Expected 3 columns: nom, date, cumulative_confirmed_cases.');
    end
    T.Properties.VariableNames = {'nom', 'date', 'cum'};

    %% CLEAN ----------------------------------------------------------
    key = clean_key(T.nom, alias_from, alias_to);
    key(ismissing(T.nom)) = "na";          % a missing name = unassigned cases

    % Date: pull out yyyy-mm-dd (rescues entries such as "2026-06-25]")
    d = datetime(regexp(T.date, '\d{4}-\d{2}-\d{2}', 'match', 'once'), ...
        'InputFormat', 'yyyy-MM-dd');

    % Count: keep the digits only ("17-" becomes 17; "ND" = not determined -> NaN)
    c = str2double(regexp(T.cum, '\d+', 'match', 'once'));

    %% REPORT DATES UP TO THE CUTOFF -----------------------------------
    all_dates = unique(d(~isnat(d) & d <= cutoff_date));

    ituri_key = clean_key(ituri_zones(:), alias_from, alias_to);
    nkivu_key = clean_key(nkivu_zones(:), alias_from, alias_to);

    in_ituri = ismember(key, ituri_key);
    in_nkivu = ismember(key, nkivu_key);

    %% ADD UP THE HEALTH ZONES IN EACH PROVINCE ------------------------
    % A zone missing on a given date keeps its last reported value, and a
    % zone that has not reported yet counts as 0.
    ituri = sum_zones(key, d, c, all_dates, ituri_key, true);
    nkivu = sum_zones(key, d, c, all_dates, nkivu_key, true);

    % Cases not yet linked to a health zone ("NA" in the file). These are not
    % carried forward: they disappear once the cases are assigned to zones.
    unassigned = sum_zones(key, d, c, all_dates, "na", false);

    % Cumulative counts should never fall. Where a province total dips
    % (because a zone's count was corrected downwards), keep the highest
    % value reached so far.
    info.decreases_corrected_ituri = sum(diff(ituri) < 0);
    info.decreases_corrected_nkivu = sum(diff(nkivu) < 0);
    ituri = cummax(ituri);
    nkivu = cummax(nkivu);

    %% OUTPUT ---------------------------------------------------------
    P = table(all_dates, ituri, nkivu, unassigned, ...
        'VariableNames', {'Date', 'Ituri', 'NorthKivu', 'Unassigned'});

    info.cutoff_date       = cutoff_date;
    info.last_date_in_file = max(d);
    info.rows_bad_date     = sum(isnat(d));        % rows whose date could not be read
    info.rows_no_count     = sum(isnan(c));        % rows with "ND" or no number
    info.excluded_names    = unique(T.nom(~in_ituri & ~in_nkivu));
    info.ituri_never_seen  = ituri_zones(~ismember(ituri_key, key));
    info.nkivu_never_seen  = nkivu_zones(~ismember(nkivu_key, key));

    fprintf('Data: %d report dates from %s to %s (cutoff %s; last date in file %s)\n', ...
        numel(all_dates), char(all_dates(1), 'yyyy-MM-dd'), ...
        char(all_dates(end), 'yyyy-MM-dd'), char(cutoff_date, 'yyyy-MM-dd'), ...
        char(info.last_date_in_file, 'yyyy-MM-dd'));
    fprintf('Last date: Ituri = %d, North Kivu = %d, unassigned = %d\n', ...
        ituri(end), nkivu(end), unassigned(end));
    fprintf('Decreases corrected: Ituri %d, North Kivu %d\n', ...
        info.decreases_corrected_ituri, info.decreases_corrected_nkivu);
    fprintf('Rows with an unreadable date: %d; rows with no numeric count: %d\n', ...
        info.rows_bad_date, info.rows_no_count);
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
% Total cumulative cases over a set of zones, for every date in all_dates.
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
    M(isnan(M)) = 0;                              % not yet reported = 0 cases
    total = sum(M, 2);
end
