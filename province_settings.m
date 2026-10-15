function S = province_settings()
%PROVINCE_SETTINGS  Population and demographic rates for each province.
%
%   Fill in the numbers in the INPUTS block below. They are converted to the
%   units used by the model (people, and rates PER DAY). This is the only
%   place these values are set.
%
%   Returns a struct array S (one element per province) with fields:
%     name   : display name
%     column : column name in the table from load_province_cases
%     N      : total population (people)
%     b      : birth rate (per person per day)
%     mu     : natural death rate (per person per day)

    %% INPUTS (fill in; keep the order Ituri, North Kivu) --------------
    pop        = [4600499, 9314854];   % TODO: total population N (people)
    birth_rate = [40.88, 40.88];   % TODO: births per 1000 people per YEAR
    death_rate = [8.42, 8.42];   % TODO: deaths per 1000 people per YEAR

    per_people   = 1000;       % the rates above are per 1000 people;
                               % use 1 if your rates are per person per year
    days_in_year = 365;

    % If your source gives a life expectancy L (years) instead of a death
    % rate, the per-day death rate is 1/(L*365).

    %% CONVERT TO PER-DAY RATES ----------------------------------------
    labels  = {'Ituri', 'North Kivu'};
    columns = {'Ituri', 'NorthKivu'};

    for q = 1:2
        S(q).name   = labels{q};
        S(q).column = columns{q};
        S(q).N      = pop(q);
        S(q).b      = birth_rate(q) / per_people / days_in_year;
        S(q).mu     = death_rate(q) / per_people / days_in_year;
    end

    if any(isnan([S.N, S.b, S.mu]))
        error('province_settings:missing', ...
            'Fill in pop, birth_rate and death_rate in province_settings.m.');
    end

    fprintf('Province settings (rates per day):\n');
    for q = 1:2
        fprintf('  %-10s N = %.0f, b = %.3e, mu = %.3e\n', ...
            S(q).name, S(q).N, S(q).b, S(q).mu);
    end
end
