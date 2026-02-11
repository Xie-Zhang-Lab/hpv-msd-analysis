function run_msd_analysis()
% msd_piecewise_multi_spots
% TrackMate Spots CSV -> clean tracks -> ensemble MSD -> power-law fits
%
% Outputs (current folder):
%   - msd_analysis_powerlaw_fits.png
%   - msd_analysis_powerlaw_fit_params.csv
%
% Requires: msdanalyzer on MATLAB path.

%% =========================== USER INPUTS ================================

fnames = { ...
    'RKDA3-4.nd2 - C=1-1_spots', ...    % Brownian
    %'RKD-A3-D.nd2 - C=1-1_spots', ...  % Confined
    %'RKD-A4-C.nd2 - C=1-4_spots', ...  % Directed
    %'C2-Ctrl-A1-1-1_spots', ...        % Control 1
    %'Ctrl-A1-4-1_spots' ...            % Control 2
    %'C2-Ctrl-A2-A-2_spots'             % Control 3 - Adjusted X & Y axis
    
};

frameInterval_s = 5;

isPixels     = false;
pixelSize_um = 0.108;

ROUND_TIME_TO_GRID   = true;
ROUND_GRID_CANDIDATE = [5, 10];

MAX_LAG_TIME_S = [];   % [] disables

MAX_SEGMENTS = 3;
MIN_SEG_LEN  = 5;
SLOPE_WIN    = 5;
SLOPE_TOL    = 0.3;

MIN_TRACKS_PER_LAG = 20;

USE_LOGLOG       = false;
units_space      = 'um';
units_time       = 's';
fontSize         = 18;
AXIS_LINE_WIDTH  = 6.0;
MAIN_MARKER_SIZE = 20;
STAGE_FIT_WIDTH  = 7.5;

colors = [ ...
    %58/255 9/255 235/255   % blue
    %0.20 0.60 0.25         % green
    1, 0, 0;                % red
];

USE_CUSTOM_MSD_AXES = true;

MSD_X_LIM   = [0 30];
MSD_Y_LIM   = [0 1];
MSD_X_TICKS = 0:5:30;
MSD_Y_TICKS = [0 0.5 1.0];

APPLY_CUSTOM_TICK_FORMAT = true;
X_TICK_FORMAT            = '%.0f';
Y_TICK_FORMAT            = '%.1f';

% ---- Output filenames (as requested) ----
OUT_FIG_PNG   = 'msd_analysis_powerlaw_fits.png';
OUT_PARAMSCSV = 'msd_analysis_powerlaw_fit_params.csv';

%% ================== READ AND POOL SPOTS FILES ============================

if isempty(fnames)
    error('fnames is empty. Add at least one TrackMate Spots CSV (with or without .csv).');
end

TRACK_all   = [];
X_all       = [];
Y_all       = [];
Tsec_all    = [];
trackOffset = 0;

for f = 1:numel(fnames)
    fname = fnames{f};
    if ~contains(fname, '.csv')
        fname = [fname '.csv'];
    end

    fprintf('Reading %s\n', fname);

    raw = readlines(fname);
    if numel(raw) < 6
        error('File %s too short.', fname);
    end

    if contains(raw(2), sprintf('\t'))
        delim = sprintf('\t');
    else
        delim = ',';
    end

    varNames = ["Label","SpotID","TrackID","Quality","X","Y","Z","T","Frame","Radius", ...
                "Visibility","ManualSpotColor","MeanIntensityCh1","MedianIntensityCh1", ...
                "MinIntensityCh1","MaxIntensityCh1","SumIntensityCh1","StdIntensityCh1", ...
                "ContrastCh1","SignalNoiseRatioCh1"];

    opts = detectImportOptions(fname,'FileType','text','Delimiter',delim);
    opts.VariableNames     = varNames;
    opts.VariableNamesLine = 2;
    opts.DataLines         = [5 Inf];

    T = readtable(fname, opts);

    TRACK = toNum(T.TrackID);
    X     = toNum(T.X);
    Y     = toNum(T.Y);

    hasT     = ismember('T',     T.Properties.VariableNames);
    hasFrame = ismember('Frame', T.Properties.VariableNames);

    if hasT
        Tsec_col = toNum(T.T);
    else
        Tsec_col = nan(height(T),1);
    end

    if hasFrame
        Frame_col = toNum(T.Frame);
    else
        Frame_col = nan(height(T),1);
    end

    valid = ~isnan(TRACK) & ~isnan(X) & ~isnan(Y) & ...
            ((hasT & ~isnan(Tsec_col)) | (hasFrame & ~isnan(Frame_col)));

    TRACK = TRACK(valid);
    X     = X(valid);
    Y     = Y(valid);
    Tsec  = Tsec_col(valid);
    FRAME = Frame_col(valid);

    if hasT && any(~isnan(Tsec))
        % Use T as-is (seconds)
    else
        if ~(hasFrame && any(~isnan(FRAME)))
            error('No valid time column in %s', fname);
        end
        Tsec = FRAME * frameInterval_s;
    end

    if ROUND_TIME_TO_GRID && numel(Tsec) > 1
        tuniq = unique(Tsec(:));
        if numel(tuniq) > 1
            dts   = diff(tuniq);
            dtMed = median(dts);
            [~, idxClosest] = min(abs(ROUND_GRID_CANDIDATE - dtMed));
            dt_round = ROUND_GRID_CANDIDATE(idxClosest);
            fprintf('  Estimated Δt ~ %.6g s, rounding to grid of %.3g s\n', dtMed, dt_round);
            Tsec = round(Tsec / dt_round) * dt_round;
        end
    end

    if isPixels
        X = X * pixelSize_um;
        Y = Y * pixelSize_um;
    end

    TRACK = TRACK + trackOffset;
    trackOffset = max(TRACK);

    TRACK_all = [TRACK_all; TRACK]; %#ok<AGROW>
    X_all     = [X_all; X]; %#ok<AGROW>
    Y_all     = [Y_all; Y]; %#ok<AGROW>
    Tsec_all  = [Tsec_all; Tsec]; %#ok<AGROW>
end

TRACK = TRACK_all;
X     = X_all;
Y     = Y_all;
Tsec  = Tsec_all;

%% ======================== BUILD CLEAN TRACKS ==============================

uTracks = unique(TRACK);
tracks  = cell(0,1);

for k = 1:numel(uTracks)
    idx = (TRACK == uTracks(k));
    tk  = Tsec(idx);
    xk  = X(idx);
    yk  = Y(idx);

    [tk, ord] = sort(tk);
    xk = xk(ord);
    yk = yk(ord);

    [tkU, ia] = unique(tk, 'stable');
    xkU = xk(ia);
    ykU = yk(ia);

    if numel(tkU) < 2
        continue;
    end

    tracks{end+1,1} = [tkU, xkU, ykU]; %#ok<AGROW>
end

fprintf('Tracks after cleaning: %d\n', numel(tracks));
if isempty(tracks)
    error('No valid tracks with >=2 points.');
end

%% ========================= ENSEMBLE MSD ===================================

ma = msdanalyzer(2, units_space, units_time);
ma = ma.addAll(tracks);
ma = ma.computeMSD;

props = properties(ma);
if any(strcmpi(props,'msd'))
    msdCell = ma.msd;
else
    msdCell = ma.MSD;
end
msdCell = msdCell(~cellfun(@isempty, msdCell));
nTracks = numel(msdCell);

all_t = [];
for k = 1:nTracks
    Mk = msdCell{k};
    all_t = [all_t; Mk(:,1)]; %#ok<AGROW>
end
all_t    = all_t(isfinite(all_t) & all_t > 0);
t_unique = sort(unique(all_t(:)), 'ascend');

S = NaN(numel(t_unique), nTracks);
for k = 1:nTracks
    Mk = msdCell{k};
    tk = Mk(:,1);
    yk = Mk(:,2);

    good = isfinite(tk) & isfinite(yk) & tk > 0 & yk > 0;
    tk = tk(good);
    yk = yk(good);

    [tk, ord] = sort(tk);
    yk = yk(ord);

    [~, ia, ib] = intersect(round(t_unique,12), round(tk,12), 'stable');
    if ~isempty(ia)
        S(ia,k) = yk(ib);
    end
end

mu        = mean(S, 2, 'omitnan');
n_per_lag = sum(isfinite(S), 2);  % used for filtering only (not exported)

valid = isfinite(mu) & isfinite(t_unique) & t_unique > 0 & mu > 0;
t = t_unique(valid);
y = mu(valid);
n_per_lag = n_per_lag(valid);

if ~isempty(MAX_LAG_TIME_S) && isfinite(MAX_LAG_TIME_S)
    mask = (t <= MAX_LAG_TIME_S);
    t = t(mask);
    y = y(mask);
    n_per_lag = n_per_lag(mask);
end

lastIdx = find(n_per_lag >= MIN_TRACKS_PER_LAG, 1, 'last');
if isempty(lastIdx)
    warning('No lag has >= %d tracks; using all lags within max lag cutoff.', MIN_TRACKS_PER_LAG);
else
    t = t(1:lastIdx);
    y = y(1:lastIdx);
end

if isempty(t)
    error(['After applying MAX_LAG_TIME_S and MIN_TRACKS_PER_LAG, no valid lag points remain. ' ...
           'Try lowering MIN_TRACKS_PER_LAG or increasing MAX_LAG_TIME_S.']);
end

%% ======================== SEGMENTATION ===================================

logt = log10(t);
logy = log10(y);

localSlope = local_loglog_slope(logt, logy, SLOPE_WIN);
bounds     = auto_segment_bounds(localSlope, MIN_SEG_LEN, MAX_SEGMENTS, SLOPE_TOL);
segments   = build_segments(numel(t), bounds, MIN_SEG_LEN);

%% ============================= PLOTTING ===================================

f = figure('Color','w','Position',[150 150 650 480]);
ax = gca; hold(ax,'on');

plot(ax, t, y, 'o', ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'MarkerSize', MAIN_MARKER_SIZE);

if USE_LOGLOG
    set(ax,'XScale','log','YScale','log');
end

set(ax,'Color','w','XColor','k','YColor','k','FontSize',fontSize, ...
    'LineWidth',AXIS_LINE_WIDTH,'TickDir','out');

xlabel(ax, sprintf('\\Delta t (%s)', units_time), 'FontSize',fontSize,'Color','k');
ylabel(ax, sprintf('MSD (%s^2)',    units_space), 'FontSize',fontSize,'Color','k');

fitHandles  = gobjects(0);
stageLabels = {};
paramRows   = [];

for s = 1:min(numel(segments), size(colors,1))
    idxS = segments{s};
    tS   = t(idxS);
    yS   = y(idxS);

    p      = polyfit(log10(tS), log10(yS), 1);
    alphaS = p(1);
    AS     = 10.^p(2);

    tfit = linspace(min(tS), max(tS), 200);
    yfit = AS .* (tfit.^alphaS);

    h = plot(ax, tfit, yfit, '-', ...
        'Color', colors(s,:), ...
        'LineWidth', STAGE_FIT_WIDTH);
    fitHandles(end+1) = h; %#ok<AGROW>

    stageLabels{end+1} = sprintf('Stage %d (\\alpha=%.2f)', s, alphaS); %#ok<AGROW>

    paramRows = [paramRows; s, idxS(1), idxS(end), AS, alphaS]; %#ok<AGROW>
end

if USE_CUSTOM_MSD_AXES
    if ~isempty(MSD_X_LIM),   xlim(ax, MSD_X_LIM);   end
    if ~isempty(MSD_Y_LIM),   ylim(ax, MSD_Y_LIM);   end
    if ~isempty(MSD_X_TICKS), xticks(ax, MSD_X_TICKS); end
    if ~isempty(MSD_Y_TICKS), yticks(ax, MSD_Y_TICKS); end
end

if ~isempty(fitHandles)
    legend(fitHandles, stageLabels, 'Location','southeast','Box','off');
end

if APPLY_CUSTOM_TICK_FORMAT && ~USE_LOGLOG
    xt = xticks(ax);
    xticklabels(ax, arrayfun(@(v) sprintf(X_TICK_FORMAT, v), xt, 'UniformOutput', false));
    yt = yticks(ax);
    yticklabels(ax, arrayfun(@(v) sprintf(Y_TICK_FORMAT, v), yt, 'UniformOutput', false));
end

% Save PNG only (no PDFs)
print(f, '-dpng', OUT_FIG_PNG, '-r300');

% Save parameters CSV only
if ~isempty(paramRows)
    paramTbl = array2table(paramRows, ...
        'VariableNames',{'stage','idx_start','idx_end','A','alpha'});
    writetable(paramTbl, OUT_PARAMSCSV);
end

fprintf('Saved %s and %s\n', OUT_FIG_PNG, OUT_PARAMSCSV);

end % ========================= END MAIN =====================================


%% ========================== HELPER FUNCTIONS ==============================

function v = toNum(c)
if isnumeric(c)
    v = double(c);
elseif isstring(c)
    v = str2double(c);
elseif iscellstr(c)
    v = str2double(string(c));
else
    v = str2double(string(c));
end
end

function slope = local_loglog_slope(logt, logy, win)
logt = logt(:); logy = logy(:);
L = numel(logt);
slope = nan(L,1);
if L < win
    return;
end
for i = 1:(L-win+1)
    idx = i:(i+win-1);
    p = polyfit(logt(idx), logy(idx), 1);
    slope(i+floor(win/2)) = p(1);
end
first = find(isfinite(slope),1,'first');
last  = find(isfinite(slope),1,'last');
if ~isempty(first), slope(1:first-1) = slope(first); end
if ~isempty(last),  slope(last+1:end) = slope(last);  end
end

function bounds = auto_segment_bounds(slope, minLen, maxSeg, tol)
L = numel(slope);
if L < 2*minLen
    bounds = [];
    return;
end
ds   = abs(diff(slope));
cand = find(ds > tol) + 1;

bounds = [];
last   = minLen;

for c = cand(:)'
    if c - last >= minLen && numel(bounds) < (maxSeg-1)
        bounds(end+1) = c; %#ok<AGROW>
        last = c;
    end
end

if ~isempty(bounds) && (L - bounds(end) < minLen)
    bounds(end) = [];
end
end

function segments = build_segments(L, bounds, minLen)
if isempty(bounds)
    segments = {1:L};
else
    segments = {};
    prev = 1;
    for b = bounds
        segments{end+1} = prev:b; %#ok<AGROW>
        prev = b+1;
    end
    segments{end+1} = prev:L;
end
segments = segments(cellfun(@numel,segments) >= minLen);
end

