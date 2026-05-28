%% Conditional Saddle Path Figures — Bianchi (2011)
% 'Dancing on the Saddles: A Geometric Framework for Stochastic Equilibrium Dynamics' (Lee, 2026)

% 1. Setup Data from RTM Workspace
valid = BURNIN:(pathlength-1);

% Define the two extreme shocks to plot
% Ny=4 means 16 total states. 
% State 1 = Lowest (yT, yN), State 16 = Highest (yT, yN)
shocks_to_plot = [2, Ny*Ny-1]; 
shock_names = {'Low Shock (Crisis)', 'High Shock (Boom)'};

% Style Definitions
colorRec    = [0.80, 0.10, 0.10];       % red   (Crisis)
colorExp    = [0.27, 0.50, 0.71];       % blue  (Boom)
colorArrow  = [0, 0, 0];                % black arrows
markerSize  = 18;
dotSize     = 10;
fontSize    = 14;
fontSizeTick = 11;
nBins       = 150;

db = struct(); % Data structure to hold smoothed paths

%% 2. Build Data and Compute Conditional Steady States (CSS)
for i = 1:2
    target_state = shocks_to_plot(i);
    
    % Find indices where this shock occurred
    idx = valid(vStatePath(valid) == target_state);
    
    % Extract raw scatter data
    B_raw  = vB(idx);
    C_raw  = vC_T(idx);
    Bp_raw = vB(idx + 1); % b_{t+1}
    
    db(i).B_scatter = B_raw;
    db(i).C_scatter = C_raw;
    
    % Bin-average for smooth interpolation
    B_lo = min(B_raw);  B_hi = max(B_raw);
    
    % Safety: If min == max, we don't have a spread. Nudge it.
    if B_hi - B_lo < 1e-6
        B_lo = B_lo - 1e-4; B_hi = B_hi + 1e-4;
    end
    
    edges   = linspace(B_lo, B_hi, nBins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    C_bin  = nan(nBins, 1);
    Bp_bin = nan(nBins, 1);
    
    for ib = 1:nBins
        % Use <= on the right side to catch edge cases
        in_bin = (B_raw >= edges(ib)) & (B_raw <= edges(ib+1));
        
        
        if sum(in_bin) >= 1 
            C_bin(ib)  = mean(C_raw(in_bin));
            Bp_bin(ib) = mean(Bp_raw(in_bin));
        end
    end
    
    keep = ~isnan(C_bin);
    B_smooth  = centers(keep)';
    C_smooth  = C_bin(keep);
    Bp_smooth = Bp_bin(keep);
    
    % 
    if length(B_smooth) < 2
        fprintf('Warning: Not enough spread in data for %s. Using simple mean for CSS.\n', shock_names{i});
        Bcss = mean(B_raw);
        Ccss = mean(C_raw);
    else
        % Compute CSS: iterate b'(b) forward from a guess
        B_curr = mean(B_raw); 
        for ii = 1:5000
            B_curr = max(min(B_curr, B_hi), B_lo); 
            B_next = interp1(B_smooth, Bp_smooth, B_curr, 'linear', 'extrap');
            
            if abs(B_next - B_curr) < 1e-6
                break; 
            end
            B_curr = B_next;
        end
        Bcss = B_curr;
        Ccss = interp1(B_smooth, C_smooth, Bcss, 'linear', 'extrap');
    end
    
    db(i).B_smooth  = B_smooth(:);
    db(i).C_smooth  = C_smooth(:);
    db(i).Bp_smooth = Bp_smooth(:);
    db(i).Bcss = Bcss;
    db(i).Ccss = Ccss;
    
    fprintf('CSS %s: b_css = %.4f, c_css = %.4f\n', shock_names{i}, Bcss, Ccss);
end

%% 3. Plot the Saddle Paths
Bcss_rec = db(1).Bcss;  Ccss_rec = db(1).Ccss;
Bcss_exp = db(2).Bcss;  Ccss_exp = db(2).Ccss;

all_B = [db(1).B_scatter; db(2).B_scatter; Bcss_rec; Bcss_exp];
all_C = [db(1).C_scatter; db(2).C_scatter; Ccss_rec; Ccss_exp];
xpad = 0.05 * (max(all_B) - min(all_B));
ypad = 0.05 * (max(all_C) - min(all_C));

fig = figure('Units','inches','Position',[0 0 8 6], 'Name', 'Bianchi Saddle Paths');
hold on;

% Scatter dots
scatter(db(1).B_scatter, db(1).C_scatter, dotSize, ...
    colorRec, 'filled', 'MarkerFaceAlpha', 0.4);
scatter(db(2).B_scatter, db(2).C_scatter, dotSize, ...
    colorExp, 'filled', 'MarkerFaceAlpha', 0.4);

% CSS pentagram markers
plot(Bcss_rec, Ccss_rec, 'pentagram', 'MarkerSize', markerSize, ...
    'MarkerFaceColor', colorRec, 'MarkerEdgeColor', 'none');
plot(Bcss_exp, Ccss_exp, 'pentagram', 'MarkerSize', markerSize, ...
    'MarkerFaceColor', colorExp, 'MarkerEdgeColor', 'none');

% Vertical dashed lines at CSS
xline(Bcss_rec, '--', 'Color', colorRec, 'LineWidth', 0.8, 'HandleVisibility', 'off');
xline(Bcss_exp, '--', 'Color', colorExp, 'LineWidth', 0.8, 'HandleVisibility', 'off');

% Set axis limits before drawing arrows (annotation needs final axes)
xlim([min(all_B) - xpad, max(all_B) + xpad]);
ylim([min(all_C) - ypad, max(all_C) + ypad]);

% Convergence arrows (SaddleDance styling)
% -1 puts arrows below the curve, +1 puts them above
addSaddleArrows(gca, db(1).B_smooth, db(1).C_smooth, Bcss_rec, colorArrow, -1);
addSaddleArrows(gca, db(2).B_smooth, db(2).C_smooth, Bcss_exp, colorArrow, +1);

% Formatting
xlabel('Current Debt ($b_t$)', 'FontSize', fontSize, 'Interpreter', 'latex');
ylabel('Tradable Consumption ($c^T_t$)', 'FontSize', fontSize, 'Interpreter', 'latex');
title('Bianchi (2011) - Conditional Saddle Paths', 'FontSize', fontSize+1);
legend('Low Shock (Crisis)', 'High Shock (Boom)', ...
       'CSS: Low', 'CSS: High', ...
       'Location', 'northwest', 'FontSize', 10, 'Interpreter', 'latex');
set(gca, 'FontSize', fontSizeTick, 'TickDir', 'out');
box(gca, 'off');
grid on;
hold off;


%% ========================================================================
% LOCAL FUNCTIONS FOR SADDLEDANCE ARROWS
% =========================================================================

function addSaddleArrows(ax, Nv, Cv, Ncss, colArrow, side)
    nPts = length(Nv);
    if nPts < 6, return; end

    xl = get(ax, 'XLim');
    yl = get(ax, 'YLim');
    rx = xl(2) - xl(1);
    ry = yl(2) - yl(1);

    % resample to uniform x-spacing
    xq = linspace(min(Nv), max(Nv), 200)';
    yq = interp1(Nv, Cv, xq, 'linear');
    vis = yq >= yl(1) & yq <= yl(2) & xq >= xl(1) & xq <= xl(2);
    xq = xq(vis);  yq = yq(vis);
    if length(xq) < 6, return; end

    % split into left arm (x < CSS) and right arm (x > CSS)
    leftIdx  = find(xq < Ncss);
    rightIdx = find(xq > Ncss);

    hasLeft  = length(leftIdx)  >= 3;
    hasRight = length(rightIdx) >= 3;

    if hasLeft && hasRight
        placeSaddleArrow(ax, xq, yq, leftIdx,  Ncss, +1, rx, ry, xl, yl, colArrow, side, 0.10);
        placeSaddleArrow(ax, xq, yq, rightIdx, Ncss, -1, rx, ry, xl, yl, colArrow, side, 0.10);
    elseif hasLeft && ~hasRight
        placeSaddleArrow(ax, xq, yq, leftIdx, Ncss, +1, rx, ry, xl, yl, colArrow, side, 0.35);
        placeSaddleArrow(ax, xq, yq, leftIdx, Ncss, +1, rx, ry, xl, yl, colArrow, side, 0.10);
    elseif ~hasLeft && hasRight
        placeSaddleArrow(ax, xq, yq, rightIdx, Ncss, -1, rx, ry, xl, yl, colArrow, side, 0.35);
        placeSaddleArrow(ax, xq, yq, rightIdx, Ncss, -1, rx, ry, xl, yl, colArrow, side, 0.10);
    end
end

function placeSaddleArrow(ax, xq, yq, armIdx, xss, dir, rx, ry, xl, yl, colArrow, side, distFrac)
    arrFrac  = 0.06;    
    offFrac  = 0.03;    

    targetX = xss - dir * distFrac * rx;
    [~, bestIdx] = min(abs(xq(armIdx) - targetX));
    ii = armIdx(max(2, min(bestIdx, length(armIdx)-1)));

    span = max(1, round(length(armIdx) * 0.15));
    i_lo = max(armIdx(1), ii - span);
    i_hi = min(armIdx(end), ii + span);
    slope = (yq(i_hi) - yq(i_lo)) / (xq(i_hi) - xq(i_lo) + 1e-15);

    sn = slope * rx / ry;
    mag = sqrt(1 + sn^2);
    dxn = dir / mag;
    dyn = dir * sn / mag;

    ox = -dyn * rx * offFrac * side;
    oy =  dxn * ry * offFrac * side;

    cx = xq(ii) + ox;
    cy = yq(ii) + oy;

    tail_x = cx - dxn * rx * arrFrac / 2;
    tail_y = cy - dyn * ry * arrFrac / 2;
    tip_x  = cx + dxn * rx * arrFrac / 2;
    tip_y  = cy + dyn * ry * arrFrac / 2;

    tipGap = 0.015 * rx;
    if dir == +1 && tip_x > xss - tipGap
        tip_x = xss - tipGap;
    elseif dir == -1 && tip_x < xss + tipGap
        tip_x = xss + tipGap;
    end

    ax_pos = get(ax, 'Position');
    tx_fig = ax_pos(1) + (tail_x - xl(1)) / rx * ax_pos(3);
    ty_fig = ax_pos(2) + (tail_y - yl(1)) / ry * ax_pos(4);
    hx_fig = ax_pos(1) + (tip_x  - xl(1)) / rx * ax_pos(3);
    hy_fig = ax_pos(2) + (tip_y  - yl(1)) / ry * ax_pos(4);

    annotation('arrow', [tx_fig, hx_fig], [ty_fig, hy_fig], ...
        'Color', colArrow, 'LineWidth', 1.8, ...
        'HeadWidth', 10, 'HeadLength', 8, ...
        'HeadStyle', 'vback2');
end