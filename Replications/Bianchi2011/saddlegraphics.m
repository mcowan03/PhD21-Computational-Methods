%% BIANCHI (2011) - CE vs SPP SADDLE PATH OVERLAY
clear; clc; close all;

% 1. Load Data
ce  = load('prelimresults.mat');
spp = load('prelimresultsSPP.mat');

% 2. Setup Shocks to Plot 
% Using 1 for extreme crisis, and Ny*Ny for extreme boom
Ny = ce.Ny;
shocks_to_plot = [2, Ny*Ny-1]; 
nBins = 150;

% Process Data using helper function (defined at bottom)
db_ce_crisis = get_saddle_data(ce, shocks_to_plot(1), nBins, 'CE Low');
db_ce_boom   = get_saddle_data(ce, shocks_to_plot(2), nBins, 'CE High');
db_spp_crisis = get_saddle_data(spp, shocks_to_plot(1), nBins, 'SPP Low');
db_spp_boom   = get_saddle_data(spp, shocks_to_plot(2), nBins, 'SPP High');

% 3. Style Definitions
col_CE_Crisis  = [0.80, 0.10, 0.10]; % Red
col_CE_Boom    = [0.10, 0.30, 0.80]; % Dark Blue
col_SPP_Crisis = [0.90, 0.50, 0.10]; % Orange
col_SPP_Boom   = [0.20, 0.70, 0.80]; % Teal/Cyan
col_Arrow      = [0, 0, 0];          % Black arrows

markerSize = 16;
dotSize    = 8;
alphaLevel = 0.5; % Lowered alpha because 4 clouds will overlap
fontSize   = 12;

%% 4. Plot the Overlay
fig = figure('Units','inches','Position',[1 1 10 7], 'Name', 'CE vs SPP Saddle Paths');
hold on; grid on;

% A. Scatter dots (The RTM Clouds)
scatter(db_ce_crisis.B_scatter, db_ce_crisis.C_scatter, dotSize, col_CE_Crisis, 'filled', 'MarkerFaceAlpha', alphaLevel);
scatter(db_spp_crisis.B_scatter, db_spp_crisis.C_scatter, dotSize, col_SPP_Crisis, 'filled', 'MarkerFaceAlpha', alphaLevel);
scatter(db_ce_boom.B_scatter, db_ce_boom.C_scatter, dotSize, col_CE_Boom, 'filled', 'MarkerFaceAlpha', alphaLevel);
scatter(db_spp_boom.B_scatter, db_spp_boom.C_scatter, dotSize, col_SPP_Boom, 'filled', 'MarkerFaceAlpha', alphaLevel);

% B. CSS pentagram markers
h1 = plot(db_ce_crisis.Bcss, db_ce_crisis.Ccss, 'pentagram', 'MarkerSize', markerSize, 'MarkerFaceColor', col_CE_Crisis, 'MarkerEdgeColor', 'k');
h2 = plot(db_spp_crisis.Bcss, db_spp_crisis.Ccss, 'pentagram', 'MarkerSize', markerSize, 'MarkerFaceColor', col_SPP_Crisis, 'MarkerEdgeColor', 'k');
h3 = plot(db_ce_boom.Bcss, db_ce_boom.Ccss, 'pentagram', 'MarkerSize', markerSize, 'MarkerFaceColor', col_CE_Boom, 'MarkerEdgeColor', 'k');
h4 = plot(db_spp_boom.Bcss, db_spp_boom.Ccss, 'pentagram', 'MarkerSize', markerSize, 'MarkerFaceColor', col_SPP_Boom, 'MarkerEdgeColor', 'k');

% C. Vertical dashed lines at CSS
xline(db_ce_crisis.Bcss, '--', 'Color', col_CE_Crisis, 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(db_spp_crisis.Bcss, '--', 'Color', col_SPP_Crisis, 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(db_ce_boom.Bcss, '--', 'Color', col_CE_Boom, 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(db_spp_boom.Bcss, '--', 'Color', col_SPP_Boom, 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Axis Limits
all_B = [db_ce_crisis.B_scatter; db_ce_boom.B_scatter; db_spp_crisis.B_scatter; db_spp_boom.B_scatter];
all_C = [db_ce_crisis.C_scatter; db_ce_boom.C_scatter; db_spp_crisis.C_scatter; db_spp_boom.C_scatter];
xpad = 0.05 * (max(all_B) - min(all_B));
ypad = 0.05 * (max(all_C) - min(all_C));
xlim([min(all_B) - xpad, max(all_B) + xpad]);
ylim([min(all_C) - ypad, max(all_C) + ypad]);

% D. Convergence arrows
% Offsetting the 'side' parameter so arrows don't overlap between CE and SPP
addSaddleArrows(gca, db_ce_crisis.B_smooth, db_ce_crisis.C_smooth, db_ce_crisis.Bcss, col_Arrow, -1);
addSaddleArrows(gca, db_spp_crisis.B_smooth, db_spp_crisis.C_smooth, db_spp_crisis.Bcss, col_Arrow, +1);
addSaddleArrows(gca, db_ce_boom.B_smooth, db_ce_boom.C_smooth, db_ce_boom.Bcss, col_Arrow, -1);
addSaddleArrows(gca, db_spp_boom.B_smooth, db_spp_boom.C_smooth, db_spp_boom.Bcss, col_Arrow, +1);

% E. Formatting
xlabel('Current Debt (b_t)', 'FontSize', fontSize);
ylabel('Tradable Consumption (c^T_t)', 'FontSize', fontSize);
title('Bianchi (2011) - CE vs SPP Conditional Saddle Paths', 'FontSize', fontSize+2);

% Clean Legend (Only showing the CSS markers to prevent clutter)
legend([h1, h2, h3, h4], {'CE: Low Shock', 'SPP: Low Shock', 'CE: High Shock', 'SPP: High Shock'}, ...
       'Location', 'northwest', 'FontSize', 11);
set(gca, 'TickDir', 'out');
box(gca, 'off');
hold off;


%% ========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================
function out = get_saddle_data(data, target_state, nBins, name_label)
    % Extracts data and calculates CSS for a specific struct and shock
    valid = data.BURNIN:(data.pathlength-1);
    idx = valid(data.vStatePath(valid) == target_state);
    
    B_raw  = data.vB(idx);
    C_raw  = data.vC_T(idx);
    Bp_raw = data.vB(idx + 1);
    
    out.B_scatter = B_raw;
    out.C_scatter = C_raw;
    
    B_lo = min(B_raw);  B_hi = max(B_raw);
    if B_hi - B_lo < 1e-6
        B_lo = B_lo - 1e-4; B_hi = B_hi + 1e-4;
    end
    
    edges   = linspace(B_lo, B_hi, nBins+1);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    C_bin  = nan(nBins, 1);
    Bp_bin = nan(nBins, 1);
    
    for ib = 1:nBins
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
    
    if length(B_smooth) < 2
        Bcss = mean(B_raw);
        Ccss = mean(C_raw);
    else
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
    
    out.B_smooth = B_smooth(:);
    out.C_smooth = C_smooth(:);
    out.Bcss = Bcss;
    out.Ccss = Ccss;
end

function addSaddleArrows(ax, Nv, Cv, Ncss, colArrow, side)
    nPts = length(Nv);
    if nPts < 6, return; end
    xl = get(ax, 'XLim');
    yl = get(ax, 'YLim');
    rx = xl(2) - xl(1);
    ry = yl(2) - yl(1);
    xq = linspace(min(Nv), max(Nv), 200)';
    yq = interp1(Nv, Cv, xq, 'linear');
    vis = yq >= yl(1) & yq <= yl(2) & xq >= xl(1) & xq <= xl(2);
    xq = xq(vis);  yq = yq(vis);
    if length(xq) < 6, return; end
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