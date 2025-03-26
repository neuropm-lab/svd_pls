function myPLSplots2(PLSout, plotComp, pctX, pctY, Xlabels, Ylabels)
%
% USAGE:
%   myPLSplots2(PLSout, plotComp, pctX, pctY, Xlabels, Ylabels)
%
% DESCRIPTION:
%   This function generates a suite of high-quality visualizations for interpreting
%   Partial Least Squares (PLS) correlation analysis results, particularly from
%   the NeuroPM Lab's PLScorr pipeline. It includes plots of explained variance,
%   permutation-derived p-values, saliences (bootstrap ratios and weights with 95% CIs),
%   and LV score correlations.
%   the old plotting function had a logic error where the upper bar was plotted as
%   plotIn - plotCI(:,2), but the correct upper bar should be upper_bound - plotIn. Corrected this?
%
% INPUTS:
%   PLSout    - Structure output from PLScorr.m containing LV scores, weights,
%               bootstrap estimates, permutation results, and explained variance.
%   plotComp  - Integer specifying the latent variable (LV) component to visualize (default = 1).
%   pctX      - Proportion (0–1) of top-weighted X features to plot (default = 0.2).
%   pctY      - Proportion (0–1) of top-weighted Y features to plot (default = 0.2).
%   Xlabels   - Cell array of labels for X features (optional). If omitted, features are indexed numerically.
%   Ylabels   - Cell array of labels for Y features (optional). If omitted, features are indexed numerically.
%
% OUTPUTS:
%   Generates multiple figures:
%     - LV explained variance and permutation p-values
%     - Top X and Y bootstrap ratios
%     - X and Y bootstrap means with 95% confidence intervals
%     - Scatter plot of X and Y LV scores with correlation and significance annotation
%
% AUTHOR:
%   Adapted by Joon Hwan Hong, based on original by TRB (NeuroPM Lab, MNI, August 2020)

% Define color palette for consistent visualization
cblue  = [83, 125, 255]./255;
cred   = [255, 84, 83]./255;
cgreen = [44, 200, 77]./255;
cgrey  = [0.3 0.3 0.3];

% Set default X and Y labels if not provided
if ~exist('Xlabels', 'var') || isempty(Xlabels)
    Xlabels = arrayfun(@num2str, 1:PLSout.plotVars.nXfeat, 'UniformOutput', false);
end
if ~exist('Ylabels', 'var') || isempty(Ylabels)
    Ylabels = arrayfun(@num2str, 1:PLSout.plotVars.nYfeat, 'UniformOutput', false);
end

% Set default latent variable component to plot
if ~exist('plotComp', 'var') || isempty(plotComp)
    plotComp = 1;
end

% Set default proportion of top features to visualize
if ~exist('pctX', 'var')
    pctX = 0.2;
end
if ~exist('pctY', 'var')
    pctY = 0.2;
end

% Determine number of latent variable components
nComps = length(PLSout.explVarLVs);

% Plot explained variance and p-values for all LVs
fig1 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
ax1 = gca;
yyaxis left
plot(1:nComps, PLSout.explVarLVs*100, 'o-', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'Color', cgrey, 'LineWidth', 2)
ylabel('Explained Covariance (%)', 'FontSize', 12, 'FontName', 'DejaVu Sans')

yyaxis right
plot(1:nComps, PLSout.perm.myLVpvals, 's-', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'Color', cgreen, 'LineWidth', 2)
ylabel('p-value', 'FontSize', 12, 'FontName', 'DejaVu Sans')

% Draw horizontal reference line at p = 0.05
href = refline(0, 0.05);
href.Color = cred;
href.LineWidth = 1.5;
xlim([0 nComps])

xlabel('Component', 'FontSize', 12, 'FontName', 'DejaVu Sans')
grid on
set(ax1, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

% Extract top-weighted features for X and Y
outX = myPLSgetLargestWeights(PLSout, plotComp, pctX);
outY = myPLSgetLargestWeights(PLSout, plotComp, pctY);

% Plot bootstrap ratios for X features
fig2 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
bar(outX.Ubr_top, 'FaceColor', cblue)
grid on
xlabel('X variables', 'FontSize', 12, 'FontName', 'DejaVu Sans');
ylabel('X bootstrap ratios', 'FontSize', 12, 'FontName', 'DejaVu Sans');
title(['LV' num2str(plotComp) ' - X Bootstrap Ratios'], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'DejaVu Sans');
set(gca, 'XTick', 1:length(outX.Ubr_top), 'XTickLabel', Xlabels(outX.Ubr_topIdx),...
    'XTickLabelRotation', 30, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

% Plot bootstrap means and 95% CI for X features
fig3 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
bar(outX.Ubmean_top, 'FaceColor', cblue); hold on;
X_lower_err = outX.Ubmean_top - outX.Ubmean_topCI(:,1);
X_upper_err = outX.Ubmean_topCI(:,2) - outX.Ubmean_top;
errorbar(1:length(outX.Ubmean_top), outX.Ubmean_top, X_lower_err, X_upper_err, '.',...
    'Color', cgrey, 'CapSize', 6, 'LineWidth', 1.5);
grid on
xlabel('X variables', 'FontSize', 12, 'FontName', 'DejaVu Sans');
ylabel('X bootstrap weights', 'FontSize', 12, 'FontName', 'DejaVu Sans');
title(['LV' num2str(plotComp) ' - X Bootstrap Weights + 95% CI'], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'DejaVu Sans');
set(gca, 'XTick', 1:length(outX.Ubmean_top), 'XTickLabel', Xlabels(outX.Ubmean_topIdx),...
    'XTickLabelRotation', 30, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

% Plot bootstrap ratios for Y features
fig4 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
bar(outY.Vbr_top, 'FaceColor', cred)
grid on
xlabel('Y variables', 'FontSize', 12, 'FontName', 'DejaVu Sans');
ylabel('Y bootstrap ratios', 'FontSize', 12, 'FontName', 'DejaVu Sans');
title(['LV' num2str(plotComp) ' - Y Bootstrap Ratios'], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'DejaVu Sans');
set(gca, 'XTick', 1:length(outY.Vbr_top), 'XTickLabel', Ylabels(outY.Vbr_topIdx),...
    'XTickLabelRotation', 30, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

% Plot bootstrap means and 95% CI for Y features
fig5 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
bar(outY.Vbmean_top, 'FaceColor', cred); hold on;
Y_lower_err = outY.Vbmean_top - outY.Vbmean_topCI(:,1);
Y_upper_err = outY.Vbmean_topCI(:,2) - outY.Vbmean_top;
errorbar(1:length(outY.Vbmean_top), outY.Vbmean_top, Y_lower_err, Y_upper_err, '.',...
    'Color', cgrey, 'CapSize', 6, 'LineWidth', 1.5);
grid on
xlabel('Y variables', 'FontSize', 12, 'FontName', 'DejaVu Sans');
ylabel('Y bootstrap weights', 'FontSize', 12, 'FontName', 'DejaVu Sans');
title(['LV' num2str(plotComp) ' - Y Bootstrap Weights + 95% CI'], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'DejaVu Sans');
set(gca, 'XTick', 1:length(outY.Vbmean_top), 'XTickLabel', Ylabels(outY.Vbmean_topIdx),...
    'XTickLabelRotation', 30, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

% Plot correlation of X and Y LV scores
fig6 = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
plot(PLSout.Lx(:, plotComp), PLSout.Ly(:, plotComp), 'o', 'Color', cgrey, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on;
lsline % add least-squares regression line
grid on
xlabel(['X scores LV' num2str(plotComp)], 'FontSize', 12, 'FontName', 'DejaVu Sans');
ylabel(['Y scores LV' num2str(plotComp)], 'FontSize', 12, 'FontName', 'DejaVu Sans');
[r,p] = corr(PLSout.Lx(:,plotComp),PLSout.Ly(:,plotComp));
title(sprintf('Covariance: %.2f%%, p-perm: %.4f, r=%.2f, p=%.4f', PLSout.explVarLVs(plotComp)*100, PLSout.perm.myLVpvals(plotComp), r, p),...
    'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'DejaVu Sans');
set(gca, 'FontSize', 12, 'FontName', 'DejaVu Sans', 'LineWidth', 1, 'TickDir', 'out')

end



