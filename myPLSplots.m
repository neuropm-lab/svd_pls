function myPLSplots(PLSout, plotComp, pct, Xlabels, Ylabels)
%
% function to generate plot based on PLScorr script
% USAGE: myPLSplots(PLSout, plotComp, pct, Xlabels, Ylabels)
%
% INPUT: PLSout   = output structure from PLScorr.m
%        plotComp = significant component to plot
%                   (default = 1)
%        pct      = optional - top n-th percentage of weights to plot
%                   indicated in [0 1] range, (default - top 20%)
%        Xlabels  = optional - cell with labels of feautres in PLSout.plotVars.X
%        Ylabels  = optional - cell with labels of feautres in PLSout.plotVars.Y
%                   if no labels provided, features will be numbered by
%                   default
%
% generates the following plots:
% - plot with PLSout.plotVars.Xplained Variance per LV and associated
%   p-values
% - plots for correlation of LVs, bootstrap ratios, as well as averaged 
%   bootstrap weights + CI for X and Y features
%
% TRB, NeuroPM lab, MNI, August2020

cblue  = [83, 125, 255]./255;
cred   = [255, 84, 83]./255;
cgreen = [44, 200, 77]./255;
cgrey  = [0.3 0.3 0.3];


if ~exist('Xlabels', 'var') || isempty(Xlabels)
    Xlabels = [1:1:PLSout.plotVars.nXfeat];
end

if ~exist('Ylabels', 'var') || isempty(Ylabels)
    Ylabels = [1:1:PLSout.plotVars.nYfeat];
end

% which component to plot?
if ~exist('plotComp', 'var') || isempty(plotComp)
    plotComp = 1;
end

% plot only features of X and Y above threshold
if ~exist('pct', 'var') 
    pct = 0.2;
end

% reference line for bootstrapping ratios at the following threshold
thr = 2.3;


%===================== Explained variance and p-values ====================
nComps = length(PLSout.explVarLVs);

figure,
yyaxis left
plot(1:nComps, PLSout.explVarLVs*100, 'o-', 'markerfacecolor', 'w', 'markersize', 10, 'color', cgrey, 'linewidth', 2)
ylabel('explained covariance [%]')
yyaxis right
plot(1:nComps, PLSout.perm.myLVpvals, 's-', 'markerfacecolor', 'w', 'markersize', 10, 'color', cgreen, 'linewidth', 2)
grid on
ylabel('p-value')
xlabel('component')
href = refline(0, 0.05);
href.Color = cred;
xlim([0 nComps])

ax = gca;
ax.YAxis(1).Color = cgrey;
ax.YAxis(2).Color = cgreen;


% ================= plot top n-th percentage of weights ===================
out = myPLSgetLargestWeights(PLSout, plotComp, pct);

% ==== plot bootstrap ratios Ubr
plotIn  = out.Ubr_top;
plotIdx = out.Ubr_topIdx;

figure,
bar(plotIn, 'facecolor', cblue)
grid on

% appearance
xlim([0 length(plotIn)+1])
xlabel('X variables');
ylabel('X saliences');
title(['LV' num2str(plotComp) ' - X bootstrap ratios']);    
set(gca, 'XTick', [1:1:length(plotIn)], 'XTickLabel', Xlabels(plotIdx), 'XTickLabelRotation', 30, 'TickLabelInterpreter','none')
set(gca, 'fontsize', 10)


% ==== plot bootsrap means + 95% CI for U
plotIn  = out.Ubmean_top;
plotIdx = out.Ubmean_topIdx;
plotCI  = out.Ubmean_topCI;

figure,
bar(plotIn, 'facecolor', cblue)
hold on
errorbar(1:length(plotIn), plotIn, plotIn-plotCI(:,1), plotIn-plotCI(:,2), '.', 'color', cgrey)
grid on

% appearance
xlim([0 length(plotIn)+1])
xlabel('X variables');
ylabel('X saliences');
title(['LV' num2str(plotComp) ' - X bootstrap weights + 95% CI']);    
set(gca, 'XTick', [1:1:length(plotIn)], 'XTickLabel', Xlabels(plotIdx), 'XTickLabelRotation', 30, 'TickLabelInterpreter','none')
set(gca, 'fontsize', 10)




% ==== plot bootstrap ratios Vbr
plotIn  = out.Vbr_top;
plotIdx = out.Vbr_topIdx;

figure,
bar(plotIn, 'facecolor', cred)
grid on

% appearance
xlim([0 length(plotIn)+1])
xlabel('Y variables');
ylabel('Y saliences');
title(['LV' num2str(plotComp) ' - Y bootstrap ratios']);    
set(gca, 'XTick', [1:1:length(plotIn)], 'XTickLabel', Ylabels(plotIdx), 'XTickLabelRotation', 30, 'TickLabelInterpreter','none')
set(gca, 'fontsize', 10)


% ==== plot bootsrap means + 95% CI for V
plotIn  = out.Vbmean_top;
plotIdx = out.Vbmean_topIdx;
plotCI  = out.Vbmean_topCI;

figure,
bar(plotIn, 'facecolor', cred)
hold on
errorbar(1:length(plotIn), plotIn, plotIn-plotCI(:,1), plotIn-plotCI(:,2), '.', 'color', cgrey)
grid on

% appearance
xlim([0 length(plotIn)+1])
xlabel('Y variables');
ylabel('Y saliences');
title(['LV' num2str(plotComp) ' - Y bootstrap weights + 95% CI']);    
set(gca, 'XTick', [1:1:length(plotIn)], 'XTickLabel', Ylabels(plotIdx), 'XTickLabelRotation', 30, 'TickLabelInterpreter','none')
set(gca, 'fontsize', 10)



% ======================= plot correlation of LVs =========================

for iter_lv=plotComp   %1:PLSout.numSignifLVs
    this_lv = PLSout.perm.mySignifLVs(iter_lv);
    
    % plot correlation of current LV
    figure,
    plot(PLSout.Lx(:, this_lv), PLSout.Ly(:, this_lv), 'o', 'color', [0.2 0.2 0.2], 'markerfacecolor', 'w', 'markersize', 10);
    
    % appearance
    xlim([min(PLSout.Lx(:, this_lv))-1 max(PLSout.Lx(:, this_lv))+1])
    ylim([min(PLSout.Ly(:, this_lv))-1 max(PLSout.Ly(:, this_lv))+1])
    lsline
    grid on
    
    xlabel(['X scores LV' num2str(this_lv)]);
    ylabel(['Y scores LV' num2str(this_lv)]);

    [r,p] = corr(PLSout.Lx(:,this_lv),PLSout.Ly(:,this_lv));
    fprintf('explained covariance of LV%d: %.2f%%, p-permutation: %.4f, correlation of LVs r:%.2f, p:%.4f \n Average bootstrap covariance:%.2f \n', ...
        this_lv, PLSout.explVarLVs(this_lv)*100, PLSout.perm.myLVpvals(this_lv), r, p, mean(PLSout.boot.explVarLVs(this_lv,:),2)*100)
    %disp(['r = ' num2str(r,'%0.2f') ' p =' num2str(p,'%0.5f')]);
    title( sprintf('Covar:%.2f%%, p-perm:%.4f, r:%.2f, p:%.4f', PLSout.explVarLVs(this_lv)*100, PLSout.perm.myLVpvals(this_lv), r, p) )
end
    



end
