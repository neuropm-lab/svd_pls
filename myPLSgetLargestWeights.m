function out = myPLSgetLargestWeights(plsout, comp, pct)
% function to extract top n-th percent weights from plsout structure
% INPUT: plsout = output structure from PLS_corr
%        comp   = integer of component to extract weights from
%        pct    = percentage of weights to keep specified from [0 1]
%
% OUTPUT: out = structure with outputfields similar to plsout.boot but now
%               only the top n-th percent of weights are included


fNames = fieldnames(plsout.boot);
fieldsToExtract = {'Ubr', 'Ubmean', 'Vbr', 'Vbmean'};

fIdx = find(ismember(fNames, fieldsToExtract));

for i=1:length(fIdx)
    
    % get weights of all structure fields for desired component
    w = plsout.boot.(fNames{fIdx(i)})(:,comp);
    % get cutoff for n-th percentage of weights
    nFeat = round(length(w)*pct);

    % sort weights
    [w_sort, sortIdx] = sort(abs(w), 'descend');
    % keep only the largests n-th percent
%     out.([fNames{fIdx(i)} '_top']) = w_sort(1:nFeat);
    out.([fNames{fIdx(i)} '_topIdx']) = sortIdx(1:nFeat);
    out.([fNames{fIdx(i)} '_top'])    = w(sortIdx(1:nFeat));
        
    % if Ubmean or Vbmean also extract the corresponding CI
    if ~isempty(strmatch(fNames{fIdx(i)}, 'Ubmean')) 
        tmpCI = squeeze(plsout.boot.Uci(sortIdx(1:nFeat),comp,:));
        % in cases where only one feature is displayed, make sure that CI
        % is column vector
        if size(tmpCI,2)==1
            tmpCI = tmpCI';
        end
        out.([fNames{fIdx(i)} '_topCI']) = tmpCI;
%         out.([fNames{fIdx(i)} '_topCI']) = squeeze(plsout.boot.Uci(sortIdx(1:nFeat),comp,:));
    end
    
    if ~isempty(strmatch(fNames{fIdx(i)}, 'Vbmean')) 
        clear tmpCI
        tmpCI = squeeze(plsout.boot.Vci(sortIdx(1:nFeat),comp,:));
        % in cases where only one feature is displayed, make sure that CI
        % is column vector
        if size(tmpCI,2)==1
            tmpCI = tmpCI';
        end
        out.([fNames{fIdx(i)} '_topCI']) = tmpCI;
%         out.([fNames{fIdx(i)} '_topCI']) = squeeze(plsout.boot.Vci(sortIdx(1:nFeat),comp,:));
    end
    
end


