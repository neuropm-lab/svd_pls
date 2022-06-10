function PLSout = PLScorr(X, Y, opts)
%
% USAGE: PLSout = PLScorr(X, Y, opts)
%        function to perform SVD on covariance matrix between X and Y,
%        provides options for data normalisation (z-score or boxcox),
%        includes permutation tests to assess significance of components,
%        as well as boostrapping to examine robustness of features and
%        includes corrections for axis rotation with procrustes
%        transformation
%
% INPUT: X    = n x d dimensional matrix (n=subjects, d=features)
%        Y    = n x l dimensional matrix (n=subjects, l=features) 
%        opts = structure with following options (optional)
%         .nPerm   = number of permutations to asses significance 
%                    (default = 500)
%         .nBoot   = number of bootstrap iterations to assess robustness of
%                    features (default = 500)
%         .norm    = normalisation option with strings 'zscore' or 'boxcox'
%                    (default = 'zscore')
%         .Xremove = covariables to remove from X
%         .Yremove = covariables to remove from Y
%
% OUTPUT: plsout = structure containing the following fields
%          .Lx, Ly     = scores of X and Y (transformed data) 
%                        ( n x min(d,l) )
%          .U          = saliences, or weights, matrix for features in X 
%                        ( d x min(d,l) )
%          .V          = saliences, or weights, matrix for features in Y
%                        ( l x min(d,l) )
%          .loadLxX    = loadings for features in X
%                        (correlations of original features in X with Lx)
%          .loadLyY    = loadings for features in Y
%                        (correlations of original features in Y with Ly)
%          .loadLxY    = cross-loadings for features in X
%                        (correlations between features in Y and Lx)
%          .loadLyX    = cross-loadings for features in Y
%                        (correlations between features in X and Ly)
%          .explVarLVs = 1 x min(d,l) vector with explained covariance per
%                        LV
%          .perm          = sub-structure for permutation related outputs
%           .myLVpvals    = 1 x min(d,l) vector of p-values per LV
%           .mySignifLVs  = indices of significant LVs
%           .numSignifLVs = number of significant LVs
%           .Sp           = singular values for each permutation iteration
%           .explVarLVs   = explained covariance per permutation iteration
%          .boot        = sub-structure with bootstrap related outputs
%           .Ubr        = bootstrap ratios for features in X
%           .Vbr        = bootstrap ratios for features in Y
%           .Uci        = 95% confidence intervals for weights in X
%           .Vci        = 95% confidence intervals for weights in Y
%           .Ubmean     = bootstrapped average weights for X
%           .Vbmean     = bootstrapped average weights for Y
%           .explVarLVs = explained covariance per bootstrap iteration
%
%
% based on functions from Dimitri van de Ville, and Rotman Baycrest
%
% TRB, NeuroPM lab, MNI, August2020


% ========================= check inputs ==================================
% get input sizes
global nSub nXfeat nYfeat
nSub   = size(X,1); % number of subjects
nXfeat = size(X,2); % number of features in X
nYfeat = size(Y,2); % number of features in Y


if size(X,1) ~= size(Y,1)
    error('Error, inputs must have same number of rows!')
end

if any(isnan(X(:))) || any(isinf(X(:))) == 1
    error('Matrix X has NaN or Inf')
end

if any(isnan(Y(:))) || any(isinf(Y(:))) == 1
    error('Matrix Y has NaN or Inf')
end


% declare empty options structure if no options provided
if nargin<3
    opts=[];
end
    
    % assign default values for nPerm and nBoot 
    if ~isfield(opts, 'nPerm') || isempty(opts.nPerm)
        opts.nPerm=500;
    end
    if ~isfield(opts, 'nBoot') || isempty(opts.nBoot)
        opts.nBoot=500;
    end

    % assign default value for normalisation option if not provided
    if ~isfield(opts, 'norm') || isempty(opts.norm)
        opts.norm = 'zscore';
        opts.FlagNorm = 1;
    else
        % check if zscore or boxcox (1=zscore, 2=boxcox+zscore)
        if strcmp(opts.norm, 'zscore')
            opts.FlagNorm = 1;
        else
            opts.FlagNorm = 2;
        end

    end

    % check if covariables should be removed
    if ~isfield(opts, 'Xremove') || isempty(opts.Xremove)
        opts.FlagXrem = 0;
    else
        opts.FlagXrem = 1;

        % check if sizes match
        if nSub ~= size(opts.Xremove,1)
            error('Error, X and Xremove need to have same number of rows')
        end
        % check for missing data
        if any(isnan(opts.Xremove(:))) || any(isinf(opts.Xremove(:))) == 1
            error('Matrix Xremove has NaN or Inf')
        end
    end

    if ~isfield(opts, 'Yremove') || isempty(opts.Yremove)
        opts.FlagYrem = 0;
    else
        opts.FlagYrem = 1;

        % check if sizes match
        if nSub ~= size(opts.Yremove,1)
            error('Error, Y and Yremove need to have same number of rows')
        end
        % check for missing data
        if any(isnan(opts.Yremove(:))) || any(isinf(opts.Yremove(:))) == 1
            error('Matrix Yremove has NaN or Inf')
        end
    end


% ============================ start analysing ============================

% Save original matrices
X0 = X; 
Y0 = Y;

% adjust and normalise input data
if opts.FlagXrem == 1
    Xa = removeCvar(X, opts.Xremove);
    X  = PLS_norm(Xa, opts.FlagNorm);
else
    X  = PLS_norm(X, opts.FlagNorm);
end
if opts.FlagYrem == 1
    Ya = removeCvar(Y, opts.Yremove);
    Y  = PLS_norm(Ya, opts.FlagNorm);
else
    Y  = PLS_norm(Y, opts.FlagNorm);
end



%========================= Cross-covariance matrix ========================

clear R
% R = X'*Y;
[R, FlagDim] = PLS_cov(X,Y);


%====================== Singular value decomposition ======================

clear U S V

[U,S,V] = svd(R,'econ');
% reassign weights so that U always belongs to X and V always to Y in case
% of transposing R
if FlagDim ==1
    tmp1 = U;
    tmp2 = V;
    clear U V
    U = tmp2;
    V = tmp1;
end
    

% number of latent variables
nLVs = min(size(S));

% Explained covariance by each LV
explVarLVs = (diag(S).^2) / sum(diag(S.^2)); 


% ICA convention: turn latent variables (LVs) such that max is positive
for iter_lv = 1 : nLVs
    [~,idx] = max(abs(V(:,iter_lv)));
    if sign(V(idx,iter_lv)) < 0
        V(:,iter_lv) = -V(:,iter_lv);
        U(:,iter_lv) = -U(:,iter_lv);
    end
end

%======================== PLS Scores and Loadings =========================

[Lx, Ly, Lx_X_corr, Ly_Y_corr, Lx_Y_corr, Ly_X_corr] = ...
            PLS_scores_loadings(X,Y,U,V);
        

%================= Permutation testing for LV significance ================

[Sp, explVarLVPerm] = PLS_perm(X,Y,S,U,V,opts);

% assess significance of LVs with permutations by counting how often the
% singular value of permutation is greater than original singular value
S_mat = repmat(diag(S), 1, opts.nPerm);
sp = sum(Sp >= S_mat, 2);

LVpvals = (sp + 1) ./ (opts.nPerm + 1);

SignifLVs = find(LVpvals<0.05); % index of significant LVs
numSignifLVs = size(SignifLVs,1); % number of significant LVs

% Display significant LVs
disp([num2str(numSignifLVs) ' significant LV(s)']);
for iter_lv = 1:numSignifLVs
    this_lv = SignifLVs(iter_lv);
    disp(['LV' num2str(this_lv) ' - p=' num2str(LVpvals(this_lv),'%0.3f') ]);
end


% ============================== Bootstrapping ============================

if opts.FlagXrem == 1 
    Xb_in = Xa;
else
    Xb_in = X;
end
if opts.FlagYrem == 1
    Yb_in = Ya;
else 
    Yb_in = Y;
end

[Ub_mean, Vb_mean, Ubr, Vbr, Uci, Vci, explLVBoot, Ub_boots, Vb_boots] = PLS_boot(Xb_in,Yb_in,X0,Y0,U,V,S,opts);



%====================== collect output variables ==========================

PLSout.Lx           = Lx;
PLSout.Ly           = Ly;
PLSout.U            = U;
PLSout.V            = V;
PLSout.loadLxX      = Lx_X_corr';
PLSout.loadLyY      = Ly_Y_corr';
PLSout.loadLyX      = Ly_X_corr';
PLSout.loadLxY      = Lx_Y_corr';
PLSout.explVarLVs   = explVarLVs;

PLSout.perm.myLVpvals    = LVpvals;
PLSout.perm.mySignifLVs  = SignifLVs;
PLSout.perm.numSignifLVs = numSignifLVs;
PLSout.perm.Sp           = Sp;
PLSout.perm.explVarLVs   = explVarLVPerm;

PLSout.boot.Ubr        = Ubr;
PLSout.boot.Vbr        = Vbr;
PLSout.boot.Ubmean     = Ub_mean;
PLSout.boot.Vbmean     = Vb_mean;
PLSout.boot.Uci        = Uci;
PLSout.boot.Vci        = Vci;
PLSout.boot.explVarLVs = explLVBoot;
PLSout.boot.Ub_boots   = Ub_boots;
PLSout.boot.Vb_boots   = Vb_boots;

PLSout.plotVars.nSub    = nSub;
PLSout.plotVars.nXfeat  = nXfeat;
PLSout.plotVars.nYfeat  = nYfeat;
PLSout.plotVars.X       = X;
PLSout.plotVars.Y       = Y;
PLSout.plotVars.S       = S;


% shut down parpool
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'))
end

end











% #########################################################################
% =========================== REQUIRED FUNCTIONS ==========================
% #########################################################################

%==========================================================================
function Xadj = removeCvar(X,cvar)
% function to remove effects of covariates from data matrix
% USAGE:  Xadj  = removeCvar(X,cvar)
% INPUT:  X     = input data matrix
%         cvar  = matrix of covariates
% OUTPUT: Xadj  = adjusted input matrix


% check input sizes
if size(X,1) ~= size(cvar,1)
    error('Input matrices need to have same number of rows')
end

% do robust regression, estimating each feature in X with covariates, then
% take Intercept + residuals as covarate adjusted data
% --pre-allocate Xadj
Xadj = zeros(size(X));
for ik=1:size(X,2)
    [b, stats] = robustfit(cvar, X(:,ik));
    Xadj(:,ik) = b(1) + stats.resid;
end 

end
%==========================================================================



%==========================================================================
function X = PLS_norm(X, FlagNorm)
% normalise the data according to normOpt
% FlagNorm = integer specifying if zscore or boxcox normalisation should be
% used

if FlagNorm==1
    
    X = zscore(X);
    
elseif FlagNorm==2
    
    % for large datasets boxcox can be very slow. try using parallel pool
    % to speed up things
    if isempty(gcp('nocreate'))
        parpool(8);
    end
    
    parfor in=1:size(X,2)
        X(:,in) = boxcox(X(:,in) - min(X(:,in)) + eps);
    end
    
    X = zscore(X);
else
    error('Not a valid normalisation option!')
end

end
%==========================================================================



%==========================================================================
function rotatemat=rri_bootprocrust(origlv,bootlv,origlv2,bootlv2)
%syntax rotatemat=rri_bootprocrust(origlv,bootlv)
% - from PLS toolbox

%define coordinate space between original and bootstrap LVs
temp=origlv'*bootlv;

%orthogonalze space
[V,W,U]=svd(temp);

%determine procrustean transform
rotatemat1=U*V';


clear tmp V W U
%define coordinate space between original and bootstrap LVs
temp=origlv2'*bootlv2;
%orthogonalze space
[V,W,U]=svd(temp);
%determine procrustean transform
rotatemat2=U*V';

% average rotatemat
rotatemat = (rotatemat1+rotatemat2)./2;

end
%==========================================================================


%==========================================================================
function [R, FlagDim] = PLS_cov(X,Y)
% function to compute cross-covariance matrix of X and Y. function always
% returns skinny matrix

n1 = size(X,2);
n2 = size(Y,2);

if n2<n1
    
    R = X'*Y;
    
    FlagDim = 0;
else
    
    R = Y'*X;
    
    FlagDim = 1;
end

end
%==========================================================================


%==========================================================================
function [Lx, Ly, Lx_X_corr, Ly_Y_corr, Lx_Y_corr, Ly_X_corr] = ...
            PLS_scores_loadings(X,Y,U,V)
% function to compute the LV scores (projected original features onto new space in which Lx and Ly
% maximally covary for each component. Also computes loadings of original features, similar to loadings in CCA.

% INPUT: X = feature matrix X
%        Y = feature matrix Y
%        U = singular vector matrix U from SVD
%        V = singular vector matrix V from SVD
%
% OUTPUT: Lx = PLS Scores to X (projected original features onto new space)
%         Ly = PLS Scores to Y (projected original features onto new space)
%         Lx_X_corr = PLS Loadings to X (correlation of original feautres
%                     in X with Lx)
%         Ly_Y_corr = PLS Loadings to Y (correlation of original feautres
%                     in Y with Ly)
%         Ly_X_corr = PLS Cross-Loadings to X (correlation of original feautres
%                     in X with Ly)
%         Lx_Y_corr = PLS Cross-Loadings to Y (correlation of original feautres
%                     in Y with Lx)


% compute PLS Scores
Lx = X * U;
Ly = Y * V;

% compute PLS loadings
Lx_X_corr = corr(Lx, X);
Ly_Y_corr = corr(Ly, Y);

% compute PLS cross-loadings
Lx_Y_corr = corr(Lx, Y);
Ly_X_corr = corr(Ly, X);

end
%==========================================================================


%==========================================================================
function [Sp_perm, explVarLVPerm] = PLS_perm(X,Y,S,U,V,opts)
% function to run SVD on permuted data 

disp('... Permutations ...')
rng(1);

for iter_perm = 1:opts.nPerm
    
    % Display number of permutations 
    if mod(iter_perm,100) == 0, disp(num2str(iter_perm)); end
    
    % Leave X unchanged (no need to permute both X and Y matrices)
    Xp = X; % X is already normalized
    
    % Permute Y by shuffling rows (subjects)
    perm_order = randperm(size(X,1));
    Yp = Y(perm_order,:);
         
    % Cross-covariance matrix between X and permuted Y
    [Rp, FlagDim] = PLS_cov(Xp,Yp);
    
    % SVD of Rp
    [Up,Sp,Vp] = svd(Rp,'econ');
    if FlagDim ==1
        tmp1 = Up;
        tmp2 = Vp;
        clear Up Vp
        Up = tmp2;
        Vp = tmp1;
    end
        
    procmat = rri_bootprocrust(V, Vp, U, Up);
    Vp = Vp * Sp * procmat; 
    Sp = sqrt(sum(Vp.^2));
    
    Sp_perm(:,iter_perm) = Sp;
    
    % caculate explained variance for each permuation run
    explVarLVPerm(:,iter_perm) = ((Sp_perm(:,iter_perm)).^2) ./ sum((Sp_perm(:,iter_perm).^2)); 

end

end
%==========================================================================



%========================================================================== 
function [Ub_mean, Vb_mean, Ubr, Vbr, Uci, Vci, explLVBoot, Ub_boots, Vb_boots] = PLS_boot(Xb_in,Yb_in,X0,Y0,U,V,S,opts)

% set random number generator (so results are the same every time)
rng(1);
disp('... Bootstrapping ...');

for iter_boot = 1:opts.nBoot

    % Display number of bootstraps 
    if mod(iter_boot,100) == 0, disp(num2str(iter_boot)); end

    % Get bootstrap subject sampling
    [~, bootIdx] = datasample(Xb_in,size(Xb_in,1));

    % normalise data based on selected samples
    %--Bootstrap of X
    Xb = Xb_in(bootIdx,:);
    Xb = PLS_norm(Xb, opts.FlagNorm);
    %--Bootstrap of Y
    Yb = Yb_in(bootIdx,:);
    Yb = PLS_norm(Yb, opts.FlagNorm);

    % Bootstrap version of R
    [Rb, FlagDim] = PLS_cov(Xb,Yb);

    % SVD of Rb
    [Ub,Sb,Vb] = svd(Rb,'econ');
    if FlagDim ==1
        tmp1 = Ub;
        tmp2 = Vb;
        clear Ub Vb
        Ub = tmp2;
        Vb = tmp1;
    end


    % Procrustes
    procmatU = rri_bootprocrust(V, Vb, U, Ub);
    Ub = Ub * Sb * procmatU;  
    Vb = Vb * Sb * procmatU;
    Sb(:,iter_boot) = sqrt(sum(Vb.^2));      

    % normalise by original singular values
    Vb = Vb./repmat(diag(S)',size(Vb,1),1);
    Ub = Ub./repmat(diag(S)',size(Ub,1),1);

    
    % store all singular values to see explained Covar in each
    % iteration
    Sboot(:,iter_boot) = Sb(:,iter_boot);
    explVarLVsBoot(:,iter_boot) = ((Sb(:,iter_boot)).^2) ./ sum((Sb(:,iter_boot).^2)); 

    % Online computing of mean and variance
    if iter_boot == 1
        u_sum = Ub;
        v_sum = Vb;
        u_sqr = Ub.^2;
        v_sqr = Vb.^2;
    else            
        u_sum = u_sum+Ub;
        v_sum = v_sum+Vb; 
        u_sqr = u_sqr+Ub.^2;
        v_sqr = v_sqr+Vb.^2;
    end

    % store bootstrap values
    Ub_boots(:,:,iter_boot) = Ub;
    Vb_boots(:,:,iter_boot) = Vb;
end

% compute standard error and bootstrap ratios
u_sum2 = (u_sum.^2)/opts.nBoot;
u_se = sqrt( (u_sqr - u_sum2) / (opts.nBoot-1) );
u_btr = U./u_se;

v_sum2 = (v_sum.^2)/opts.nBoot;
v_se = sqrt( (v_sqr - v_sum2) / (opts.nBoot-1) );
v_btr = V./v_se;

% compute bootstrap 95% CI assuming normal distribution
alpha = 0.05;
inter = [alpha/2 1-(alpha/2)];
ts = tinv(inter, opts.nBoot-1);

for k=1:min([size(Xb,2), size(Yb,2)])
    uci(:,k,:) = repmat(u_sum(:,k)/opts.nBoot, 1, 2) + ts.*u_se(:,k);
    vci(:,k,:) = repmat(v_sum(:,k)/opts.nBoot, 1, 2) + ts.*v_se(:,k);
end

% collect results
Ub_mean = u_sum/opts.nBoot;
Vb_mean = v_sum/opts.nBoot;
Ubr     = u_btr;
Vbr     = v_btr;
Uci     = uci;
Vci     = vci;
explLVBoot = explVarLVsBoot;
    
    

end
% =========================================================================




% =========================================================================
