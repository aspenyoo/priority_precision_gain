%%=====================================================================
%               univariate sustained activation
%======================================================================

%% calculating mean betas on delay period

clear all

glmname = 'shareddelay';

filepath = '/data/Pri_quad';
% filepath = '/Users/blobface/mnt/aspen@mys/DATA/data/Pri_quad';
load('plottingsettings.mat')

use_trial.CC(1:12)=[]; % delete first run for CC (was ignored in GLM)

meanbeta = nan(nSubj,nROIs);
for isubj = 1:nSubj
    subjid = subjidVec{isubj}
    
    % load pRF data
    data_pRF = niftiRead(sprintf('%s/%s/pRFs/RF_ss5-fFit.nii.gz',filepath, subjid));
    
    %  ==== initial mask of relevant voxels (based on ROI and pRF parameters) ====
    
    % threshold data by variance explained
    threshold = 0.1;
    vedim = 2; % dimension of data corresponding to variance explained
    idx = (data_pRF.data(:,:,:,vedim) >= threshold);
    
    % get coordinates of center estimates for all relevant voxels (in dva?)
    xdim = 6; % dimension of data corresponding to x coordinate
    ydim = 7; % dimension of data corresponding to y coordinate
    xdata = data_pRF.data(:,:,:,xdim);
    ydata = data_pRF.data(:,:,:,ydim);
    
    % update idx to include only prfs with min =< centers <= max dva
    dva_min = 4; % lower bound of voxels to include, based on prf center
    dva_max = 20;
    prfcenterdva = sqrt(xdata.^2 + ydata.^2);
    idx = idx & (prfcenterdva>=dva_min) & (prfcenterdva<=dva_max);

    % load glm data
    glmdata = niftiRead(sprintf('%s/%s/glm/glm_%s/Results/%s_beta.nii.gz',filepath,subjid,glmname,glmname));
    
    for iROI = 1:nROIs
        ROI = ROIVec{iROI};
        
        % load ROI
        funcdata_roi = niftiRead(sprintf('%s/%s/ROIs/bilat.%s.nii.gz',filepath,subjid,ROI));
        idxx = idx & logical(funcdata_roi.data); % which voxels
        
        meanglmbeta = squeeze(mean(glmdata.data(:,:,:,:,use_trial.(subjid)),5));
        meanbeta(isubj,iROI) = mean(meanglmbeta(idxx));
    end
end

save('shareddelay_meanbetas.mat','meanbeta','subjidVec','ROIVec')

%% bootstrapped significance test

clear all
load('plottingsettings.mat')
load('shareddelay_meanbetas.mat')

nTimes = 1e3;

pVec = nan(1,nROIs);
for iROI = 1:nROIs
    
    currROI_betas = meanbeta(:,iROI);
    
    nulldist = mean(currROI_betas(randi(nSubj,nSubj,nTimes)),1);
    pVec(iROI) = mean(nulldist <= 0);
end

%%=====================================================================
%               ITEM-SPECIFIC DELAY-PERIOD ACTIVATION
%======================================================================

%% CREATE TXT FILES FOR ANOVA

%% average betas
% save txt for rm ANOVA: BOLD ~ priority + subject

clear all

ROI = 'sPCS';
weightedby = 'all';
weightingmethod = 'weightthresholded';
glmname = 'shareddelay';

load('plottingsettings.mat')


[bold,subj,pri] = deal(nan(nSubj,nPriorities));
for isubj = 1:nSubj
    subjid = subjidVec{isubj}
    
    if strcmp(subjid,'CC')
        use_trial.CC(1:12) = [];
    end
    % get prf weighted data
    load(sprintf('weighted_averages/%s_beta_pRF%s_%s_%s.mat',glmname,weightingmethod,weightedby,subjid));

    bold(isubj,:) = nanmean(beta.(ROI)(logical(use_trial.(subjid)),:));
    subj(isubj,:) = isubj;
    pri(isubj,:) = priorityVec;
end

% save txt for rm ANOVA: BOLD ~ priority
BOLD = bold(:);
Priority = pri(:);
Subject = subj(:);
t = table(BOLD, Priority, Subject);
writetable(t,sprintf('txt_forANOVAs/%s_%s_%s_beta_priority_%s.txt',glmname,ROI,weightedby,weightingmethod));

%% REGRESSION: 
% 

clear all
clc

glmname = 'shareddelay';
load('plottingsettings.mat')

[b,p] = deal(nan(1,nROIs));
stats = cell(1,nROIs);
for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    xx = readtable(sprintf('%s_%s_all_beta_priority_weightthresholded.txt',glmname,ROI));
    xx = table2array(xx);
    
    % reshape variables
    nSubj = length(unique(xx(:,3)));
    bold = reshape(xx(:,1),nSubj,4);
    priority = reshape(xx(:,2),nSubj,4);
    
    % center bold around 0.0 error
    bold = bsxfun(@minus,bold,mean(bold,2));
    
    % rename variables
    y = bold(:);
    X = [ones(numel(priority),1) priority(:)];
    
    % regressions
    [bb,bint,r,rint,stats{iROI}] = regress(y,X); % stats: R^2, F, p, error variance
    b(iROI) = bb(2);
    p(iROI) = stats{iROI}(3);
end
% plot(bold','o-')

figure
plot(b,'k-','LineWidth',3)
hold on;
plot([0.5 nROIs+0.5],[0 0],'k-')
set(gca,'XLim',[0.5 nROIs+0.5],'XTick',1:nROIs,'XTickLabel',ROIVec)
defaultplot