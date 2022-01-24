%% =====================================================================
%                         ANALYSIS SCRIPTS
% ======================================================================
%
% This file contains analysis code used in the following publication: 
%     Yoo et al., 2022. JoCN
%
% Experimental data and processed eyetracking and neural data are available
% on OSF: XX. 
% Some scripts cannot be run without neuroimaging data, which is available
% upon request. Please email for necessary files. 
% 
% ===== LIST OF SECTIONS =====
% - xxx
% - xxx



%% =====================================================================
%               univariate sustained activation
%======================================================================
%%  calculate trial-wise BOLD activity
% for one subject, ROI, and run, separate all trials lined up at trial start

clear all

filepath = '.';

subjid = 7;

ROIVec = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','iPCS','sPCS'};
nROIs = length(ROIVec);

switch subjid
    case {1:5}
        dim_firstitem = 8; % index of first dimension of first item location (dva) in outputmatrix
        dim_priorities = 24:27; % which columns in behavior data corresponds to item priorities
        
        switch subjid
            case 4
                subjnum = 2;
            case 5
                subjnum = 3;
            case 1
                subjnum = 4;
            case 2
                subjnum = 6;
            case 3
                subjnum = 7;
        end
        
    case {6:11}
        ppd = 39.9334;
        dim_firstitem = 9; % index of x-coord of first item location (pixels) in designMat (note that y-coord is dim_firstitem+4)
        dim_priorities = 1:4; % which columns in behavior data corresponds to item priorities
end

nRuns = 10;
switch subjid
    case 9
        nRuns = 12;
    case 7
        nRuns = 13;
    case 1
        nRuns = 14;
    case 10
        nRuns = 11;
end

nTrials = 12;
nTRs = 17;

priorityVec = [0.6 0.3 0.1 0];
nPriorities = length(priorityVec);

% ====== PRF ======

% load pRF data
data_pRF = niftiRead(sprintf('%s/%d/pRFs/RF_ss5-fFit.nii.gz',filepath, subjid));

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

for iROI = 1:nROIs % for each ROI (or section of ROIs)
% iROI = 1;
    ROI = ROIVec{iROI}
    
        % load ROI
        funcdata_roi = niftiRead(sprintf('%s/%d/ROIs/bilat.%s.nii.gz',filepath,subjid, ROI));
        idxx = idx & logical(funcdata_roi.data);
    
    % load behavioral data across all runs if current data set up
    if any(subjid==6:11)
        load(sprintf('%s/%d/experiment/%d_designMat.mat', filepath, subjid, subjid));
        
        nTRsperTrial = designMat(:,20);
        nTRsperTrial(nTRsperTrial == 8.8)  = 17;
        nTRsperTrial(nTRsperTrial == 10.1) = 18;
        nTRsperTrial(nTRsperTrial == 11.4) = 19;
    end
    
    funcdata_alltrials = [];
    for irun = 1:nRuns; % for each run
        
        % load functional data
        %         funcdata_all = niftiRead(sprintf('%s%s/functional/Pri1/surf_volreg_detrend%02d.nii.gz',filepath, subjid, irun));
        funcdata_all = niftiRead(sprintf('%s/%d/functional/Pri1/surf_volreg_normPctDet%02d.nii.gz',filepath, subjid, irun));
        funcdata_all2d = reshape(funcdata_all.data,[prod(funcdata_all.dim(1:3)) funcdata_all.dim(4)]);
        
        % obtain data for relevant voxels (nVoxels x nTimePts)
        funcdata_all = funcdata_all2d(idxx,:);
        clear funcdata_all2d
        
        % zscore all voxels
        std_idx = std(funcdata_all,[],2)==0; % for some ps, (e.g., KD), there are some voxels that are always 0, so standardizing leads to nans
        funcdata_all = bsxfun(@minus, funcdata_all, mean(funcdata_all,2));
        funcdata_all = bsxfun(@rdivide, funcdata_all, std(funcdata_all,[],2));
        funcdata_all(std_idx,:) = 0; % turn the nans bck to 0s. ASPEN: this probably reflects a bug in ROI definition
        
        % ========== PRIORITY PER QUADRANT INFORMATION =========
        
        % load run-wise behav data if Alfredo set up
        if any(subjid==1:5)
            outputMatrix = load(sprintf('%s/%d/experiment/outptMatrx_P%d_r%02d.mat', filepath, subjid, subjnum, irun));
            outputMatrix = outputMatrix.outputMatrix(end-nTrials+1:end,:);
            behavdata_priorities = outputMatrix(:,dim_priorities);
        end
        
        % how many TRs per trial
        switch subjid
            case {1,2,3,5}
                nTRsperTrial = round(outputMatrix(:,28)/1.3);
                nTRsperTrial = nTRsperTrial(end-11:end);
                cumTRs = cumsum([2; nTRsperTrial]);
            case 4
                cumTRs = [2; round(outputMatrix(end-11:end,28)/1.3)-round(outputMatrix(end,28)/1.3)+218];
                nTRsperTrial = diff(cumTRs);
            case 6:11
                currtrials = ((irun-1)*nTrials+1):(irun*nTrials);
                cumTRs = cumsum([2; nTRsperTrial(currtrials)]);
                behavdata_priorities = designMat(:,dim_priorities);
                if (subjid==7)
                    if cumTRs(end) ~= 218
                        fprintf('participnt 7 run %d is %d TRs \n',irun,cumTRs(end))
                    end
                end
        end
       
            funcdata_trial = nan(nTrials,nTRs); % 19 is the max TRs per trial
            for itrial = 1:nTrials % for each trial
                trial_start = cumTRs(itrial)+1;
                trial_end = cumTRs(itrial)+nTRs; %cumTRs(itrial+1);
                try
                funcdata_trial(itrial,1:nTRs) = mean(funcdata_all(:,trial_start:trial_end));
                catch
                    xx = mean(funcdata_all(:,trial_start:end));
                    funcdata_trial(itrial,1:length(xx)) = xx;
                end
            end
            funcdata_alltrials = [funcdata_alltrials; funcdata_trial];
    end
    
    data.(ROI) = funcdata_alltrials;
end

save(sprintf('unweighted averages/trialdata_%d.mat',subjid),'subjid','data')


%% calculating mean betas on delay period

clear all

glmname = 'shareddelay';

filepath = '/data/Pri_quad';
% filepath = '/Users/blobface/mnt/aspen@mys/DATA/data/Pri_quad';
load('plottingsettings.mat')

use_trial{7}(1:12)=[]; % delete first run for subject (was ignored in GLM)

meanbeta = nan(nSubj,nROIs);
for isubj = 1:nSubj
    subjid = subjidVec(isubj);
    
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

%% =====================================================================
%               ITEM-SPECIFIC DELAY-PERIOD ACTIVATION
%======================================================================

%% create txt files for anova
% calculate pRF-weighted betas for each participant and priority

clear all

ROI = 'sPCS';
weightedby = 'all';
weightingmethod = 'weightthresholded';
glmname = 'shareddelay';

load('plottingsettings.mat')


[bold,subj,pri] = deal(nan(nSubj,nPriorities));
for isubj = 1:nSubj
    subjid = subjidVec(isubj);
    
    if (subjid==7)
        use_trial{7}(1:12) = [];
    end
    % get prf weighted data
    load(sprintf('weighted_averages/%s_beta_pRF%s_%s_%d.mat',glmname,weightingmethod,weightedby,subjid));

    bold(isubj,:) = nanmean(beta.(ROI)(logical(use_trial{subjid}),:));
    subj(isubj,:) = isubj;
    pri(isubj,:) = priorityVec;
end

% save txt for rm ANOVA: BOLD ~ priority
BOLD = bold(:);
Priority = pri(:);
Subject = subj(:);
t = table(BOLD, Priority, Subject);
writetable(t,sprintf('txt_forANOVAs/%s_%s_%s_beta_priority_%s.txt',glmname,ROI,weightedby,weightingmethod));

%% regression: effect of priority on item-specific delay-period BOLD activity

% see final_plots.m for Figure 6B