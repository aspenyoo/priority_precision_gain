%% =====================================================================
%                         PLOTTING SCRIPTS
% ======================================================================
%
% This file contains all necessary code to generate the plots used in the
% following publication: 
%     Yoo et al., 2022. JoCN
% 
% If plots are in the main manuscript, the corresponding figure number is
% labelled. Plotting scripts without a figure label are supplementary, and
% do not exist in the main manuscript. 
%
% ===== LIST OF SECTIONS =====
% - Figure 1B & 1C: behavioral main effects across participants
% - univariate BOLD time courses, single participant across ROIs
% - Figure 4A & 5A: univariate BOLD time courses, across participants and ROIs
% - Figure 4B & 5B: item-specific BOLD time courses, single participant across ROIs
% - item-specific BOLD time courses, across participants and ROIs
% - Figure 4C & 5C: indvl and group pRF-weighted betas for each ROI
% - Figure 6A: summary plot - univariate delay-period BOLD activity
% - Figure 6B: summary plot - item-specific delay-period BOLD activity
% - Figure 6C: summary plot - correlation between neural and behavioral
%   effects of priority for ROI V1. inset: summary plot - correlation 
%   between neural and behavioral effects of priority across ROIs
% - plot different pRF sigma weighting strategies, across participants

%% Figure 1B & 1C: behavioral main effects across participants
clear all

load('maineffects_behavioral.mat')
load('plottingsettings.mat')

% rename for easier plotting
xx{1} = fErr;
xx{2} = pRT*1000;
titleVec = {'error', 'RT'};

figure;
for ii = 1:2
    subplot(1,2,ii); hold on;
    plot(fliplr(xx{ii})','Color',0.8*ones(1,3))
    
    m = mean(xx{ii});
    plot(fliplr(m),'k-')
    for ipriority = 1:3
        plot(4-ipriority,m(ipriority),'.','MarkerSize',24,'Color',colorMat(ipriority,:));
    end
    if ii==1
        ylim([1 5])
        ylabel('error (dva)')
    else
        ylim([300 700])
        ylabel('RT (ms)')
    end
    xlabel('priority')
    xlim([0.5 3.5])
    set(gca,'XTick',1:3,'XTickLabel',[0.1 0.3 0.6])
    title(titleVec{ii})
    defaultplot
end

% % ===== look at diff in subject based on exlusion criteria =====
% badsubjid = [2 6 9 11];
% goodsubjid = 1:11;
% 
% figure;
% 
% subplot(1,2,1);
% plot(fliplr(xx{2}(goodsubjid,:))','k-')
% hold on;
% plot(fliplr(xx{2}(badsubjid,:))','r-')
% defaultplot
% ylim([1 5])
% ylabel('error (dva)')
% xlim([0.5 3.5])
% set(gca,'XTick',1:3,'XTickLabel',[0.1 0.3 0.6])
% title(titleVec{2})
% 
% subplot(1,2,2);
% plot(fliplr(xx{3}(goodsubjid,:))','k-')
% hold on;
% plot(fliplr(xx{3}(badsubjid,:))','r-')
% defaultplot
% ylim([300 700])
% ylabel('RT (ms)')
% xlim([0.5 3.5])
% set(gca,'XTick',1:3,'XTickLabel',[0.1 0.3 0.6])
% title(titleVec{3})

%% univariate BOLD time courses, single participant across ROIs

clear all
filepath = 'unweighted_averages';
load('plottingsettings.mat')
nTRs = 17;
TR = 1.3;

delaytime = ([1.6 11.7]./TR);
figure('Position',[200 400 1200 400]);

% load data
isubj = 2;
load(sprintf('%s/trialdata_%d.mat',filepath,isubj))

for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    activityMat = data.(ROI)(logical(use_trial{isubj}),1:nTRs);
    
    activityMat = bsxfun(@minus,activityMat(:,1:nTRs),activityMat(:,1));
    m = mean(activityMat);
    sem = std(activityMat)./sqrt(nSubj);
    
    subplot(2,5,iROI);
    fill([1:nTRs, nTRs:-1:1]-1, [m+sem fliplr(m-sem)],0.8*ones(1,3),'EdgeColor','none')
    hold on;
    plot(0:(nTRs-1),m,'Color',0.3*ones(1,3));
    plot([0 nTRs*TR],[0 0],'k-')
    xlim([0 nTRs-1])
    secVec = [0 5 10 15 20];
    set(gca,'XTick',secVec./TR,'XTickLabel',secVec);%,'YTick',[])
    ylim([-0.5 1.3])
    plot([delaytime; delaytime],[[-0.5 1]' [-0.5 1]'],'Color',0.7*ones(1,3));
    title(ROI)
    xlabel('time (sec)')
    defaultplot
end


%% Figure 4A & 5A: univariate BOLD time courses, across participants and ROIs

clear all
filepath = 'unweighted_averages';
load('plottingsettings.mat')
nTRs = 17;
TR = 1.3;

delaytime = ([1.6 11.7]./TR);
figure('Position',[200 400 1200 400]);
for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    
    activityMat = nan(nSubj,17);
    for isubj = 1:nSubj
        subjid = subjidVec(isubj);
        
        load(sprintf('%s/trialdata_%d.mat',filepath,subjid))
        activityMat(isubj,:) = nanmean(data.(ROI)(logical(use_trial{subjid}),1:nTRs));
    end
    activityMat = bsxfun(@minus,activityMat(:,1:nTRs),activityMat(:,1));
    m = mean(activityMat);
    sem = std(activityMat)./sqrt(nSubj);
    
    subplot(2,5,iROI);
    fill([1:nTRs, nTRs:-1:1]-1, [m+sem fliplr(m-sem)],0.8*ones(1,3),'EdgeColor','none')
    hold on;
    plot(0:(nTRs-1),m,'Color',0.3*ones(1,3));
    plot([0 nTRs*TR],[0 0],'k-')
    xlim([0 nTRs-1])
    secVec = [0 5 10 15 20];
    set(gca,'XTick',secVec./TR,'XTickLabel',secVec);%,'YTick',[])
    ylim([-0.5 1])
    plot([delaytime; delaytime],[[-0.5 1]' [-0.5 1]'],'Color',0.7*ones(1,3));
    title(ROI)
    xlabel('time (sec)')
    defaultplot
    
end

%% Figure 4B & 5B: item-specific BOLD time courses, single participant across ROIs

clear all

subjid = 2;
weightedBy = 'all';

datafilepath = 'weighted_averages';
load(sprintf('%s/data_pRFweighted_%s_%d.mat',datafilepath,weightedBy,subjid))
load('plottingsettings.mat')

% trial exclusion
use_trial = use_trial{subjid};
TRend = 17;
TR = 1.3; % seconds

delaytime = ([1.6 11.7]./TR);
figure('Position',[200 400 1200 400]);
for iROI = 1:nROIs
    ROI = ROIVec{iROI};

    currdat = data.(ROI)(:,:,use_trial);
    currdat = bsxfun(@minus,currdat,mean(currdat(:,1,:),1)); % making it start on average 0
    currdat = currdat./sum(nanmean(nanmean(currdat,3)));      % normalizing
    
    subplot(2,5,iROI);
    hold on
    for ipriority = 1:nPriorities
        M = nanmean(currdat(ipriority,:,:),3);
        SEM = nanstd(currdat(ipriority,:,:),[],3)/sqrt(sum(use_trial));
        M = M(1:TRend);
        SEM = SEM(1:TRend);
        
        fill([1:TRend, TRend:-1:1]-1,[M-SEM, fliplr(M+SEM)],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha','0.3')
        plot(0:(TRend-1),M,'Color',colorMat(ipriority,:));
    end
    plot([0 TRend],[0 0],'k-')
    title(ROI)
    defaultplot
    xlim([0 TRend-1])

    ylims = get(gca,'YLim');
    plot([delaytime; delaytime],[ylims' ylims'],'Color',0.7*ones(1,3));
    ylim(ylims);
    secVec = [0 5 10 15 20];
    set(gca,'XTick',secVec./TR,'XTickLabel',secVec);
    xlabel('time (sec)')
    
end

%% item-specific BOLD time courses, across participants and ROIs

clear all

load('plottingsettings.mat') % load stuff consistent across plots
nSubj = length(subjidVec);

weightingmethod = 'weighted';
weightedBy = 'all';
TRend = 17;
TR = 1.3; % seconds
secVec = [0 5 10 15 20];

TRs_period = 1:TRend;


figure('Position',[200 400 1200 400]);
means = nan(nSubj,nROIs,nPriorities,length(TRs_period));

for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    
    for isubj = 1:nSubj
        subjid = subjidVec(isubj);
        
        % load weighted average
        load(sprintf('weighted_averages/data_pRF%s_%s_%d.mat',weightingmethod,weightedBy,subjid))

        data_relTRs = data.(ROI)(:,TRs_period,:); % relevant TRs of data
        use_trial{subjid} = logical(use_trial{subjid});
        means(isubj,iROI,:,:) = nanmean(data_relTRs(:,:,use_trial{subjid}),3);
    end
    
    % demean
    means = bsxfun(@minus,means,mean(means(:,:,:,1),3));
    
    % plot
    subplot(2,5,iROI); hold on
    for ipriority = 1:nPriorities
        M = squeeze(nanmean(means(:,iROI,ipriority,:),1))';
        SEM = squeeze(nanstd(means(:,iROI,ipriority,:)))'./sqrt(nSubj);
        
        fill([1:TRend, TRend:-1:1],[M-SEM, fliplr(M+SEM)],colorMat(ipriority,:),'EdgeColor','none','FaceAlpha','0.4')
        plot(TRs_period,M,'Color',colorMat(ipriority,:));
    end
    defaultplot
    ylims = get(gca,'YLim');
    plot([2 9; 2 9],[ylims' ylims'],'--','Color',0.7*ones(1,3));
    xlim([1 TRend])
    set(gca,'XTick',secVec./TR,'XTickLabel',secVec);%,'YTick',[])
    xlabel('time (sec)')
end

%% Figure 4C & 5C: indvl and group pRF-weighted betas for each ROI

clear all

load('plottingsettings.mat')
filedir='txt_forANOVAs';

weightingmethod = 'weightthresholded';
glmname = 'shareddelay'; % [];%
weightedBy = 'all';
specialsubjid = 2;


figure('Position',[200 400 1200 400]);
for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    r = importdata(sprintf('%s/%s_%s_%s_beta_priority_%s.txt',filedir,glmname,ROI,weightedBy,weightingmethod));
    beta = reshape(r.data(:,1), [nSubj 4]);
    
    subplot(2,5,iROI); hold on
    plot(fliplr(beta)','Color',0.7*ones(1,3))
    plot(fliplr(beta(specialsubjid,:)),'Color',0.5*ones(1,3),'LineWidth',2)
    plot(fliplr(mean(beta)),'k-','LineWidth',2)
    
    for ipriority=1:nPriorities
        plot(5-ipriority,mean(beta(:,ipriority)),'.','Color',colorMat(ipriority,:),'MarkerSize',24)
    end
    
    xlim([0.5 4.5])
    set(gca,'XTick',1:4,'XTickLabel',fliplr(priorityVec))
    if any(iROI == [1 2 3 4])
        ylim([-0.15 0.15]) 
    else
        ylim([0 0.6])
    end
    title(ROI)
    defaultplot
end


%% Figure 6A: summary plot - univariate delay-period BOLD activity
% indvl subject grey, average black
clear all

load('plottingsettings.mat')
load('shareddelay_meanbetas.mat')

figure;
plot(meanbeta','Color',0.8*ones(1,3)); hold on;
plot(mean(meanbeta),'k-','LineWidth',3)
plot([1 10],[0 0],'k-')
set(gca,'XTick',1:10,'XTickLabel',ROIVec)
defaultplot
xlabel('ROI')
ylabel('univariate delay period activity')

% % looking at diff between good and bad subjects (bad = high exclusion
% % criteria)
% figure;
% badsubjid = [2 6 9 11];
% goodsubjid = 1:11;
% plot(meanbeta(goodsubjid,:)','k-')
% hold on
% plot(meanbeta(badsubjid,:)','r-')
% defaultplot
% set(gca,'XTick',1:10,'XTickLabel',ROIVec)
% ylabel('univariate delay period activity')
% xlabel('ROI')

%% Figure 6B: summary plot - item-specific delay-period BOLD activity

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

%% Figure 6C: 
% - summary plot - correlation between neural and behavioral
%   effects of priority for ROI V1
% - inset: summary plot - correlation between neural and
%   behavioral effects of priority across ROIs

clear all

load('plottingsettings.mat')
filedir='txt_forANOVAs';

weightingmethod = 'weightthresholded';
glmname = 'shareddelay'; % [];%
weightedBy = 'all';

x = priorityVec(1:3);

% ========== get behavioral data ========
load('maineffects_behavioral.mat','fErr','pRT')

% rename for easier plotting
behav_error = fErr;
behav_RT = pRT*1000;

% get behavioral slopes
for isubj = 1:nSubj
    y = behav_error(isubj,:);
    p = polyfit(x,y,1);
    behav_slopes(isubj,1) = p(1);
    
    y = behav_RT(isubj,:);
    p = polyfit(x,y,1);
    behav_slopes(isubj,2) = p(1);
end

% ======= get neural slopes =======
for iROI = 1:nROIs
    ROI = ROIVec{iROI}
%     sprintf('%s/%s%s_%s_beta_priority_%s.txt',filedir,file_prefix,ROI,weightedBy,weightingmethod)
    r = importdata(sprintf('%s/%s_%s_%s_beta_priority_%s.txt',filedir,glmname,ROI,weightedBy,weightingmethod));
    neural_data = reshape(r.data(:,1), [nSubj 4]);
    
    for isubj = 1:nSubj
        y = neural_data(isubj,1:3);
        p = polyfit(x,y,1);
        neural_slopes(isubj,iROI) = p(1);
    end
end

[rhoMat,pMat,CIMat] = deal(nan(2,nROIs));
for iROI = 1:nROIs
    [rho,p,lb,ub] = corrcoef([behav_slopes(:,1) neural_slopes(:,iROI)]);
    rhoMat(1,iROI) = rho(2);
    pMat(1,iROI) = p(2);
    CIMat(:,iROI) = [lb(2); ub(2)];
    
end

figure;
iROI = 1;
ROI = ROIVec{iROI}; % 'V1';
plot(behav_slopes(:,1), neural_slopes(:,iROI),'k.','MarkerSize',24)
title(ROI);
ylabel('neural slope')
xlabel('behavioral slope')
defaultplot

% % for all ROIs
% figure('Position',[200 400 1200 400]);
% for iROI = 1:nROIs
%     ROI = ROIVec{iROI}
%     
%     subplot(2,5,iROI)
% %     behav_norm = behav_slopes(:,1) - mean(behav_slopes(:,1));
% 
%     plot(behav_slopes(:,1), neural_slopes(:,iROI),'ko')
%     title(ROI);
%     if mod(iROI,5) == 1
%         ylabel('neural slope')
%     end
%     if ceil(iROI/5)==2
%         xlabel('behavioral slope')
%     end
%     defaultplot
% end

figure;
errorbar(1:10,rhoMat(1,:),CIMat(1,:)-rhoMat(1,:),CIMat(2,:)-rhoMat(1,:),'k-')
set(gca,'Xtick',1:10,'XTickLabel',ROIVec); hold on
plot([0.5 10.5],[0 0],'k-')
xlim([0.5 10.5])
xlabel('ROI')
ylabel('correlation between behavioral and neural priority effect')
defaultplot

%% plot different pRF sigma weighting strategies, across participants

clear all

% variables to change
weightingmethod = 'weightthresholded';  % how to weight pRF info: 'weighted','weightthresholded','none'
weightedBy = 'all';                     % how much of visual field to weight: 'all', 'hemifield', 'quad'
glmname = 'shareddelay';                % which glm to use: 'shareddelay','sharedtarg','sharedcue'

sigmaweightVec = {'true','average',1,2,4,6,8};
nSigmas = length(sigmaweightVec);

weightthresh = 1e-3;

load('plottingsettings.mat')
filepath = [];

betaMat = nan(nSubj,nROIs,nSigmas); % beta's weighted by sigma
averagesigmaMat = nan(nSubj,nROIs);
for isubj = 1:nSubj
%     subjid = subjidVec{isubj};
    
    for isigma = 1:nSigmas
        sigmaweight = sigmaweightVec{isigma};
        
        % load weightedaverages
        if ischar(sigmaweight)
            switch sigmaweight
                case 'average'
                    load(sprintf('weighted_averages/%s_beta_averagepRFsigma_%s_%s_%d.mat',glmname,weightingmethod,weightedBy,isubj),'beta','BVec')
                    averagesigmaMat(isubj,:) = BVec;
                case 'true'
                    load(sprintf('weighted_averages/%s_beta_pRF%s_%s_%d.mat',glmname,weightingmethod,weightedBy,isubj),'beta')
            end
        else
            load(sprintf('weighted_averages/%s_beta_sigma%0.2f_%s_%s_%d.mat',glmname,sigmaweight,weightingmethod,weightedBy,isubj),'beta')
        end
        
        if ~exist('priMat')
            priMat = repmat([0.6 0.3 0.1 0],size(beta.V1,1),1);
        end
        
        for iROI = 1:nROIs
            ROI = ROIVec{iROI};
            
            b = regress(beta.(ROI)(:),[ones(numel(beta.(ROI)),1) priMat(:)]);
            betaMat(isubj,iROI,isigma) = b(2);
        end
    end
    clear priMat
end

% get averages
m_beta = squeeze(mean(betaMat)); % nROIs x nSigma
sem_beta = squeeze(std(betaMat)./sqrt(nSubj));

m_averagesigma = mean(averagesigmaMat);
sem_averagesigma = std(averagesigmaMat)./sqrt(nSubj);

figure;
for iROI = 1:nROIs
    ROI = ROIVec{iROI};
    
    subplot(2,5,iROI); hold on;
    fill([0 10 10 0], m_beta(iROI,1)+[-1 -1 1 1].*sem_beta(iROI,1),...
        [0.8 0.7 0.7],'EdgeColor','none') % true
    errorbar(m_averagesigma(iROI),m_beta(iROI,2),sem_beta(iROI,2),'b*') % average
    errorbar([1,2,4,6,8],m_beta(iROI,3:end),sem_beta(iROI,3:end),'k-')
    plot([0 10],[0 0],'k-')
    
    if (iROI==5);legend({'true','average','fixed value'}); end
    if any(iROI == [1 6]); ylabel('slope'); end
    if (iROI >5); xlabel('pRF sigma value'); end
    title(ROI)
    defaultplot
end
