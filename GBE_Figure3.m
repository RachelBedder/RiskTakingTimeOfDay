function [statistics] = GBE_Figure3(allData,nSample)


%% 1. Model-Free Statistics and Plots

%...extract the vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay;
yearDay         =       allData.yearDay;
runData         =       [allData.gambles allData.betas];
version         =       allData.version;

runData{6}      =       log(runData{6});    %...log transform lambda as standard in many studies

%...make all the different groups
i.PARTS     =       length(allData.timeOfDay);
i.UK        =   1;      kLabel{i.UK}    = 'UK';            kIdx{i.UK}     =  find(true(i.PARTS,1) & location==1);
i.US        =   2;      kLabel{i.US}    = 'USA';           kIdx{i.US}     =  find(location==2);
i.All       =   3;      kLabel{i.All}   = 'All';           kIdx{i.All}     =  find(location<3);

%%...identify all different data locations and what other data to compare
%%in the hypothesis tests
i.gain                =       1;                
i.loss                =       2;          
i.mix                 =       3;       
i.riskGain            =       4;        
i.riskLoss            =       5;           
i.lambda              =       6;      
i.mu                  =       7;       

compare(i.gain)       =       [i.loss]; 
compare(i.loss)       =       [i.loss];    
compare(i.mix)        =       [i.loss];
compare(i.riskGain)   =       [i.riskLoss];
compare(i.riskLoss)   =       [i.riskLoss]; 
compare(i.lambda)     =       [i.riskLoss];
compare(i.mu)         =       [i.riskLoss];
 
%...this is set for when there are any hypothesis tests comparing who data
%vectors, the data vector the current running vector needs to be compare
%to is run prior. If you are comparing everything to loss trials/alpha- you
%don't need to change this. There is an example order below to compare all
%to mixed trials and lambda.
runOrder = [i.loss i.gain i.mix i.riskLoss i.riskGain i.lambda i.mu];
%runOrder = [i.mix i.gain i.loss i.lambda i.riskGain i.riskLoss i.mu];

%...set the colours to be plotted in
colours{i.loss}       =       [[1 .4 .4];[1 0 0];[.6 0 0]];
colours{i.gain}       =       [[.4 1 .4];[0 1 0];[0 .4 0]];
colours{i.mix}        =       [[.4 .7 1];[0 0 1];[0 .3 .6]];
colours{i.riskLoss}   =       [[1 .4 .4];[1 0 0];[.6 0 0]];
colours{i.riskGain}   =       [[.4 1 .4];[0 1 0];[0 .4 0]];
colours{i.lambda}     =       [[.4 .7 1];[0 0 1];[0 .3 .6]];
colours{i.mu}         =       [[0 0 0];[0 0 0];[0 0 0]];
   
rand('seed',sum(100*clock));

%...make the 3.5 hour daily time difference bins 
timeBins          =     [0 .1449 .2908 .4367 .5828];

%%...add the time parameters
earlyTime         =     0.25+(0.0417*2);
lateTime          =     0.75+(0.0417*4);

%% ...for each participant find the first two eligible plays

withinTime  =  nan(i.PARTS,2);
withinGam   =  nan(i.PARTS,2);

for count = runOrder; %...run in order for datasets given above (i.e loss trials, gain trials lambda fits)
    
for part = 1:i.PARTS %...run for each participant
       
pTime           =       timeOfDay(part,:);
pData           =       runData{count}(part,:);
pDay            =       yearDay(part,:);
pVerse          =       version(part,:);
 
%% ...find the first eligible play
earlyPlay =  find(pTime >earlyTime & pTime <lateTime,1,'first'); 

if ~isempty(earlyPlay) %...only run for those with a first eligible play
 
timeWithin(1)  =       pTime(earlyPlay);
gamWithin(1)   =       pData(earlyPlay);   

%% ...find the next eligible play

%...transform all the data relative to first play
diffTime      =       pTime-timeWithin;
diffData      =       pData-gamWithin;
diffDay       =       pDay-pDay(earlyPlay);

%...find all other eligible plays
timeIdx       =       find(pTime>earlyTime & pTime<lateTime); %...must be later than earliest time (8am) and earlier than the latest time (10pm)
dayIdx        =       find(diffDay>=1); %...must be at least one day after the first play
verseIdx      =       find(pVerse==pVerse(earlyPlay)); %...must be the same version as the first play

allEligible   =       intersect(timeIdx,dayIdx); %...find the plays that are both a eligible time and a different day
allEligible   =       intersect(allEligible,verseIdx); %...of those above plays find those that haven't changed version

if ~isempty(allEligible)
    
    if timeWithin<pTime(min(allEligible)); %...if the first play is an earlier time, return in order [first eligible, second eligible]
    
withinTime(part,:) =     [timeWithin pTime(min(allEligible))];    
withinGam(part,:)  =     [gamWithin pData(min(allEligible))];  

    else %...if the second play is an earlier time, return in order [second eligible, first eligible]
        
withinTime(part,:) =   [pTime(min(allEligible)) timeWithin];    
withinGam(part,:)  =   [pData(min(allEligible)) gamWithin];  
        
    end

end
 
else %...don't run for those without eligible plays
    
end

end

newGam(:,count)  = withinGam(:,2)-withinGam(:,1); %...for each data set, get the gambling difference

end

newTime(:,1)         = withinTime(:,2)-withinTime(:,1); %...get the time difference 


%...add the new indexes to the subgroups to preserve the structures
i.UK        =   1;      kLabel{i.UK}    = 'UK';            kIdx{i.UK}     =  intersect(kIdx{i.UK},find(~isnan(newTime)));
i.US        =   2;      kLabel{i.US}    = 'USA';           kIdx{i.US}     =  intersect(kIdx{i.US},find(~isnan(newTime)));
i.All       =   3;      kLabel{i.All}   = 'All';           kIdx{i.All}    =  intersect(kIdx{i.All},find(~isnan(newTime)));

%% ...bootstrapping and permutation tests

for k = 1:length(kIdx); %...run for each subgroup of participants

    rSample{k}     =   nan(nSample,3);
    reBoot{k}      =   nan(nSample,3); 
    idx            =   kIdx{k};
    
 for count = runOrder; 
     
for n = 1:nSample; %...iterate for n samples for permutation tests
 
    %...randomly permutate the gamble data for each participant (preserving
    %each participants data together)
    reData                  =       newGam(idx(randperm(length(idx))),:);
    rSample{k}(n,count)     =       corr(newTime(idx),reData(:,i.riskGain),'type','spearman');
    
    diffSample{k}(n,count)=       corr(newTime(idx),reData(:,count))-corr(newTime(idx),reData(:,compare(count)),'type','spearman');
    
 clear reData
   
    %...sample without replacement for bootstrapped effect sizes
    bootParts                  =      datasample(idx,length(idx));
    reBoot{k}(n,count)         =      corr(newTime(bootParts),newGam(bootParts,count),'type','spearman'); clear bootParts

end



[Rho(k,count) pValue(k,count)]    =       corr(newTime(idx),newGam(idx,count),'type','spearman');%...get real effect size and p value            
pPerm(k,count)                    =       sum(Rho(k,count)>rSample{k}(:,count))/nSample;%...generate a permutated p value

sortBoot                          =       sort(reBoot{k}(:,count));
lowerBound(k,count)               =       Rho(k,count) - sortBoot(ceil(nSample*0.025));%...get the bootstrapped lower confidence bound
upperBound(k,count)               =       sortBoot(ceil(nSample*0.975)) - Rho(k,count);%...get the bootstrapped upper confidence bound
    
clear sortBoot  

%...generate means and standard errors for each time of day bin to plot
for tb = 1:length(timeBins)-1;
    
    todIdx                         =      intersect(idx,find(newTime>=timeBins(tb) & newTime<=timeBins(tb+1)));
    meanData{k}(tb,count)          =      mean(newGam(todIdx,count));
    seData{k}(tb,count)            =      std(newGam(todIdx,count))/sqrt(length(idx));
    
end

%...get the actual difference in Rhos
diffR(k,count)                  =       Rho(k,count)-Rho(k,compare(count));

%...get the permutated p value for the difference in Rhos

pDiffPerm(k,count)              =       sum(diffR(k,count)<=diffSample{k}(:,count))/nSample;
%% THIS MIGHT NEED CHANGING FOR LAMBDA!!!!!!! So currently use the real p value as a marker

end; clear count
end; clear k idx

%...make tables of all relevant statistics for each gender group
stats               =   {'Rho','P value','Permu P value','diff Perma P value','BS Lower Bound','BS Upper Bound', 'N Participants'}';
for k = 1:length(kIdx)
    
    riskGain        =       [Rho(k,i.riskGain) pValue(k,i.riskGain) pPerm(k,i.riskGain) pDiffPerm(k,i.riskGain) lowerBound(k,i.riskGain) upperBound(k,i.riskGain) length(kIdx{k})]';
    riskLoss        =       [Rho(k,i.riskLoss) pValue(k,i.riskLoss) pPerm(k,i.riskLoss) pDiffPerm(k,i.riskLoss) lowerBound(k,i.riskLoss) upperBound(k,i.riskLoss) length(kIdx{k})]';
    lambda          =       [Rho(k,i.lambda) pValue(k,i.lambda) pPerm(k,i.lambda) pDiffPerm(k,i.lambda) lowerBound(k,i.lambda) upperBound(k,i.lambda) length(kIdx{k})]';
    mu              =       [Rho(k,i.mu) pValue(k,i.mu) pPerm(k,i.mu) pDiffPerm(k,i.mu) lowerBound(k,i.mu) upperBound(k,i.mu) length(kIdx{k})]';
    gainPlays       =       [Rho(k,i.gain) pValue(k,i.gain) pPerm(k,i.gain) pDiffPerm(k,i.gain) lowerBound(k,i.gain) upperBound(k,i.gain) length(kIdx{k})]';
    lossPlays       =       [Rho(k,i.loss) pValue(k,i.loss) pPerm(k,i.loss) pDiffPerm(k,i.loss) lowerBound(k,i.loss) upperBound(k,i.loss) length(kIdx{k})]';
    mixPlays        =       [Rho(k,i.mix) pValue(k,i.mix) pPerm(k,i.mix) pDiffPerm(k,i.mix) lowerBound(k,i.mix) upperBound(k,i.mix) length(kIdx{k})]';
    
    
statistics.(strcat('WithinSubjects_',kLabel{k})) = table(stats,riskGain,riskLoss,lambda,mu,lossPlays,gainPlays,mixPlays);clear riskGain riskLoss lambda mu
 
end; clear k


%% ...generate figure 3a

figure;

plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4],meanData{i.All}(:,i.lambda),seData{i.All}(:,i.lambda),'Color',colours{i.riskLoss}(2,:),'LineWidth',3);
% e2 = errorbar([1:4],meanData{i.All}(:,i.loss),seData{i.All}(:,i.loss),'Color',colours{i.riskLoss}(1,:),'LineWidth',3);
set(gca,'xtick',[1:4],'xticklabel',{'3.5 hours','7 hours','10.5 hours','14 hours'});xtickangle(45)
set(gca,'ytick',[-0.1:0.02:0.04])
ylim([-0.1 0.04]);ylabel('\alpha- difference (late minus early play')
xlim([0 5]);xlabel('Time of day different (up to 14 hours)');
axis square
legend([e1],{'\alpha-'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])

figure;

plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4],meanData{i.All}(:,i.loss),seData{i.All}(:,i.loss),'Color',colours{i.loss}(1,:),'LineWidth',3);
set(gca,'xtick',[1:4],'xticklabel',{'3.5 hours','7 hours','10.5 hours','14 hours'});xtickangle(45)
set(gca,'ytick',[-0.04:0.01:0.08])
ylim([-0.04 0.08]);ylabel('% Loss chosen difference (late minus early play')
xlim([0 5]);xlabel('Time of day different (up to 14 hours)');
axis square
legend([e1],{'Loss'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])


%% ...generate figure 3c
figure; hold on
bar([1],Rho(i.All,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([2],Rho(i.All,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([3],Rho(i.All,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rho(i.All,[i.riskGain i.riskLoss i.lambda]),lowerBound(i.All,[i.riskGain i.riskLoss i.lambda]),upperBound(i.All,[i.riskGain i.riskLoss i.lambda]),'k','LineStyle', 'none')
ylabel('Within-subject effect size');ylim([-0.12 0.06])
set(gca,'xtick',[1 2 3],'xticklabel',{'\alpha+','\alpha-','log(\lambda)'},'ytick',-.12:0.02:0.06);xlim([0 4])
axis square
legend({'\alpha+','\alpha-','log(\lambda)'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])

%% ...generate figure 3c
figure; hold on
bar([1],Rho(i.All,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([2],Rho(i.All,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([3],Rho(i.All,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rho(i.All,[i.gain i.loss i.mix]),lowerBound(i.All,[i.gain i.loss i.mix]),upperBound(i.All,[i.gain i.loss i.mix]),'k','LineStyle', 'none')
ylabel('Within-subject effect size');ylim([-0.06 0.12])
set(gca,'xtick',[1 2 3],'xticklabel',{'Gain','Loss','Mix'},'ytick',-.06:0.02:0.12);xlim([0 4])
axis square
legend({'Gain','Loss','Mix'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])



save('GBE_Figure3_Data.mat')

end

