function [statistics] = GBE_Figure2(allData,nSample)

%% 1. Model-Free Statistics and Plots

%...extract the vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay(:,1);
betas           =       allData.betas;
gender          =       allData.gender;
age             =       allData.age;


%...make all the different groups
i.PARTS     =           length(allData.timeOfDay);
i.UK        =   1;      kLabel{i.UK}    = 'UK';            kIdx{i.UK}     =  find(true(i.PARTS,1) & location==1);
i.Male      =   2;      kLabel{i.Male} 	= 'Male';          kIdx{i.Male}     =  find(gender==1 & location==1);
i.Female    =   3;      kLabel{i.Female}= 'Female';        kIdx{i.Female}     =  find(gender==2 & location==1);
i.Young     =   4;      kLabel{i.Young} = 'Young';         kIdx{i.Young}     =  find(age<=3 & location==1);
i.Old       =   5;      kLabel{i.Old}   = 'Old';           kIdx{i.Old}     =  find(age>=4 & location==1);
i.US        =   6;      kLabel{i.US}    = 'USA';           kIdx{i.US}     =  find(location==2);
i.All       =   7;      kLabel{i.All}   = 'All';           kIdx{i.All}     =  find(location<3);

%%...identify all different data locations and what other data to compare
%%in the hypothesis tests
      
i.riskGain            =       1;        
i.riskLoss            =       2;           
i.lambda              =       3;      
i.mu                  =       4;       

compare(i.riskGain)   =       [i.riskLoss];
compare(i.riskLoss)   =       [i.riskLoss]; 
compare(i.lambda)     =       [i.riskLoss];
compare(i.mu)         =       [i.riskLoss];
 
%...this is set for when there are any hypothesis tests comparing who data
%vectors, the data vector the current running vector needs to be compare
%to is run prior
runOrder = [i.riskLoss i.riskGain i.lambda i.mu];

colours{i.riskLoss}   =       [[1 .4 .4];[1 0 0];[.6 0 0]];
colours{i.riskGain}   =       [[.4 1 .4];[0 1 0];[0 .4 0]];
colours{i.lambda}     =       [[.4 .7 1];[0 0 1];[0 .3 .6]];
colours{i.mu}         =       [[0 0 0];[0 0 0];[0 0 0]];
   
%...make the 4 hour daily time bins
timeBins        =       {[0.2500 0.5000],[0.5000 0.7500],[0.7500 1.0000],[1 1.2500]};
   
rand('seed',sum(100*clock));

betaData        =       [betas{i.riskGain}(:,1) betas{i.riskLoss}(:,1) log(betas{i.lambda}(:,1)) betas{i.mu}(:,1)];

for count = runOrder

for k = 1:length(kIdx); %...run for each subgroup of participants

    rSample{k}     =   nan(nSample,3);
    reBoot{k}      =   nan(nSample,3); 
    idx            =   kIdx{k};
  
for n = 1:nSample; %...iterate for n samples for permutation tests
    
    %...randomly permutate the gamble data for each participant (keeping
    %their loss, gain and mixed percentage together)
    reBeta                     =       betaData(idx(randperm(length(idx))),:);
    rSample{k}(n,count)        =       corr(timeOfDay(idx),reBeta(:,count),'type','spearman'); clear reGam
    
    %...find the p value for the differences between each parameter and the
    %alpha minus
    diffSample{k}(n,count)     =       corr(timeOfDay(idx),reBeta(:,count))-corr(timeOfDay(idx),reBeta(:,compare(count)),'type','spearman');

    %...sample without replacement for bootstrapped effect sizes
    bootParts                  =      datasample(idx,length(idx));
    reBoot{k}(n,count)         =      corr(timeOfDay(bootParts),betaData(bootParts,count),'type','spearman');

    end
     
end; clear k idx

%...generate the new statistics and add them to a table
for k = 1:length(kIdx);
    
    idx            =   kIdx{k};

     [Rho(k,count) pValue(k,count)] =       corr(timeOfDay(idx),betaData(idx,count),'type','spearman');%...get real effect size and p value            
     pPerm(k,count)                  =       sum(Rho(k,count)>rSample{k}(:,count))/nSample;%...generate a permutated p value

    sortBoot                         =       sort(reBoot{k}(:,count));
    lowerBound(k,count)              =       Rho(k,count) - sortBoot(ceil(nSample*0.025));%...get the bootstrapped lower confidence bound
    upperBound(k,count)              =       sortBoot(ceil(nSample*0.975)) - Rho(k,count);%...get the bootstrapped upper confidence bound
    
clear sortBoot  

%...generate means and standard errors for each time of day bin to plot
for tod = 1:4
    
    todIdx                          =      intersect(idx,find(timeOfDay>=timeBins{tod}(1) & timeOfDay<=timeBins{tod}(2)));
    meanBeta{k}(tod,count)          =      mean(betaData(todIdx,count));
    seBeta{k}(tod,count)            =      std(betaData(todIdx,count))/sqrt(length(idx));
    
end

%...get the actual difference in Rhos
diffR(k,count)                  =       Rho(k,count)-Rho(k,compare(count));

%...get the permutated p value for the difference in Rhos
pDiffPerm(k,count)              =       sum(diffR(k,count)<=diffSample{k}(:,count))/nSample;

end; clear kidx
end; clear count

%...make tables of all relevant statistics for each gender group
stats               =   {'Rho','P value','Permu P value','diff Perma P value','BS Lower Bound','BS Upper Bound'}';
for k = 1:length(kIdx)
    
    riskGain        =       [Rho(k,i.riskGain) pValue(k,i.riskGain) pPerm(k,i.riskGain) pDiffPerm(k,i.riskGain) lowerBound(k,i.riskGain) upperBound(k,i.riskGain)]';
    riskLoss        =       [Rho(k,i.riskLoss) pValue(k,i.riskLoss) pPerm(k,i.riskLoss) pDiffPerm(k,i.riskLoss) lowerBound(k,i.riskLoss) upperBound(k,i.riskLoss)]';
    lambda          =       [Rho(k,i.lambda) pValue(k,i.lambda) pPerm(k,i.lambda) pDiffPerm(k,i.lambda) lowerBound(k,i.lambda) upperBound(k,i.lambda)]';
    mu              =       [Rho(k,i.mu) pValue(k,i.mu) pPerm(k,i.mu) pDiffPerm(k,i.mu) lowerBound(k,i.mu) upperBound(k,i.mu)]';
    
statistics.(strcat('BetaTimeofDay_',kLabel{k})) = table(stats,riskGain,riskLoss,lambda,mu);clear riskGain riskLoss lambda mu
 
end; clear k


%% ...generate figure 1b

100*meanBeta{i.UK}(:,i.riskLoss)/meanBeta{i.UK}(:,i.riskLoss)

figure;
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4],100*(meanBeta{i.UK}(:,i.riskLoss)./meanBeta{i.UK}(1,i.riskLoss)-1),seBeta{i.UK}(:,i.riskLoss)*100,'Color',colours{i.riskLoss}(1,:),'LineWidth',3);
e2 = errorbar([1:4],100*(meanBeta{i.US}(:,i.riskLoss)./meanBeta{i.US}(1,i.riskLoss)-1),seBeta{i.US}(:,i.riskLoss)*100,'Color',colours{i.riskLoss}(3,:),'LineWidth',3);
% e3 = errorbar([1:4],meanBeta{i.All}(:,i.lambda)-meanBeta{i.All}(1,i.lambda),seBeta{i.All}(:,i.lambda),'Color',colours{i.lambda}(2,:),'LineWidth',3);
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'})
% set(gca,'ytick',[-0.01:0.01:0.05],'yticklabel',{'-1%','0%','1%','2%','3%','4%','5%'})
ylim([-8 2]);ylabel('\alpha- parameter estimate (percent relative to 6am)')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2],{'UK','US'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])

figure
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4],100*(meanBeta{i.UK}(:,i.lambda)./meanBeta{i.UK}(1,i.lambda)-1),seBeta{i.UK}(:,i.lambda)*100,'Color',colours{i.lambda}(1,:),'LineWidth',3);
e2 = errorbar([1:4],100*(meanBeta{i.US}(:,i.lambda)./meanBeta{i.US}(1,i.lambda)-1),seBeta{i.US}(:,i.lambda)*100,'Color',colours{i.lambda}(3,:),'LineWidth',3);
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'})
% set(gca,'ytick',[-0.01:0.01:0.05],'yticklabel',{'-1%','0%','1%','2%','3%','4%','5%'})
ylim([-2 20]);ylabel('log(\lambda) parameter estimate (percent relative to 6am)')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2],{'UK','US'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])


%% ...generate figure 1c
figure; hold on
bar([1],Rho(i.UK,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([2],Rho(i.UK,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([3],Rho(i.UK,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rho(i.UK,[1:3]),lowerBound(i.UK,[1:3]),upperBound(i.UK,[1:3]),'k','LineStyle', 'none')
ylabel('Effect sizes for time of day');ylim([-0.06 0.04])
set(gca,'xtick',[1 2 3],'xticklabel',{'\alpha+','\alpha-','log(\lambda)'},'ytick',-0.06:0.02:0.04);xtickangle(45);xlim([0 4])
axis square
legend({'\alpha+','\alpha-','log(\lambda)'})
title([kLabel{i.UK},' N = ',num2str(length(kIdx{i.UK}))])

%% ...generate figure 1c
figure; hold on
bar([1],Rho(i.UK,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([2],Rho(i.UK,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([3],Rho(i.UK,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
bar([5],Rho(i.US,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([6],Rho(i.US,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([7],Rho(i.US,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rho(i.UK,[1:3]),lowerBound(i.UK,[1:3]),upperBound(i.UK,[1:3]),'k','LineStyle', 'none')
errorbar([5:7],Rho(i.US,[1:3]),lowerBound(i.US,[1:3]),upperBound(i.US,[1:3]),'k','LineStyle', 'none')
ylabel('Effect sizes for time of day');ylim([-0.06 0.06])
set(gca,'xtick',[2 6],'xticklabel',{[kLabel{i.UK},' N = ',num2str(length(kIdx{i.UK}))],[kLabel{i.US},' N = ',num2str(length(kIdx{i.US}))]},...
    'ytick',-0.06:0.02:0.06);xtickangle(45);xlim([0 8])
axis square
legend({'\alpha+','\alpha-','log(\lambda)'})




save('GBE_Figure2_Data.mat')

end

% figure;
% 
% plot([0 5],[0 0],'k'); hold on
% e1 = errorbar([1:4],meanBeta{i.All}(:,i.riskGain)-meanBeta{i.All}(1,i.riskGain),seBeta{i.All}(:,i.riskGain),'Color',colours{i.riskGain}(2,:),'LineWidth',3);
% e2 = errorbar([1:4],meanBeta{i.All}(:,i.riskLoss)-meanBeta{i.All}(1,i.riskLoss),seBeta{i.All}(:,i.riskLoss),'Color',colours{i.riskLoss}(2,:),'LineWidth',3);
% e3 = errorbar([1:4],meanBeta{i.All}(:,i.lambda)-meanBeta{i.All}(1,i.lambda),seBeta{i.All}(:,i.lambda),'Color',colours{i.lambda}(2,:),'LineWidth',3);
% set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'})
% % set(gca,'ytick',[-0.01:0.01:0.05],'yticklabel',{'-1%','0%','1%','2%','3%','4%','5%'})
% ylim([-0.1 0.1]);ylabel('% Gambles chosen (percent relative to 6am)')
% xlim([0 5]);xlabel('Time of day');
% axis square
% legend([e1 e2 e3],{'\alpha+','\alpha-','log(\lambda)'})
% title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])



% %% ... generate the error bar time of day plot for various group comparisons
% plotData{1} = [2 3]; %Male and Female
% plotData{2}  = [4 5]; %Ratio and Uncorrelated
% plotData{3}  = [6 7]; %Old and Young
% plotData{4}  = [1 8]; %UK and US
% 
% plotData = [2 3 4 5 6 7];
% % plotData = [1 8];
% 
% figure;
% 
% plotMean        =       nan(4,10);      plotSE      =       nan(4,10);
% plotColors      =       repmat([0 0 0;colours{i.riskLoss}([1 3],:)],length(plotData)/2,1);
% plotMean(:,[2 3 5 6 8 9])  =   meanBeta{i.riskLoss}(:,plotData)-meanBeta{i.riskLoss}(1,plotData);
% plotSE(:,[2 3 5 6 8 9])    =   seBeta{i.riskLoss}(:,plotData);
% xPositions      =       [[1:4];[1:4];[1:4];[1:4];[6:9];[6:9];[6:9];[11:14];[11:14];[11:14]]';
% 
% c = axes;c.ColorOrder = plotColors; c.NextPlot = 'add';
% c = errorbar(xPositions,plotMean,plotSE,'LineWidth',5);
% plot([0 15],[0 0],'k')
% set(gca,'xtick',[1:4],'xticklabel',{'6am','12pm','6pm','12pm'})
% xlim([0 5]);xlabel('Time of Day (Start of bin)');
% pbaspect([2,1,1])
% lgd = legend(kLabel{plotData{k}},'location','southwest');lgd.FontSize = 20;
