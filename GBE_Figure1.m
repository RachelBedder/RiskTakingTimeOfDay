function [statistics] = GBE_Figure1(allData,nSample)

%% 1. Model-Free Statistics and Plots

%...extract the vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay(:,1);
gambles         =       allData.gambles;
gender          =       allData.gender;
age             =       allData.age;

%...make all the different groups
i.PARTS     =           length(allData.timeOfDay);
i.UK        =   1;      kLabel{i.UK}    = 'UK';            kIdx{i.UK}     =  find(true(i.PARTS,1) & location==1);
i.Male      =   2;      kLabel{i.Male} 	= 'Male';          kIdx{i.Male}   =  find(gender==1 & location==1);
i.Female    =   3;      kLabel{i.Female}= 'Female';        kIdx{i.Female} =  find(gender==2 & location==1);
i.Young     =   4;      kLabel{i.Young} = 'Young';         kIdx{i.Young}  =  find(age<=3 & location==1);
i.Old       =   5;      kLabel{i.Old}   = 'Old';           kIdx{i.Old}    =  find(age>=4 & location==1);
i.US        =   6;      kLabel{i.US}    = 'USA';           kIdx{i.US}     =  find(location==2);
i.All       =   7;      kLabel{i.All}   = 'All';           kIdx{i.All}    =  find(location<3);


%%...identify all different data locations and what other data to compare
%%in the hypothesis tests

i.gain                =       1;                
i.loss                =       2;          
i.mix                 =       3;       

compare(i.gain)       =       [i.loss]; 
compare(i.loss)       =       [i.loss];    
compare(i.mix)        =       [i.loss];

%...this is set for when there are any hypothesis tests comparing who data
%vectors, the data vector the current running vector needs to be compare
%to is run prior
runOrder = [i.loss i.gain i.mix];

colours{i.loss}       =       [[1 .4 .4];[1 0 0];[.6 0 0]];
colours{i.gain}       =       [[.4 1 .4];[0 1 0];[0 .4 0]];
colours{i.mix}        =       [[.4 .7 1];[0 0 1];[0 .3 .6]];

%...make the 4 hour daily time bins
timeBins        =       {[0.2500 0.5000],[0.5000 0.7500],[0.7500 1.0000],[1 1.2500]};
   
rand('seed',sum(100*clock));

gamData         =       [gambles{i.gain}(:,1) gambles{i.loss}(:,1) gambles{i.mix}(:,1)];

for count = runOrder;

for k = 1:length(kIdx); %...run for each subgroup of participants

    rSample{k}     =   nan(nSample,3);
    reBoot{k}      =   nan(nSample,3); 
    idx            =   kIdx{k};
  
for n = 1:nSample; %...iterate for n samples for permutation tests
    
    %...randomly permutate the gamble data for each participant (keeping
    %their loss, gain and mixed percentage together)
    reGam                      =       gamData(idx(randperm(length(idx))),:);       %to check this, sort(idx(randperm(length(idx)))) should be identical to idx
    rSample{k}(n,count)        =       corr(timeOfDay(idx),reGam(:,count),'type','spearman'); clear reGam
    
    %...sample without replacement for bootstrapped effect sizes
    bootParts                  =      datasample(idx,length(idx));
    reBoot{k}(n,count)         =      corr(timeOfDay(bootParts),gamData(bootParts,count),'type','spearman');
    
    %...sample without replacement for bootstrapped means for each time of day bin
    if k == 3; 
    for tod = 1:length(timeBins)
        
        todIdx                         =      find(timeOfDay>=timeBins{tod}(1) & timeOfDay<=timeBins{tod}(2));
        
        bootParts                      =      datasample(todIdx,length(todIdx));
        meanBoot{tod}(n,count)         =      mean(gamData(bootParts,count));clear bootParts todIdx
        
    end
    end
    
    end
     
end; clear k idx


for k = 1:length(kIdx);
    
    idx            =   kIdx{k};
    
    [Rho(k,count) pValue(k,count)]  =       corr(timeOfDay(idx),gamData(idx,count),'type','spearman');            %...get real effect size and p value
    pPerm(k,count)                  =       (sum(Rho(k,count)<rSample{k}(:,count))/nSample)  ;    %...get the permutated p value
                                           
    sortBoot                        =       sort(reBoot{k}(:,count));
    lowerBound(k,count)             =       Rho(k,count) - sortBoot(ceil(nSample*0.025));     %...get the bootstrapped lower confidence bound
    upperBound(k,count)             =       sortBoot(ceil(nSample*0.975)) - Rho(k,count);     %...get the bootstrapped upper confidence bound

    
clear sortBoot  

%...generate means and standard errors for each time of day bin to plot
for tod = 1:length(meanBoot)
    
    todIdx                          =      find(timeOfDay>=timeBins{tod}(1) & timeOfDay<=timeBins{tod}(2));
    idx                             =      intersect(todIdx,kIdx{k});
    meanGam{k}(tod,count)           =      mean(gamData(idx,count));
    seGam{k}(tod,count)             =      std(gamData(idx,count))/sqrt(length(idx));

end

end; clear k idx



end; clear count

%...make tables of all relevant statistics for each gender group
stats               =   {'Rho','P value','Permu P value','BS Lower Bound','BS Upper Bound','N Participants'}';
for k = 1:length(kIdx);
    
    gainPlays       =       [Rho(k,i.gain) pValue(k,i.gain) pPerm(k,i.gain) lowerBound(k,i.gain) upperBound(k,i.gain) length(kIdx{k})]';
    lossPlays       =       [Rho(k,i.loss) pValue(k,i.loss) pPerm(k,i.loss) lowerBound(k,i.loss) upperBound(k,i.loss) length(kIdx{k})]';
    mixPlays        =       [Rho(k,i.mix) pValue(k,i.mix) pPerm(k,i.mix) lowerBound(k,i.mix) upperBound(k,i.mix) length(kIdx{k})]';
    
    
statistics.(strcat('GambingTimeofDay_',kLabel{k})) = table(stats,gainPlays,lossPlays,mixPlays);clear LossPlays GainPlays MixPlays
 
end; clear k



%% ...generate figure 1b

figure;
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4],100*(meanGam{i.UK}(:,i.gain)./meanGam{i.UK}(1,i.gain)-1),seGam{i.All}(:,i.gain)*100,'Color',colours{i.gain}(2,:),'LineWidth',3);
e2 = errorbar([1:4],100*(meanGam{i.UK}(:,i.loss)./meanGam{i.UK}(1,i.loss)-1),seGam{i.All}(:,i.loss)*100,'Color',colours{i.loss}(2,:),'LineWidth',3);
e3 = errorbar([1:4],100*(meanGam{i.UK}(:,i.mix)./meanGam{i.UK}(1,i.mix)-1),seGam{i.All}(:,i.mix)*100,'Color',colours{i.mix}(2,:),'LineWidth',3);
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'})
% set(gca,'ytick',[-0.01:0.01:0.05],'yticklabel',{'-1%','0%','1%','2%','3%','4%','5%'})
set(gca,'ytick',[-2:2:10])
ylim([-2 10]);ylabel('% Gambles chosen (percent relative to 6am)')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2 e3],{'Gain','Loss','Mix'})
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])


%% ...generate figure 1c
figure; hold on
bar([1],Rho(i.UK,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([2],Rho(i.UK,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([3],Rho(i.UK,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
bar([5],Rho(i.US,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([6],Rho(i.US,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([7],Rho(i.US,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rho(i.UK,:),lowerBound(i.UK,:),upperBound(i.UK,:),'k','LineStyle', 'none')
errorbar([5:7],Rho(i.US,:),lowerBound(i.US,:),upperBound(i.US,:),'k','LineStyle', 'none')
ylabel('Effect sizes for time of day');ylim([-0.02 0.06])
set(gca,'xtick',[2 6],'xticklabel',{[kLabel{i.UK},' N = ',num2str(length(kIdx{i.UK}))],[kLabel{i.US},' N = ',num2str(length(kIdx{i.US}))]...
    },'ytick',-0.02:0.02:0.06);xtickangle(45);xlim([0 8])
axis square
legend({'Gain','Loss','Mix'})
title([kLabel{i.US},' N = ',num2str(length(kIdx{i.US}))])

% %% ...generate figure 1c
% figure; hold on
% bar([1],Rho(i.UK,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
% bar([2],Rho(i.UK,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
% bar([3],Rho(i.UK,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
% errorbar([1:3],Rho(i.UK,:),lowerBound(i.UK,:),upperBound(i.UK,:),'k','LineStyle', 'none')
% ylabel('Effect sizes for time of day');ylim([-0.02 0.06])
% set(gca,'xtick',[1 2 3],'xticklabel',{'Gain','Loss','Mix'},'ytick',-0.02:0.02:0.06);xtickangle(45);xlim([0 4])
% axis square
% legend({'Gain','Loss','Mix'})
% title([kLabel{i.UK},' N = ',num2str(length(kIdx{i.UK}))])
% 



save('GBE_Figure1_Data.mat')

end
