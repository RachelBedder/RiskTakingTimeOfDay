%% Generate all the statistics and plots for the Great Brain Experiment Time of Day Paper

%Rachel Bedder 2020 (rachel.bedder.15@ucl.ac.uk) 

%To run the script from scratch, load the following and run the code block below.
load('UKUSdata.mat'); 

%If you do not want to run the script but just look at the statistics tables that correspond to the paper, load
%the following, this is with 10,000 iterations for bootstrapping and
%permutation tests
load('GBE_PaperStatistics.mat');

%Or if you want to look at the data generated by each figure function, load the
%following. This includes the statistics structure for this function only.
load('GBE_Figure1_Data.mat')
%or
load('GBE_Figure2_Data.mat')
%or
load('GBE_Figure3_Data.mat')

%% ...set the number of permutations for the p statistics and bootstrapping tests
nSample                  =       [10];

%%...generate figures
[statistics.fig1] = GBE_Figure1(mData,nSample);  %...make figure 1: Model Free Effect Sizes for Gambling and Time of Day
[statistics.fig2] = GBE_Figure2(mData,nSample);  %...make figure 2: Model Based Effect Sizes for Prospect Theory Betas and Time of Day
[statistics.fig3] = GBE_Figure3(mData,nSample);  %...make figure 2: Within Subjects model based and model free
