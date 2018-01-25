clearvars; clc;

load('TEST_data.mat');
pooledRaw_TEST =[freeSwim(1).boutData;freeSwim(2).boutData];
pooledMean_TEST = nanmean(pooledRaw_TEST);
pooledMedian_TEST = nanmedian(pooledRaw_TEST);
close();

histogram(pooledRaw_TEST(:,3),100);