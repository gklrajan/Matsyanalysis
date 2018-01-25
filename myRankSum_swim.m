clearvars; clc;

load('AST_swim.mat');
pooledRaw_MT =[freeSwim(1).boutData;freeSwim(2).boutData];
pooledMedian_MT = nanmedian(pooledRaw_MT);
pooledMAD_MT = mad(pooledRaw_MT);
close();

load('WT_swim.mat');
pooledRaw_WT =[freeSwim(1).boutData;freeSwim(2).boutData; freeSwim(3).boutData];
pooledMedian_WT = nanmedian(pooledRaw_WT);
pooledMAD_WT = mad(pooledRaw_WT);
close();

for ii=1:size(pooledRaw_MT,2)
[p(ii),h(ii),stats(ii)] = ranksum(pooledRaw_MT(:,ii),pooledRaw_WT(:,ii),'alpha',0.01);
end

% for ii=1:size(pooledRaw_MT,2)
% [h(ii),p(ii)] = ttest2(pooledRaw_MT(:,ii),pooledRaw_WT(:,ii));
% end

jugaad = [3,6,7,8,9,12,13];

figure1 = figure;
for gg = 1:length(jugaad)
subplot(3,3,gg);
histogram(pooledRaw_WT(:,jugaad(gg)),100);
end

figure2 = figure;
for gg = 1:length(jugaad)
subplot(3,3,gg);
histogram(pooledRaw_MT(:,jugaad(gg)),100);
end
