clearvars; clc;

load('AST_TAP_data.mat');
pooledRaw_MT =[tapEscape(1).boutData;tapEscape(2).boutData];
pooledMean_MT = nanmean(pooledRaw_MT);
close();

load('WT_TAP_data.mat');
pooledRaw_WT =[tapEscape(1).boutData;tapEscape(2).boutData; tapEscape(3).boutData];
pooledMean_WT = nanmean(pooledRaw_WT);
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