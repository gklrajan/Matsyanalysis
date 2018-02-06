clearvars; clc;

load('Danionella_040218.mat');
pooledRaw_Danio =[freeSwim(1).boutData;freeSwim(2).boutData;freeSwim(3).boutData;freeSwim(4).boutData;freeSwim(5).boutData;freeSwim(6).boutData;freeSwim(7).boutData;freeSwim(8).boutData;freeSwim(9).boutData];
pooledMedian_MT = nanmean(pooledRaw_Danio);
pooledMAD_MT = std(pooledRaw_Danio);
close();

load('Danionella5dpf_050218.mat');
pooledRaw_WT =[freeSwim(1).boutData;freeSwim(2).boutData; freeSwim(3).boutData;freeSwim(4).boutData];
pooledMedian_WT = nanmean(pooledRaw_WT);
pooledMAD_WT = std(pooledRaw_WT);
close();

for ii=1:size(pooledRaw_Danio,2)
[p(ii),h(ii),stats(ii)] = ranksum(pooledRaw_Danio(:,ii),pooledRaw_WT(:,ii),'alpha',0.01);
end

% for ii=1:size(pooledRaw_MT,2)
% [h(ii),p(ii)] = ttest2(pooledRaw_MT(:,ii),pooledRaw_WT(:,ii));
% end

pooledRaw_Danio(:,12)=rad2deg(pooledRaw_Danio(:,12));
pooledRaw_WT(:,12)=rad2deg(pooledRaw_WT(:,12));

jugaad = [3,6,7,8,9,12,13,16];
title_jugaad = ["bout duration in ms","turn angle in deg","mean velocity in mm/sec","peak velocity in mm/sec","total distance in mm","angular velocity in deg/sec","inter-bout interval in ms","number of tail beats per bout"];
xlabel_jugaad = ["duration (ms)","turn (deg)","mean velocity(mm/sec)","peak velocity (mm/sec)","total distance (mm)","angular velocity (deg/sec)","inter-bout interval (ms)","number of tail beats/bout"];
figure1 = figure;
for gg = 1:length(jugaad)
subplot(2,8,gg);
histogram(pooledRaw_WT(:,jugaad(gg)),70);        
title(title_jugaad(gg));
xlabel(xlabel_jugaad(gg));
ylabel({'# of bouts'});

subplot(2,8,gg+8);
histogram(pooledRaw_Danio(:,jugaad(gg)),70);
title(title_jugaad(gg));
xlabel(xlabel_jugaad(gg));
ylabel({'# of bouts'});
end

% %to plot overlayed histograms of WT & mutant
% for kk = 1:length(jugaad)
% subplot(3,3,kk);
% histogram(pooledRaw_WT(:,jugaad(kk)),100,'FaceColor','g','FaceAlpha',1);
% hold on;
% histogram(pooledRaw_Danio(:,jugaad(kk)),100,'FaceColor','r','FaceAlpha',0.5)
% hold off;
% end