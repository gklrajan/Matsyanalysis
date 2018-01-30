%% Figures on/off

set(0, 'DefaultFigureVisible', 'on');

%%

%%
% Data key
% 1 frame #
% 2 camera time in microsec (for data transfer purposes microseconds instead of nanoseconds used)
% 3 system time in ms
% 4 xpos, start at 0
% 5 ypos
% 6 angle in rad
% 7(tap)
clearvars; clc;

[FileName, PathName] = uigetfile('*.bin');
cd(PathName);

% extract parameters from filename
tmp_str = strsplit(FileName, '_');

% save parameters in strings - FORMAT: Date_Time_ExperimentName_Animal#_Remark_Trial#_.bin
acquis_date = tmp_str{1, 1}; acquis_time = tmp_str{1, 2}; exp_type    = tmp_str{1, 3}; fish_num    = tmp_str{1, 4}; trial_num   = tmp_str{1, 6};

% acquisition parameters
window_size         = 0;
num_data_categories = 7+window_size^2; %6 for swim /6+5 for loom
camscale_px_per_mm  = 20.6; % px/mm
datarate_Hz         = 750;  % Hertz
NumVar = 7; %6 for swim /6+5 for loom

% Read and reorganize the bin file
h        = fopen([PathName, FileName]);
tmp_data = fread(h, inf, 'float');
fclose(h);

tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.datarate_Hz = datarate_Hz;

freeSwim.boutAnalysis.boutLength = 150; % considering the bout length to be 200ms / 750Hz acq
freeSwim.boutAnalysis.VthreshON = 1;
freeSwim.boutAnalysis.BthreshOFF = 2;

%%
% CHECK FOR TIMING PROBLEMS AND LOST FRAMES

% TIMER COUNTER: 
% time difference between frames in microseconds, based on the cameras 32bit time stamp counter (25Mhz, 40ns)

time_diff      = [0; diff(tmp_data(:, 2))];                      

% linearize the "saw-tooth function" for the timecounter
idx            = time_diff <= -2^32/1000 + 1.5*median(time_diff); % find the frames when 32bit timecounter was reset to zero
time_diff(idx) = median(time_diff);                               % reset the time difference between these frames to the median in microseconds

% camera frame length in microseconds calculated from the camera timecounter
frame_length_calc_ms = median(time_diff)/1000;

% aquisition datarate in Hertz calculated from the camera timecounter
datarate_calc_Hz = (1/(median(time_diff)).*10^6);  

% CHECK for timing problems (e.g. frames that took unusually long or are
% unusually shorter than what is expected from the used datarate)

idx_time       = abs(time_diff)/1000 >= 1.01*frame_length_calc_ms;  % searches for time differences between frames that are +-1% of the expected frame duration


% FRAME COUNTER:
% index difference between frames, based on the cameras 24bit frame counter

frame_diff = [0; diff(tmp_data(:, 1))]; 

% linearize the "saw-tooth function" for the frame counter (should not
% happen at low datarates) 
idx             = frame_diff <= -2^24 + 1.5*median(frame_diff); % find the frames when 24bit framecounter was reset to zero
frame_diff(idx) = 1;                                            % reset the frame difference between these frames to 1

% CHECK for missing frames
idx_frame  = frame_diff > 1;                         % index of missing frames
idx_lost   = find(idx_frame == 1);                   % first frame in the block of missed frames

% checks if missed timestamps coincide with missed frames, is 0 if inconsistent timestamps outside of missed frames 
isTime = isequal(idx_time, idx_frame);

% calculate total duration of video

duration = tmp_data(end,1)/datarate_Hz/60;

idx_tap = tmp_data(:,7)==1;
idx_tap = find(idx_tap==1);

% prints the above calculated values

fprintf('\nacquisition time: %s %s', acquis_date, acquis_time); 
fprintf('\nexperiment type: %s', exp_type); 
fprintf('\nfish number: %s', fish_num); 
fprintf('\ntrial number: %s', trial_num);
fprintf('\nvideo duration: %2.2f min', duration);
fprintf('\n\nfirst frame in the block of missed frames : number of frames lost\n');
fprintf('\n %d: %d',  [idx_lost, frame_diff(idx_frame)-1].');
fprintf('\n\ntiming flawed (outside of lost frames):  %d  \n', ~isTime  );


%%
% % INSERT nans for lost frames...
% 
% % define anonymous function that inserts (nxm) blocks into (oxm) matrices
insert_blocks = @(insert_block, matrix, n) cat(1,  matrix(1:n-1,:), insert_block, matrix(n:end,:) );

data_raw = tmp_data;

for ii = nnz(idx_frame):-1:1 % starts from the last row in the matrix to keep the indizes for block-insertion
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    tmp_data        = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1) + 1; % framecounter starts at 1


% % read the fish images
% 
% % Read the 3600px associated with each frame and store it row-wise
% tmp_CXP = NaN(size(tmp_data,1),size(tmp_data,2)-NumVar);
% for gg = 1:size(tmp_data,1)
%     pp = 1;
%     for rr = 7:((size(tmp_data,2))-NumVar)
%         tmp_CXP(gg,pp) = tmp_data(gg,rr);
%         pp=pp+1;
%     end
% end
% 
% tmp_CXP = tmp_CXP'; %Transpose to column-wise
% 
% % reshape each 3600 column into a 60x60 px matrix and join each new consecutive
% % frame onto it
% CXPimg = zeros(100,100,(size(tmp_CXP,2)));
% for tt = 1:size(tmp_CXP,2)
%     res = reshape(tmp_CXP(:,tt),100,[])';
%     CXPimg(:,:,tt) = res;
%     f(tt)=im2frame(uint8(CXPimg(:,:,tt)),gray(256)); % create frame for video writing
% end
% figure(3);imshow3D(CXPimg);
% 
% % video conversion
% 
% vidObj = VideoWriter('awesomeDanionella.avi');
% vidObj.FrameRate=20;
% vidObj.Quality=100;
% open(vidObj)
% writeVideo(vidObj,f);
% close(vidObj);

%%
%-------------------------------------------------------------------------------------
%
% IDENTIFICATION OF SWIM BOUTS (bout speeds, IBIs etc)
%
%-------------------------------------------------------------------------------------

xpos = tmp_data(:, 4);
ypos = tmp_data(:, 5);
tmp_ori_deg = rad2deg(tmp_data(:,6));

% plot fish position 
xx = [xpos'; xpos']; tt = [ypos'; ypos'];
z  = 1:size(xpos', 2); zz = [z; z];       % create frame-vector

figure
hs = surf(xx, tt, zz, zz, 'EdgeColor', 'interp', 'LineWidth', 2);
colormap('parula');
view([0 90 0]); 
xlabel('x-position');
ylabel('y-position');
zlabel('frames');

dxV   = [0; diff(xpos)]; % distance between two consecutive x-coordinates
dy   = [0; diff(ypos)]; % distance between two consecutive y-coordinates

tmp_dist_unfilt           = sqrt(dxV.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_Hz;         % convert to velocity in mm/s

%% 
idx_nan     = isnan(dxV);
dxV(idx_nan) = 0; % for filtering nan values need to be removed 
dy(idx_nan) = 0;

% filters used in the analysis
filterB = ones(1,100)./100; %for event detection 100 ms
filterV = ones(1,30)./30; %for more precise onset detection 30 ms
filterF = ones(1,4)./4; %for escapes 4.5 ms

freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.B = filterB;
freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.V = filterV;
freeSwim.boutAnalysis.acquis_date.exp_type.fish_num.filters.F = filterF;


%% orientation

tmp_delta_ori = [nan; diff(tmp_data(:, 6))];

% correction of discontinuties when fish turns from 0 to 2*pi or vice versa

for kk = 1: length(tmp_delta_ori)
    
    if tmp_delta_ori(kk) > pi % right turn
        tmp_delta_ori(kk) =  tmp_delta_ori(kk) - 2*pi;
        
    elseif tmp_delta_ori(kk) < -pi % left turn
        tmp_delta_ori(kk) =  2*pi + tmp_delta_ori(kk);
    end
    
end 

%% Event detection
dxB        = filtfilt(filterB, 1, dxV);
dyB        = filtfilt(filterB, 1, dy);

tmp_dist_fB = sqrt(dxB.^2 + dyB.^2);     % distance moved between iterations, in pixels
tmp_dist_fB = tmp_dist_fB./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fB  = tmp_dist_fB.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fB(idx_nan) = nan; % re-insert the nan values

%% Onset detection
dxV        = filtfilt(filterV, 1, dxV);
dyV        = filtfilt(filterV, 1, dy);

tmp_dist_fV = sqrt(dxV.^2 + dyV.^2);     % distance moved between iterations, in pixels
tmp_dist_fV = tmp_dist_fV./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fV  = tmp_dist_fV.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fV(idx_nan) = nan; % re-insert the nan values

%% Peak detection
dxF        = filtfilt(filterF, 1, dxV);
dyF        = filtfilt(filterF, 1, dy);

tmp_dist_fF = sqrt(dxF.^2 + dyF.^2);     % distance moved between iterations, in pixels
tmp_dist_fF = tmp_dist_fF./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_fF  = tmp_dist_fF.*datarate_Hz;  % convert to velocity in mm/s

tmp_vel_fF(idx_nan) = nan; % re-insert the nan values


%% find peaks
[pks,locs] = findpeaks(tmp_vel_fB,'MinPeakProminence',1,'MinPeakDistance',90); %minPeakProminence & minPeakDist as bout interval and min Velo resp

%% find escape
[escs,fast_locs] = findpeaks(tmp_vel_fF,'MinPeakHeight',90,'MinPeakDistance',90); %minPeakHeight is TH velocity peak for escape

%%
for qq = 1:size(fast_locs)
    
    if (fast_locs(qq) - 3*freeSwim.boutAnalysis.boutLength <=0)... % get rid of any half bout in the begining
                ||(fast_locs(qq) + 3*freeSwim.boutAnalysis.boutLength >=size(tmp_vel_fB,1))... % half bout in the end
                ||(any(isnan(fast_locs(qq):fast_locs(qq)+3*freeSwim.boutAnalysis.boutLength)))... % remove nan-area
                ||(any(isnan(fast_locs(qq):fast_locs(qq)-3*freeSwim.boutAnalysis.boutLength)))... % nan-area again
                ||tmp_vel_fB(fast_locs(qq))<= freeSwim.boutAnalysis.VthreshON... %check w/ Vthresh
                ||tmp_vel_fB(fast_locs(qq))<= freeSwim.boutAnalysis.BthreshOFF %check w/ Bthresh
            
            fast_locs(qq)= NaN;
    end
end


%%
%to eliminate escape-like peaks unassociated with tap events
    ttt=1;
for hhh = 1:size(fast_locs)
    
    if ((fast_locs(hhh)-idx_tap(ttt)<=3760)... %tap more than 5 secs before escape
        && (fast_locs(hhh)-idx_tap(ttt)>0)) %negative i.e. tap is after swim (next tap event)    
    ttt = ttt+1;
    
    else
    fast_locs(hhh)= NaN;
    
    end
end

%% ID swim bouts and operate on them

%bout key
%1 start frame
%2 end frame
%3 bout duration in ms
%4 orientation before bout in deg
%5 orientation after bout in deg
%6 turn in deg

%7 mean velo
%8 peak velo
%9 total distance
%10 total L yaw in deg
%11 total R yaw in deg
%12 angular velocity in deg/sec

%13 inter-bout-interval in ms

fast_locs(isnan(fast_locs)) = []; %imp - esle, error w/ for loop - sunbscript indices must be real integers or logicals

tmp_swim_bouts = zeros(size(fast_locs,1),15);

for mm = 1:size(fast_locs,1)
   
    %find start frame
    sss = 0;
    while fast_locs(mm)-sss>=2*freeSwim.boutAnalysis.boutLength...
            && tmp_vel_fV(fast_locs(mm)-sss) >= freeSwim.boutAnalysis.VthreshON % using narrow filter for onset
        sss = sss + 1;
    end   
    
    tmp_swim_bouts(mm,1) = fast_locs(mm)-sss+1; %bout start frame
    
    
    %find end frame
    eee = 0;
    while fast_locs(mm)<=size(tmp_vel_fB,1)...
            && tmp_vel_fB(fast_locs(mm)+eee) >= freeSwim.boutAnalysis.BthreshOFF % using broad filter for offset
        eee = eee + 1;
    end
    
    tmp_swim_bouts(mm,2) = fast_locs(mm)+eee-1; %bout end frame
    
    %bout duration in ms
    tmp_swim_bouts(mm,3) = (tmp_swim_bouts(mm,2) - tmp_swim_bouts(mm,1))*frame_length_calc_ms;
    
    %orientation before bout /avg of 10 frames or 13.3 ms
    tmp_swim_bouts(mm,4) = nanmean(tmp_ori_deg((tmp_swim_bouts(mm,1)-11):(tmp_swim_bouts(mm,1)-1)));
    
    %orientation after bout /avg of 10 frames or 13.3 ms
    tmp_swim_bouts(mm,5) = nanmean(tmp_ori_deg((tmp_swim_bouts(mm,2)+1):(tmp_swim_bouts(mm,2)+11)));
        
    %TURN (ori after - ori before)
    tmp_swim_bouts(mm,6) = tmp_swim_bouts(mm,5)-tmp_swim_bouts(mm,4);
    
    %correct for change in orientation
    if tmp_swim_bouts(mm,6) > pi % right turn 
        tmp_swim_bouts(mm,6) =  tmp_swim_bouts(mm,6) - 2*pi;
        
    elseif tmp_swim_bouts(mm,6) < -pi % left turn
        tmp_swim_bouts(mm,6) =  2*pi + tmp_swim_bouts(mm,6);
    end
    
    % mean velocity during bout
    tmp_swim_bouts(mm,7) = nanmean(tmp_vel_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    % max velocity during bout
    tmp_swim_bouts(mm,8) = nanmax(tmp_vel_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    % total distance covered during the bout
    tmp_swim_bouts(mm,9) = nansum(tmp_dist_fV(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2)));
    
    %head yaw during bouts
    yaws = tmp_delta_ori(tmp_swim_bouts(mm,1):tmp_swim_bouts(mm,2));
    tmp_swim_bouts(mm,10)= nansum(yaws(yaws>0)); %summation of total left yaws
    tmp_swim_bouts(mm,11)= nansum(yaws(yaws<0)); %summation of total right yaws
    
    %angular velocity
    tmp_swim_bouts(mm,12) = nansum(abs(yaws))/tmp_swim_bouts(mm,3); %deg/sec
    
    %delay
    tmp_swim_bouts(mm,13) = (tmp_swim_bouts(mm,1)-idx_tap(mm))*frame_length_calc_ms; %if there is an error --> error in tap/excape ID
end

 
% %% IBI section
% tmp_vel_IBI_fV = tmp_vel_fV;
% 
% % ID overlapping bouts and NaN them
% % very poorly written; had an indexing error and then just worked
% % around making it poorer and poorer, but it works!
% % changed to more efficient one on 22/01/18
% 
% for oo = 1:size(tmp_swim_bouts,1)
%       if (oo>=2 && tmp_swim_bouts(oo,1)<=tmp_swim_bouts(oo-1,2))...
%            || (oo<=size(tmp_swim_bouts,1)-1 && tmp_swim_bouts(oo,2)>=tmp_swim_bouts(oo+1,1))
%        
%        tmp_swim_bouts(oo, 3:end) = NaN; % convert all bout parameters corresp to this bout to NaN;
%        locs(oo) = NaN;
%        tmp_vel_IBI_fV(tmp_swim_bouts(oo,1):tmp_swim_bouts(oo,2))= NaN; % the velo vector to NaN as well - imp to ID NaN values in IBI cals  
%        
%       end
% end
% 
%             tmp_swim_bouts(isnan(tmp_swim_bouts(:,3)), : ) = [];
%             locs(isnan(locs))                              = [];
%                        
% % ID consecutive NaNs in velocity vector - imp for IBI identification
% tmp_vel_CC = bwconncomp(isnan(tmp_vel_IBI_fV));
% pixel_Idx = cellfun('prodofsize',tmp_vel_CC.PixelIdxList);
% tmp_vel_NaNs = zeros(size(tmp_vel_IBI_fV));
% 
% for zzz = 1:tmp_vel_CC.NumObjects
%     tmp_vel_NaNs(tmp_vel_CC.PixelIdxList{zzz})=pixel_Idx(zzz);
% end 
% 
% % Calculate the IBI
% tmp_swim_bouts(:,13) = NaN;
% 
% for ibi=2:size(tmp_swim_bouts,1)
%    tmp_swim_bouts(ibi,13)= (tmp_swim_bouts(ibi,1)-tmp_swim_bouts(ibi-1,2))*frame_length_calc_ms;
%    
%    if tmp_vel_NaNs(ibi)>=1.5*freeSwim.boutAnalysis.boutLength % >300 ms
%        tmp_swim_bouts(ibi,13)=NaN;
%    end
%     
% end


%% some plots
fig1 = figure;
hold on; 
plot(tmp_vel_unfilt); 
%plot(tmp_vel_fV,'LineWidth', 3);
plot(tmp_vel_fB,'LineWidth', 3);
%plot(tmp_vel_fF,'LineWidth', 3);
%plot(tmp_vel_fF,'LineWidth', 3);                

%plot(20*pks);
%plot(20*dxB-10);
%plot(dxV-10);
%plot(20*dyB-20);
%plot(dyV-20);
plot(50*tmp_delta_ori-3);
vline(tmp_swim_bouts(:,1),'m');
vline(tmp_swim_bouts(:,2),'m');
vline(idx_tap+1880,'k');
%vline(fast_locs,'c');
hold off;

% fig2 = figure; histfit(tmp_swim_bouts(:,3),50);
% fig3 = figure; histogram(tmp_swim_bouts(:,7),50);
% fig4 = figure; histogram(tmp_swim_bouts(:,8),50);
% fig5 = figure; histogram(tmp_swim_bouts(:,12),50);
% fig6 = figure; histogram(tmp_swim_bouts(:,13),50);
% fig7 = figure; histfit(tmp_swim_bouts(:,6),50);
% fig8 = figure; histfit(tmp_swim_bouts(:,12),50);

 dcm1 = datacursormode(fig1);
 set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');

%save('/Institut Curie/Lab/Projects/Scripts/ZebranalysisSystem/AST_TAP_data.mat',freeSwim);