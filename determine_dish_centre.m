function [dish_centre, flag] = determine_dish_centre(snap, inmiddlethr, date, exp_identifier, fish_no, flag_key);

figure(1), imshow(snap);

for j = 1:2
    for i = 1:2
        shg
        dcm_obj = datacursormode(1);
        set(dcm_obj,'DisplayStyle','window',...
            'SnapToDataVertex','off','Enable','on')
        waitforbuttonpress
        c_info{i,j} = getCursorInfo(dcm_obj);
        
        points{i,j} = c_info{i,j}.Position;
    end
end

close(1);

points = cell2mat(reshape(points, [4,1]));
P1=(points(1,:))'; P2=(points(2,:))'; P3=(points(3,:))'; P4=(points(4,:))';
P1P2 = (P2-P1); P3P4 = (P4-P3);

% equations:
% x_1(:,ii) = P1 + s(ii)* P1P2;
% x_2(:,ii) = P3 + t(ii)* P3P4;

% find vectors perpendicular to respective direction vectors

orth_x12 = [1 ; -P1P2(1)/P1P2(2)];
orth_x34 = [1 ; -P3P4(1)/P3P4(2)];

P0_12    = P1 + ( (P2-P1) * 0.5 );
P0_34    = P3 + ( (P4-P3) * 0.5 );

% equations:
% x_012(:,ii) = P0_12 + u(ii)* orth_x12;
% x_034(:,ii) = P0_34 + v(ii)* orth_x34;

% find intersection of the lines by solving the respective linear equation
% system u_v * M = P0_12 - P0_34

M           = [ -orth_x12 , orth_x34] ;
u_v         = inv(M) * (P0_12 - P0_34);
dish_centre = round(u_v(1)*orth_x12 + P0_12);

%clearvars -except camscale data datarate dish_centre exp_identifier fish_no fname no_unique_abl_sites free_swim_abl path_abl path_FS snap

% radius_eff  = 392; %round(32.5/2*camscale); % 27mm/2 = inner radius on the bottom of the dish
t           = linspace(0,2*pi);

figure(2); imshow(snap);
% hold on,
% plot(dish_centre(1), dish_centre(2),'x','markers',10, 'MarkerEdgeColor','k');
% plot(dish_centre(1) + inmiddlethr*cos(t), dish_centre(2) + inmiddlethr*sin(t), 'LineWidth', 2);
% plot(dish_centre(1) + 430*cos(t), dish_centre(2) + 430*sin(t), 'LineWidth', 2, 'Color', [1 1 1]); %421 should depend on scale of the camera
% % text(200, 500 , {['ablation date: ' date];...
% %     ['experiment: A' num2str(exp_identifier) ];...
% %     ['fish-no:  ' num2str(fish_no) ];...
% %    }) 
% text(700, 500 , ['dish-centre: ' '[' num2str(dish_centre(1)) ',' num2str(dish_centre(2)) ']' ]) 
% %text(800, 900 , ['to flag that fish press ' flag_key] )
% 
% waitforbuttonpress
% 
% % flag dish centres that were not fitted well
% flag = strcmp(get(gcf,'currentkey'), flag_key);

close all

end
