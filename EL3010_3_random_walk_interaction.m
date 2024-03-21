


clear;clc; close all


%% Config
dx =1;
dt =1;

N = 10000;
t_end = 1000; % time period to simulate
t_vec = 0:dt:t_end;
M = length(t_vec);

L = 10; % L-by-L 2D space

x_mat = zeros(N,M);


c_max = 3000;

%% INITIAL POSITIONS
x_mat(:,1) = unidrnd(L,N,1);

px = zeros(N,1);
for m = 2:M

    %%  Random Walks
        % generate random walk
        
        % CASE 1
        % dx_now = (2*unidrnd(2,N,1)-3)*dx;

        % CASE 2
%         x_vec = x_mat(:,m-1);
%         c = histcounts(x_vec,1:L+1);
%         c_extended = [c(1) c c(end)];
%         
%         px(c_extended(x_vec) < c_extended(x_vec+2)) = 0.75;
%         px(c_extended(x_vec) == c_extended(x_vec+2)) = 0.5;
%         px(c_extended(x_vec) > c_extended(x_vec+2)) = 0.25;
% 
%         px_now = px;
%         dx_now = 2*(rand([N,1]) < px_now)-1;
% 

        % CASE 3
        x_prev = x_mat(:,m-1);
        c = histcounts(x_prev,1:L+1);

        c_0 = nan(size(x_prev));
        c_p = nan(size(x_prev));
        c_n = nan(size(x_prev));
        p_p = nan(size(x_prev));
        p_0 = nan(size(x_prev));
        dx_now = nan(size(x_prev));


        c_0 = c(x_prev)';


        % not x = 1
        i_not_1 = (x_prev >1);
        c_n(i_not_1) = c(x_prev(i_not_1)-1); % at 1, NaN (leftover)
        
        % not x = L
        i_not_L =(x_prev <L);
        c_p(i_not_L) = c(x_prev(i_not_L)+1); % at L, NaN  (leftover)


        % mid
        i_bulk = (x_prev <L) & (x_prev >1);
        p_p(i_bulk) = c_p(i_bulk)./(c_p(i_bulk)+c_n(i_bulk));
        p_p(i_bulk & (c_p>c_max)) = 0;
        p_p(i_bulk & (c_n>c_max)) = 1;
        dx_now(i_bulk) = 2*(rand([nnz(i_bulk),1]) < p_p(i_bulk))-1;

        % x=1
        i_1 =(x_prev ==1);
        p_p(i_1) = c_p(i_1)./(c_p(i_1)+c_0(i_1));
        p_p(i_1 & (c_0>c_max)) = 1-c_max./c_0(i_1 & (c_0>c_max));
        dx_now(i_1) = (rand([nnz(i_1),1]) < p_p(i_1));

        % x=L
        i_L =(x_prev ==L);
        p_0(i_L) = c_0(i_L)./(c_0(i_L)+c_n(i_L));
        p_0(i_L & (c_0>c_max)) = c_max./c_0(i_L & (c_0>c_max));
        dx_now(i_L) = (rand([nnz(i_L),1]) < p_0(i_L))-1;

        % 
        if nnz(isnan(dx_now)) > 0
            error('random step of some particles are not calculated.')
        end

%         % crowding
%         px(c_max < c_extended(x_prev+2)) = 0;
%         px(c_max < c_extended(x_prev)) = 1;
%  
%              % edge
%              i_ub =(c_max < c_extended(x_prev+2)) & (x_prev ==L)';
%              px(i_ub) = c_max./c_extended(x_prev(i_ub));
%              i_lb =(c_max < c_extended(x_prev)) & (x_prev ==1)';
%              px(i_lb) = 1-c_max./c_extended(x_prev(i_lb))';


        % calculate position
        x_mat(:,m) = x_mat(:,m-1) + dx_now;

        % bounded space
        i_x_lb = x_mat(:,m) < 1;
        i_x_ub = x_mat(:,m) > L;

        x_mat(i_x_lb,m) = 1;
        x_mat(i_x_ub,m) = L;




end

%% Result Plot

for m = 1:1000

figure(1)


c = histcounts(x_mat(:,m),1:L+1);
plot(c); hold on
% c = histcounts(x_mat(:,m-1),1:L+1);
% plot(c); 
% c = histcounts(x_mat(:,m-2),1:L+1);
% plot(c); hold off
%ylim([0 N/L*2])
pause(0.021)

end
return

