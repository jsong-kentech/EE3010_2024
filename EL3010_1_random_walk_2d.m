%% Random walk 
% 2024.02.25
% Ju Song
% KENTECH

% Contents
% (1) Config 
% (2) 2D Random walk (to visualize)
% (3) 1D Random walk (to analyze)

clear; clc; close all
%% Config
dx_size = 1;
dy_size = 1;
dt = 1;

N = 4; % number of random walkers
t_end = 10000; % time period to simulate
t_vec = 0:dt:t_end;
M = length(t_vec); % number of time steps

dt_pause = 0.1; % animation duration (sec)

%% 2D Random walk

% initialize 
x_mat = zeros(N,M);
y_mat = zeros(N,M);


% time step loop
for m = 1:M
    
    % generate random walk for N particles
    dx_now = (2*unidrnd(2,N,1)-3)*dx_size;
    dy_now = (2*unidrnd(2,N,1)-3)*dy_size;
    
    % figure(1)
    % plot(dx_vec,'ob'); hold on 
    % plot(dy_vec,'sqr')

    % calculate position
    if m ==1
       x_mat(:,m) = zeros(N,1);
       y_mat(:,m) = zeros(N,1);
    else
        x_mat(:,m) = x_mat(:,m-1) + dx_now;
        y_mat(:,m) = y_mat(:,m-1) + dy_now;
    end
end


% visualize each path
c_mat = lines(N);
for n = 1:N
    figure(2); hold on
    plot(x_mat(n,:),y_mat(n,:),'Color',c_mat(n,:))

end

% same plot
% figure(3)
% plot(x_mat',y_mat')


% animation 
for m = 1:M

    % each time, draw up to t
    for n = 1:N
    plot(x_mat(n,1:m),y_mat(n,1:m),'-','Color',c_mat(n,:)); hold on
    plot(x_mat(n,m),y_mat(n,m),'o','MarkerFaceColor',c_mat(n,:))
    end
    hold off
    pause(dt_pause)

end
