% example code to visualize FFT solution

clear; clc; close all

N = 200;

n_plot = [1,10,100];

x_vec = linspace(0,1,201);

theta_N = 0;

for n = 1:N

    l_n = n*pi;
    si_n = sqrt(2)/l_n;
    phi_n = sqrt(2)*sin(l_n*x_vec);

    theta_N = theta_N + si_n*phi_n; % sum

    if ismember(n,n_plot)
        plot(x_vec,theta_N); hold on

    end


end

legend({'N=1','N=10','N=100'})
ylim([0 1.2])
