close all
phi = [10 100];

x_vec = 0:0.01:1;



c_mat = lines(length(phi));

for i = 1:length(phi)


c_vec = cosh(phi(i)*x_vec) - tanh(phi(i))*sinh(phi(i)*x_vec);

c_vec2 = exp(-phi(i)*x_vec);

figure(1)
plot(x_vec,c_vec,'-','color',c_mat(i,:)); hold on
plot(x_vec,c_vec2,'--','color',c_mat(i,:),'LineWidth', 3); hold on
%ylim([0 1])


end
