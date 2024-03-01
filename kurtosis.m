% kurtosis.m 
%
% This code takes a PDF and calculates the relaxation advection using the
% following formula: D = lambda*I, v(x) = (lambda/p(x))*grad(p(x)). The code then
% propagates a uniform PDF using the calculated advection and
% diffusion. Example problems are 2D gaussians and s-order super-gaussians,
% where the advection is known analytically. Other examples are numerically
% approximated, via convex/non-convex shapes of datasets and KDE. 
%
% By Benjamin L. Hanson, 2024

clear all; close all; clc; 

beta = linspace(1,50,1000);
g2 = gamma2(beta);

figure(1); clf; hold on;
xlim([beta(1),beta(end)])
xlabel('$\beta$','Interpreter','latex','FontSize',16)
ylabel('$\gamma_2(\beta)$','Interpreter','latex','FontSize',16)
plot(beta,g2,'k-','LineWidth',2)

beta = linspace(0.1,1,1000);
g2 = gamma2(beta);

figure(2); clf; hold on;
xlim([beta(1),beta(end)])
xlabel('$\beta$','Interpreter','latex','FontSize',16)
ylabel('$\gamma_2(\beta)$','Interpreter','latex','FontSize',16)
plot(beta,g2,'k-','LineWidth',2)
set(gca,'YScale','log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g2 = gamma2(beta)
    g2 = NaN(1,length(beta));
    for i = 1:length(beta)
        g2(i) = ((gamma(5/(2*beta(i)))*gamma(1/(2*beta(i))))/(gamma(3/(2*beta(i))))^2)-3;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%