% GGDAF_demo.m 
%
% This code plots three 2D Generalized Gaussian Distributions with Anisotropic
% Flatness, one where the shaping parameter is the same along the
% eigenvectors, and the other two where the shaping parameter is different.
% It also demonstrates that the area under the curve, and first and second
% central moments are as expected.
%
% By Benjamin L. Hanson, Feb 23 2024

clear all; close all; clc;

mu1 = 0; mu2 = 0; mu = [mu1; mu2]; sigma1 = 1; sigma2 = 1; Sigma = [sigma1^2 0; 0 sigma2^2];
b1 = [2 2]; b2 = [1 3]; b3 = [2 0.5]; 

dx = 0.01; dy = 0.01;
x = -2:dx:2;
y = -2:dy:2;
[X,Y] = meshgrid(x,y);
P1 = GGDAF(X,Y,mu,Sigma,b1);
P2 = GGDAF(X,Y,mu,Sigma,b2);
P3 = GGDAF(X,Y,mu,Sigma,b3);

fig = figure(1); clf; hold on; cb = colorbar(); axis square; fig.Position = [50 300 500 500];
xlim([x(1),x(end)])
ylim([y(1),y(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
contour(X,Y,P1,[linspace(0.05,max(P1,[],'all'),10)], 'LineWidth',2);
tightfig; 

fig = figure(2); clf; hold on; cb = colorbar(); axis square; fig.Position = [550 300 500 500];
xlim([x(1),x(end)])
ylim([y(1),y(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
contour(X,Y,P2,[linspace(0.05,max(P2,[],'all'),10)], 'LineWidth',2);
tightfig;

fig = figure(3); clf; hold on; cb = colorbar(); axis square; fig.Position = [1050 300 500 500];
xlim([x(1),x(end)])
ylim([y(1),y(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
contour(X,Y,P3,[linspace(0.05,max(P3,[],'all'),10)], 'LineWidth',2);
tightfig; 

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')

fprintf('GGDAF #1 \n')
vol = integral2(@(x,y) GGDAF(x,y,mu,Sigma,b1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2),

m1 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b1,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
m2 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b1,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
mean = [m1;m2],

s11 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b1,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s12 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b1,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s21 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b1,3),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s22 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b1,4),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
cov = [s11 s12; s21 s22],

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')

fprintf('GGDAF #2 \n')
vol = integral2(@(x,y) GGDAF(x,y,mu,Sigma,b2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2),

m1 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b2,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
m2 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b2,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
mean = [m1;m2],

s11 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b2,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s12 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b2,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s21 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b2,3),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s22 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b2,4),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
cov = [s11 s12; s21 s22],

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')

fprintf('GGDAF #3 \n')
vol = integral2(@(x,y) GGDAF(x,y,mu,Sigma,b3),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2),

m1 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b3,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
m2 = integral2(@(x,y) GGDAFM1(x,y,mu,Sigma,b3,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
mean = [m1;m2],

s11 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b3,1),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s12 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b3,2),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s21 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b3,3),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
s22 = integral2(@(x,y) GGDAFM2(x,y,mu,Sigma,b3,4),-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2);
cov = [s11 s12; s21 s22],

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAF(X,Y,mu,Sigma,beta)
    d = length(mu);
    if all(beta == beta(1))
        same = true; 
        Ai = A(beta(1),d,Sigma);
    else
        same = false; 
        Ai = (det(Sigma))^(1/2);
        for i = 1:length(mu)
            Ai = Ai*A(beta(i),1,Sigma);
        end
    end
    [S,D] = eig(Sigma); L = sqrtm(inv(D));
    
    P = NaN(size(X)); s = size(X); row = s(1); col = s(2); 
    for i = 1:row
        for j = 1:col
            x = [X(i,j); Y(i,j)];
            q = L*S'*(x-mu);
            sum = 0;
            for k = 1:length(mu)
                sum = sum + (Bt(beta(k),d,same)*q(k)^2)^beta(k); 
            end
            sum2 = 0;
            for k = 1:beta(1)-1
                sum2 = sum2 + nchoosek(beta(1),k)*(q(1)^2)^(beta(1)-k)*(q(2)^2)^k;
            end
            kd = beta(1)==beta(2);
            P(i,j) = Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAFM1(X,Y,mu,Sigma,beta,comp)
    d = length(mu);
    if all(beta == beta(1))
        same = true; 
        Ai = A(beta(1),d,Sigma);
    else
        same = false; 
        Ai = (det(Sigma))^(1/2);
        for i = 1:length(mu)
            Ai = Ai*A(beta(i),1,Sigma);
        end
    end
    [S,D] = eig(Sigma); L = sqrtm(inv(D));
    
    P = NaN(size(X)); s = size(X); row = s(1); col = s(2); 
    for i = 1:row
        for j = 1:col
            x = [X(i,j); Y(i,j)];
            q = L*S'*(x-mu);
            sum = 0;
            for k = 1:length(mu)
                sum = sum + (Bt(beta(k),d,same)*q(k)^2)^beta(k); 
            end
            sum2 = 0;
            for k = 1:beta(1)-1
                sum2 = sum2 + nchoosek(beta(1),k)*(q(1)^2)^(beta(1)-k)*(q(2)^2)^k;
            end
            kd = beta(1)==beta(2);
            P(i,j) = (x(comp)-mu(comp))*Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAFM2(X,Y,mu,Sigma,beta,comp)
    d = length(mu);
    if all(beta == beta(1))
        same = true; 
        Ai = A(beta(1),d,Sigma);
    else
        same = false; 
        Ai = (det(Sigma))^(1/2);
        for i = 1:length(mu)
            Ai = Ai*A(beta(i),1,Sigma);
        end
    end
    [S,D] = eig(Sigma); L = sqrtm(inv(D));
    
    P = NaN(size(X)); s = size(X); row = s(1); col = s(2); 
    for i = 1:row
        for j = 1:col
            x = [X(i,j); Y(i,j)];
            q = L*S'*(x-mu);
            sum = 0;
            for k = 1:length(mu)
                sum = sum + (Bt(beta(k),d,same)*q(k)^2)^beta(k); 
            end
            sum2 = 0;
            for k = 1:beta(1)-1
                sum2 = sum2 + nchoosek(beta(1),k)*(q(1)^2)^(beta(1)-k)*(q(2)^2)^k;
            end
            kd = beta(1)==beta(2);
            switch comp
                case 1
                    P(i,j) = (x(1)-mu(1))^2*Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
                case 2
                    P(i,j) = (x(1)-mu(1))*(x(2)-mu(2))*Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
                case 3
                    P(i,j) = (x(1)-mu(1))*(x(2)-mu(2))*Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
                case 4
                    P(i,j) = (x(2)-mu(2))^2*Ai*exp(-(sum + B(beta(1),d)^(beta(1))*kd*sum2));
                otherwise
                    error('Incorrect component number.')
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = A(beta,d,Sigma)
    f = (B(beta,d)/pi)^(d/2)*(gamma(d/2)*beta)/(gamma(d/(2*beta))*det(Sigma)^(1/2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Bt(beta,d,same)
    if ~same
        f = B(beta,1);
    else
        f = B(beta,d);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = B(beta,d)
    f = gamma((d+2)/(2*beta))/(d*gamma(d/(2*beta)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%