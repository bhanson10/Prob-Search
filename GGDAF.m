% GGDAF.m 
%
% This code plots two Generalized Gaussian Distributions with Anisotropic
% Flatness, one where the shaping parameter is the same along the
% eigenvectors, and the other where the shaping parameter is different.
%
% By Benjamin L. Hanson, Feb 23 2024

clear all; close all; clc;

mu = [0; 0]; Sigma = [1 0; 0 1];
[S,D] = eig(Sigma); L = sqrtm(D);
b1 = [2 2]; b2 = [1 3]; 

dx = 0.01; dy = 0.01;
x = -2:dx:2;
y = -2:dy:2;
[X,Y] = meshgrid(x,y);

P1 = NaN(size(X)); P2 = NaN(size(X));
rc = 1;
for i=x(end):-dx:x(1)
    cc = 1;
    for j=y(1):dy:y(end)
        P1(rc,cc) = pGGDAF(i,j,L,S,mu,b1);
        P2(rc,cc) = pGGDAF(i,j,L,S,mu,b2);
        cc = cc + 1;
    end
    rc = rc + 1; 
    
end
P1 = P1./(dx*dy*sum(P1,'all'));
P2 = P2./(dx*dy*sum(P2,'all'));

figure(1); clf; hold on; cb = colorbar(); axis square; 
xlim([x(1),x(end)])
ylim([y(1),y(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
contour(X,Y,P1,[linspace(0.01,max(P1,[],'all'),10)], 'LineWidth',2);

figure(2); clf; hold on; cb = colorbar(); axis square; 
xlim([x(1),x(end)])
ylim([y(1),y(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
contour(X,Y,P2,[linspace(0.01,max(P2,[],'all'),10)], 'LineWidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = pGGDAF(x,y,L,S,mu,b)
    bs = max(b); 
    q = L*S'*([x;y]-mu);
    sum = 0;
    for i = 1:length(mu)
        sum = sum + (q(i)^2)^b(i); 
    end
    sum2 = 0;
    for i = 1:(floor(bs)-1)
        sum2 = sum2 + nchoosek(bs,i)*(q(1)^2)^(bs-i)*(q(2)^2)^(i); 
    end
    c = max((b(1)-bs+1)*(b(2)-bs+1),0); 
    sum = sum + c*sum2;
    p = exp(-sum);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%