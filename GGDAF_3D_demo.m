% GGDAF_2D_demo.m 
%
% This code plots three 2D Generalized Gaussian Distributions with Anisotropic
% Flatness, one where the shaping parameter is the same along the
% eigenvectors, and the other two where the shaping parameter is different.
% It also demonstrates that the area under the curve, and first and second
% central moments are as expected.
%
% By Benjamin L. Hanson, 2024

clear all; clc; close all; 

mu1 = 0; mu2 = 0; mu3 = 0; p.mu = [mu1; mu2; mu3]; 
sigma1 = 1; sigma2 = 1; sigma3 = 1; p.Sigma = [sigma1^2 0 0; 0 sigma2^2 0; 0 0 sigma3^3]; 
p.b = 1; p.bt = [2 1 1]; 

dx = 0.1; dy = 0.1; dz = 0.1;
x = -5:dx:5; y = -5:dy:5; z = -5:dy:5;
[X,Y,Z] = meshgrid(x,y,z);
s = size(X); P1 = NaN(s); P2 = NaN(s); row = s(1); col = s(2); lay = s(3);
for i = 1:row
    for j = 1:col
        for k = 1:lay
            P1(i,j,k) = GGD(X(i,j,k),Y(i,j,k),Z(i,j,k),p);
            P2(i,j,k) = GGDAF(X(i,j,k),Y(i,j,k),Z(i,j,k),p);
        end
    end
end

fig = figure(1); clf; hold on; cb = colorbar(); axis square; fig.Position = [150 300 600 500];
view(-109,14); lighting phong; light('Position',[-1 0 0]); 
title(['GGD, $\beta$ = ',num2str(p.b)],'Interpreter','latex','FontSize',16)
xlim([x(1),x(end)])
ylim([y(1),y(end)])
zlim([z(1),z(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
zlabel('z','Interpreter','latex','FontSize',16)
isosurface(X,Y,Z,P1,0.01, "HandleVisibility","off");
isosurface(X,Y,Z,P1,0.001, "HandleVisibility","off"); alpha(.5);
isosurface(X,Y,Z,P1,0.0001, "HandleVisibility","off"); alpha(.3);
colormap(cool);
drawnow;
%tightfig; 

fig = figure(2); clf; hold on; cb = colorbar(); axis square; fig.Position = [750 300 600 500];
view(-109,14); lighting phong; light('Position',[-1 0 0]); 
title(['GGDAF, $\beta$ = [',num2str(p.bt(1)),', ',num2str(p.bt(2)),', ',num2str(p.bt(3)),']'],'Interpreter','latex','FontSize',16)
xlim([x(1),x(end)])
ylim([y(1),y(end)])
zlim([z(1),z(end)])
xlabel('x','Interpreter','latex','FontSize',16)
ylabel('y','Interpreter','latex','FontSize',16)
zlabel('z','Interpreter','latex','FontSize',16)
isosurface(X,Y,Z,P2,0.01, "HandleVisibility","off");
isosurface(X,Y,Z,P2,0.001, "HandleVisibility","off"); alpha(.5);
isosurface(X,Y,Z,P2,0.0001, "HandleVisibility","off"); alpha(.3);
colormap(cool);
drawnow;
%tightfig;

%{
fprintf('GGD \n')
vol = tripint(@GGD,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p),

p.comp = 1; m1 = tripint(@GGDM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 2; m2 = tripint(@GGDM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 3; m3 = tripint(@GGDM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
mean = [m1;m2;m3],

p.comp = 1; s11 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 2; s12 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 3; s13 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 4; s21 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 5; s22 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 6; s23 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 7; s31 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 8; s32 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 9; s33 = tripint(@GGDM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
cov = [s11 s12 s13; s21 s22 s23; s31 s32 s33],

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n \n')

fprintf('GGDAF \n')
vol = tripint(@GGDAF,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p),

p.comp = 1; m1 = tripint(@GGDAFM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 2; m2 = tripint(@GGDAFM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 3; m3 = tripint(@GGDAFM1,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
mean = [m1;m2;m3],

p.comp = 1; s11 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 2; s12 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 3; s13 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 4; s21 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 5; s22 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 6; s23 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 7; s31 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 8; s32 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
p.comp = 9; s33 = tripint(@GGDAFM2,-10*sigma1+mu1,10*sigma1+mu1,-10*sigma2+mu2,10*sigma2+mu2,-10*sigma3+mu3,10*sigma3+mu3,p);
cov = [s11 s12 s13; s21 s22 s23; s31 s32 s33],
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGD(X,Y,Z,p)
    d = length(p.mu); 
    
    Bi = B(p.b,d);
    Ai = A(p.b,d,p.Sigma);
    
    x = [X;Y;Z];

    P = Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDM1(X,Y,Z,p)
    d = length(p.mu); 
    
    Bi = B(p.b,d);
    Ai = A(p.b,d,p.Sigma);
    
    x = [X;Y;Z];

    P = (x(p.comp)-p.mu(p.comp))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDM2(X,Y,Z,p)
    d = length(p.mu); 

    Bi = B(p.b,d);
    Ai = A(p.b,d,p.Sigma);

    x = [X;Y;Z];

    switch p.comp
        case 1
            P = (X-p.mu(1))^2*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 2
            P = (X-p.mu(1))*(Y-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 3
            P = (X-p.mu(1))*(Z-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 4
            P = (X-p.mu(1))*(Y-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 5
            P = (Y-p.mu(2))^2*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 6
            P = (Y-p.mu(2))*(Z-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 7
            P = (X-p.mu(1))*(Z-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 8
            P = (Y-p.mu(2))*(Z-p.mu(2))*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        case 9
            P = (Z-p.mu(2))^2*Ai*exp(-(Bi*(x-p.mu)'*inv(p.Sigma)*(x-p.mu))^p.b);
        otherwise
            error("Incorrect value of component.")
    end
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAF(X,Y,Z,p)
    d = length(p.mu);
    if all(p.bt == p.bt(1))
        Ai = A(p.bt(1),d,p.Sigma);
    else
        Ai = (det(p.Sigma));
        for i = 1:length(p.mu)
            Ai = Ai*A(p.bt(i),1,p.Sigma);
        end
    end
    [S,D] = eig(p.Sigma); L = sqrtm(inv(D));
    
    x = [X; Y; Z];
    q = L*S'*(x-p.mu);
    sum1 = 0; 
    for i = 1:length(p.mu)
        sum1 = sum1 + (Bt(p.bt(i),d,p.bt)*q(i)^2)^p.bt(i); 
    end 
    sum2 = 0;
    for i = 1:3
        for j = i+1:3
            Bi = Bt(p.bt(i),d,p.bt);
            kd = p.bt(i)==p.bt(j);
            sum2_t = 0;
            for k = 1:p.bt(i)-1
                sum2_t = sum2_t + nchoosek(p.bt(i),k)*(q(i)^2)^(p.bt(i)-k)*(q(j)^2)^k;
            end
            sum2 = sum2 + Bi^(p.bt(i))*kd*sum2_t;
        end
    end
    sum3 = 0;
    Bi = B(p.bt(1),d);
    kd = (p.bt(1)==p.bt(2))&&(p.bt(2)==p.bt(3));
    if(kd~=0)
        for i = 1:p.bt(1)-2
            for j = 1:p.bt(1)-2
                for k = 1:p.bt(1)-2
                    if (i + j + k == p.bt(1))
                        n = factorial(p.bt(1))/(factorial(i)*factorial(j)*factorial(k));
                        sum3 = sum3 + n*(q(1)^2)^i*(q(2)^2)^j*(q(3)^2)^k;
                    end
                end
            end
        end
    end
    sum3 = Bi^(p.bt(1))*kd*sum3; 
    P = Ai*exp(-(sum1 + sum2 + sum3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAFM1(X,Y,Z,p)
    d = length(p.mu);
    if all(p.bt == p.bt(1))
        Ai = A(p.bt(1),d,p.Sigma);
    else
        Ai = (det(p.Sigma));
        for i = 1:length(p.mu)
            Ai = Ai*A(p.bt(i),1,p.Sigma);
        end
    end
    [S,D] = eig(p.Sigma); L = sqrtm(inv(D));
    
    x = [X; Y; Z];
    q = L*S'*(x-p.mu);
    sum1 = 0; 
    for i = 1:length(p.mu)
        sum1 = sum1 + (Bt(p.bt(i),d,p.bt)*q(i)^2)^p.bt(i); 
    end 
    sum2 = 0;
    for i = 1:3
        for j = i+1:3
            Bi = Bt(p.bt(i),d,p.bt);
            kd = p.bt(i)==p.bt(j);
            sum2_t = 0;
            for k = 1:p.bt(i)-1
                sum2_t = sum2_t + nchoosek(p.bt(i),k)*(q(i)^2)^(p.bt(i)-k)*(q(j)^2)^k;
            end
            sum2 = sum2 + Bi^(p.bt(i))*kd*sum2_t;
        end
    end
    sum3 = 0;
    Bi = B(p.bt(1),d);
    kd = (p.bt(1)==p.bt(2))&&(p.bt(2)==p.bt(3));
    if(kd~=0)
        for i = 1:p.bt(1)-2
            for j = 1:p.bt(1)-2
                for k = 1:p.bt(1)-2
                    if (i + j + k == p.bt(1))
                        n = factorial(p.bt(1))/(factorial(i)*factorial(j)*factorial(k));
                        sum3 = sum3 + n*(q(1)^2)^i*(q(2)^2)^j*(q(3)^2)^k;
                    end
                end
            end
        end
    end
    sum3 = Bi^(p.bt(1))*kd*sum3; 
    P = (x(p.comp)-p.mu(p.comp))*Ai*exp(-(sum1 + sum2 + sum3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = GGDAFM2(X,Y,Z,p)
    d = length(p.mu);
    if all(p.bt == p.bt(1))
        Ai = A(p.bt(1),d,p.Sigma);
    else
        Ai = (det(p.Sigma));
        for i = 1:length(p.mu)
            Ai = Ai*A(p.bt(i),1,p.Sigma);
        end
    end
    [S,D] = eig(p.Sigma); L = sqrtm(inv(D));
    
    x = [X; Y; Z];
    q = L*S'*(x-p.mu);
    sum1 = 0; 
    for i = 1:length(p.mu)
        sum1 = sum1 + (Bt(p.bt(i),d,p.bt)*q(i)^2)^p.bt(i); 
    end 
    sum2 = 0;
    for i = 1:3
        for j = i+1:3
            Bi = Bt(p.bt(i),d,p.bt);
            kd = p.bt(i)==p.bt(j);
            sum2_t = 0;
            for k = 1:p.bt(i)-1
                sum2_t = sum2_t + nchoosek(p.bt(i),k)*(q(i)^2)^(p.bt(i)-k)*(q(j)^2)^k;
            end
            sum2 = sum2 + Bi^(p.bt(i))*kd*sum2_t;
        end
    end
    sum3 = 0;
    Bi = B(p.bt(1),d);
    kd = (p.bt(1)==p.bt(2))&&(p.bt(2)==p.bt(3));
    if(kd~=0)
        for i = 1:p.bt(1)-2
            for j = 1:p.bt(1)-2
                for k = 1:p.bt(1)-2
                    if (i + j + k == p.bt(1))
                        n = factorial(p.bt(1))/(factorial(i)*factorial(j)*factorial(k));
                        sum3 = sum3 + n*(q(1)^2)^i*(q(2)^2)^j*(q(3)^2)^k;
                    end
                end
            end
        end
    end
    sum3 = Bi^(p.bt(1))*kd*sum3; 
    switch p.comp
        case 1
            P = (X-p.mu(1))^2*Ai*exp(-(sum1 + sum2 + sum3));
        case 2
            P = (X-p.mu(1))*(Y-p.mu(2))*Ai*exp(-(sum1 + sum2 + sum3));
        case 3
            P = (X-p.mu(1))*(Z-p.mu(3))*Ai*exp(-(sum1 + sum2 + sum3));
        case 4
            P = (X-p.mu(1))*(Y-p.mu(2))*Ai*exp(-(sum1 + sum2 + sum3));
        case 5
            P = (Y-p.mu(2))^2*Ai*exp(-(sum1 + sum2 + sum3));
        case 6
            P = (Y-p.mu(2))*(Z-p.mu(3))*Ai*exp(-(sum1 + sum2 + sum3));
        case 7
            P = (X-p.mu(1))*(Z-p.mu(3))*Ai*exp(-(sum1 + sum2 + sum3));
        case 8
            P = (Y-p.mu(2))*(Z-p.mu(3))*Ai*exp(-(sum1 + sum2 + sum3));
        case 9
            P = (Z-p.mu(3))^2*Ai*exp(-(sum1 + sum2 + sum3));
        otherwise
            error("Incorrect value of component.")
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = A(beta,d,Sigma)
    f = (B(beta,d)/pi)^(d/2)*(gamma(d/2)*beta)/(gamma(d/(2*beta))*det(Sigma)^(1/2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = Bt(beta_i,d,beta)
    if(sum(beta == beta_i)==1)
        f = B(beta_i,1);
    elseif(sum(beta == beta_i)==2)
        f = B(beta_i,2);
    else
        f = B(beta_i,d);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = B(beta,d)
    f = gamma((d+2)/(2*beta))/(d*gamma(d/(2*beta)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = tripint(fun, xmin, xmax, ymin, ymax, zmin, zmax,p)
    n = 100;
    h1 = (xmax-xmin)/n;
    h2 = (ymax-ymin)/n;
    h3 = (zmax-zmin)/n;
    
    f = 0;
    for x = xmin:h1:xmax
        for y = ymin:h2:ymax
            for z = zmin:h3:zmax
                f = f + fun(x,y,z,p);
            end
        end
    end
    f = f * (h1*h2*h3); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%