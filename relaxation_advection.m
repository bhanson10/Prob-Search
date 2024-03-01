% relaxation_advection.m 
%
% This code takes a PDF and calculates the relaxation advection using the
% following formula: D = lambda*I, v(x) = (lambda/p(x))*grad(p(x)). The code then
% propagates a uniform PDF using the calculated advection and
% diffusion. Example problems are 2D gaussians and beta-order super-gaussians,
% where the advection is known analytically. Other examples are numerically
% approximated, via convex/non-convex shapes of datasets and KDE. 
%
% By Benjamin L. Hanson, 2024

clear; close all; clc; 

IM_TYPE = 1; %0: Paper figure, 1: Presentation Figure

% Gaussian Implementation
dt = 0.005; N=51; L=6; d=L/(N-1); x=[-L/2:d:L/2]; y=[-L/2:d:L/2]; d_x = d; d_y = d;
xvbar=[0; 0]; Sigma =[1 0.4; 0.4 1]; lambda=1; [X_mesh,Y_mesh] = meshgrid(x,y); beta = 2; flag=0; % 0: analytical solution
if(IM_TYPE==1)
    title_str = "./Figures/Presentation/Gauss_" + num2str(beta)+ "/Gauss_" + num2str(beta); 
else
    title_str = "./Figures/Paper/Gauss_" + num2str(beta)+ "/Gauss_" + num2str(beta); 
end
[p,v_x,v_y]=gaussian(N,x,y,xvbar,Sigma,lambda,beta,d_x,d_y,flag,title_str,IM_TYPE);

% Kidney Bean Implementation
%{
dt = 0.001; lambda=1; N = 51; Lmin=-1; Lmax = 2.5; d=(Lmax-Lmin)/(N-1); d_x = d; d_y = d;
x=[Lmin:d:Lmax]; y=[Lmin:d:Lmax]; [X_mesh,Y_mesh] = meshgrid(x,y); flag = 1; % numerical 
[p,v_x,v_y]=kidney_bean(N,x,y,d,lambda,IM_TYPE); 
if(IM_TYPE==1)
    title_str = "./Figures/Presentation/Kidney_Bean"; 
else
    title_str = "./Figures/Paper/Kidney_Bean"; 
end
%}

% Numerical Implementation
%{
flag = 3; % 1:alpha, 2:KSD Den, 3:KSD Lake
if(IM_TYPE==1)
    if(flag==1)
        load ./Datasets/alpha.mat
        title_str = "./Figures/Presentation/Alpha/alpha"; 
    elseif(flag==2)
        load ./Datasets/kde_den.mat
        title_str = "./Figures/Presentation/Den/kde_den"; 
    elseif(flag==3)
        load ./Datasets/kde_lake.mat
        title_str = "./Figures/Presentation/Lake/kde_lake"; 
    end
else
    if(flag==1)
        load ./Datasets/alpha.mat
        title_str = "./Figures/Paper/Alpha/alpha"; 
    elseif(flag==2)
        load ./Datasets/kde_den.mat
        title_str = "./Figures/Paper/Den/kde_den"; 
    elseif(flag==3)
        load ./Datasets/kde_lake.mat
        title_str = "./Figures/Paper/Lake/kde_lake"; 
    end
end
dt = 0.001; lambda=1; N = 51; Lmin_x = -1.5; Lmax_x = 1.5; d_x=(Lmax_x-Lmin_x)/(N-1);
Lmin_y = -2; Lmax_y = 2; d_y=(Lmax_y-Lmin_y)/(N-1); 
x=[Lmin_x:d_x:Lmax_x]; y=[Lmin_y:d_y:Lmax_y]; [X_mesh,Y_mesh] = meshgrid(x,y); 
[p,v_x,v_y]=numerical(N,x,y,d_x,d_y,lambda,shp,X_mesh,Y_mesh,flag,title_str,IM_TYPE);
%}

% Propagate a new, different PDF to the relaxation goal PDF using the calculated advection field 
PDF_U = zeros(N,N); 
max_val = max(p,[],'all'); %Used for plots

% Setting PDF - Uniform PDF with area under the curve of 1
for i=2:N-1
    for j=2:N-1
        PDF_U(j,i)= 1; 
    end
end
PDF_U = PDF_U/(d_x*d_y*sum(PDF_U,'all')); 

clear F_iso; clear F_top; 
figure(5); clf;  hold on; view(30,30); %Iso View
set(gca, 'FontName' , 'Times','FontSize',16);
title('Initial Uniform PDF');
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times');
colorbar;
clim([0 max_val]);
xlim([x(1), x(end)])
ylim([y(1), y(end)])
surf(X_mesh,Y_mesh,reshape(PDF_U,[N,N]), 'EdgeColor','none'); 
F_iso(1) = getframe(gcf);
set( gcf , 'Color' , 'w' );
drawnow

figure(6); clf;  hold on; view(0,90); %Top View
set(gca, 'FontName' , 'Times','FontSize',16);
title('Initial Uniform PDF');
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
colorbar;
clim([0 max_val]);
if(flag==0)
    axis equal; 
end
xlim([x(1), x(end)])
ylim([y(1), y(end)])
if(IM_TYPE==1)
    contour(X_mesh,Y_mesh,reshape(PDF_U,[N,N]),[linspace(0.01,max_val,10)], 'LineWidth',2, 'Fill', 'on');
else
    contour(X_mesh,Y_mesh,reshape(PDF_U,[N,N]),[linspace(0.01,max_val,10)], 'LineWidth',2);
end
F_top(1) = getframe(gcf); 
set( gcf , 'Color' , 'w' );
drawnow

%Initializing F component 
F =  f(d_x,d_y,v_x,v_y,N);

B = cell(1,N); A = cell(1,N-1); C = A; % Cell array of A,B,C matrices

% Initializing block B matrices
for i=1:N
    B_i = zeros(N,N);
    B2 = zeros(1,N); 
    for j=1:N
        B2(1,j) = -F(i,j)-(((2*lambda)/(d_x^2))+((2*lambda)/(d_y^2)));
    end
    B_i = B_i + diag(B2,0);
    B{1,i} = B_i; 
end

for i=1:N
    B_i = B{1,i};
    B3 = zeros(1,N-1); B1 = B3;  
    for j=1:(N-1)
        B1(1,j) = (v_y(j+1,i)/(2*d_y))+(lambda/d_y^2);
        B3(1,j) = (-v_y(j,i)/(2*d_y))+(lambda/d_y^2);
    end
    B_i = B_i + diag(B3,1)+diag(B1,-1);
    B{1,i} = B_i; 
end

for i=1:N-1
    A_i = zeros(N,N); C_i = A_i;
    A1 = zeros(1,N); C1 = A1; 
    for j=1:(N)
        A1(1,j) = (v_x(j,i+1)/(2*d_x))+(lambda/(d_x^2));
        C1(1,j) = (-v_x(j,i)/(2*d_x))+(lambda/(d_x^2));
    end
    A_i = A_i + diag(A1,0); C_i = C_i + diag(C1,0);
    A{1,i} = A_i; C{1,i} = C_i; 
end

M = blkdiag(B{:}); % Main Block Diagonal
M_A = blkdiag(A{:}); M_A = [[zeros(N,N*N-N);M_A] zeros(N*N,N)]; % Lower Block Diagonal
M_C = blkdiag(C{:}); M_C = [zeros(N*N,N) [M_C;zeros(N,N*N-N)]]; % Upper Block Diagonal
M = sparse(M + M_A + M_C); % Complete M Matrix

PDF = reshape(PDF_U, [N*N,1]); p = reshape(p,[N*N,1]); 
eps=sum(abs(p-PDF)); diff = 1; % Initializing Difference to be above 0.01

timestep = 1; timestep_list = [timestep]; eps_list = [eps]; 
while(diff > 0.01)
  K1=M*PDF;
  K2=M*(PDF+(dt/2)*K1);
  K3=M*(PDF+(dt/2)*K2);
  K4=M*(PDF+dt*K3);
  PDF = PDF+dt*((K1/6)+(K2/3)+(K3/3)+(K4/6));
  PDF = boundary_conditions(PDF,N,d_x,d_y);
  [eps_new, F1, F2] = plot_PDF(p,PDF,timestep,dt,x,y,X_mesh,Y_mesh,N,max_val,timestep_list,eps_list,flag,title_str,IM_TYPE); 
  F_iso(timestep) = F1; F_top(timestep) = F2;
  timestep=timestep+1; diff = abs(eps-eps_new); eps = eps_new; 
  timestep_list(end+1) = timestep; eps_list(end+1) = eps_new;
end 
timestep=timestep-1; 

figure(5); clf;  hold on; view(30,30); %Iso View
set(gca, 'FontName' , 'Times','FontSize',16);
title([]);
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
colorbar;
clim([0 max_val]);
xlim([x(1), x(end)])
ylim([y(1), y(end)])
surf(X_mesh,Y_mesh,reshape(PDF,[N,N]), 'EdgeColor','none'); 
%exportgraphics(gca,title_str+'_frame_final_iso.eps','Resolution',300)
set( gcf , 'Color' , 'w' );
drawnow

figure(6); clf;  hold on; view(0,90); %Top View
set(gca, 'FontName' , 'Times','FontSize',16);
title([]);
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
colorbar;
clim([0 max_val]);
if(flag==0)
    axis equal;
end
xlim([x(1), x(end)])
ylim([y(1), y(end)])
if(IM_TYPE==1)
    contour(X_mesh,Y_mesh,reshape(PDF,[N,N]),[linspace(0.01,max_val,10)], 'LineWidth',2, 'Fill', 'on');
else
    contour(X_mesh,Y_mesh,reshape(PDF,[N,N]),[linspace(0.01,max_val,10)], 'LineWidth',2);
end

%exportgraphics(gca,title_str+'_frame_final_top.eps','Resolution',300)
set( gcf , 'Color' , 'w' );
drawnow

figure(7); clf;  hold on;
set(gca, 'FontName' , 'Times','FontSize',16);
title([char(949), '(t)']);
xlabel('timestep', 'FontSize', 16, 'FontName', 'Times')
ylabel('\epsilon', 'FontSize', 16, 'FontName', 'Times')
xlim([1,inf])
ylim([0,inf])
plot(timestep_list, eps_list);
%exportgraphics(gca,title_str+'_epsilon.eps','Resolution',300)
set( gcf , 'Color' , 'w' );
drawnow
    
%create_video(F_top, title_str + '_top_2D.mp4');
%create_video(F_iso, title_str + '_iso_2D.mp4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,v_x,v_y]=gaussian(N,x,y,xvbar,Sigma,lambda,beta,d_x,d_y,flag,title_str,IM_TYPE)
v_x = zeros(N,N); v_y =v_x; dim=2;

B = gamma((dim+2)/(2*beta))/(dim*gamma(dim/(2*beta)));
A = (B/pi)^(dim/2)*(gamma(dim/2)*beta)/(gamma(dim/(2*beta))*det(Sigma)^(-1/2));

for i=1:N
    for j=1:N
        xv=[x(i); y(j)];
        p(j,i) = A*exp(-(B*(xv-xvbar)'*inv(Sigma)*(xv-xvbar))^beta);
        v = (-2*lambda*beta)*B*(B*(xv-xvbar)'*inv(Sigma)*(xv-xvbar))^(beta-1)*inv(Sigma)*(xv-xvbar);  
        v_x(j,i) = v(1,1); v_y(j,i) = v(2,1); 
    end
end
p = reshape(p,[N,N]);

make_plots(x,y,p,v_x,v_y,title_str,flag,lambda, IM_TYPE);
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,v_x,v_y]=kidney_bean(N,x,y,d,lambda, IM_TYPE)
p = zeros(N,N); dt = 0.0005; 

for i=1:N
    for j=1:N
        if((x(i)^2+y(j)^2)^2 <= 2*(x(i)^3 + y(j)^3)) % 2D Bean Curve
            p(i,j) = 1;
        end
    end
end
p = p/(d*d*sum(p,'all'));

% Initializing Block Matrices for diffusion-only time updating
B_1 = cell(1,N); A_1 = cell(1,N-1); % Cell array of A,B,C matrices

B2 = -2*((lambda/(d^2))+(lambda/(d^2)))*ones(1,N);
B1 = (lambda/d^2)*ones(1,N-1);
B_i = zeros(N,N); B_i = B_i + diag(B2,0)+ diag(B1,-1) + diag(B1,1);

for i=1:N
    B_1{1,i} = B_i; 
end

A1 = (lambda/d^2)*ones(1,N); 
A_i = zeros(N,N); A_i = A_i + diag(A1,0); 
for i=1:N-1
    A_1{1,i} = A_i; 
end
C_1 = A_1; 

M_1 = blkdiag(B_1{:}); % Main Block Diagonal
M_A = blkdiag(A_1{:}); M_A = [[zeros(N,N*N-N);M_A] zeros(N*N,N)]; % Lower Block Diagonal
M_C = blkdiag(C_1{:}); M_C = [zeros(N*N,N) [M_C;zeros(N,N*N-N)]]; % Upper Block Diagonal
M_1 = sparse(M_1 + M_A + M_C); % Complete M Matrix
p = reshape(p,[N*N,1]); 

%RK4 Propagating the PDF
for i=1:30
  K1=M_1*p;
  K2=M_1*(p+(dt/2)*K1);
  K3=M_1*(p+(dt/2)*K2);
  K4=M_1*(p+dt*K3);
  p = p+dt*((K1/6)+(K2/3)+(K3/3)+(K4/6));
end 

p = boundary_conditions(p,N,d_x,d_y); p = reshape(p,[N,N]); 
[v_x,v_y] = gradient(p,d_x,d_y); 
v_x = (lambda).*(v_x./p); v_y = (lambda).*(v_y./p);
v_x(:,1) = 0; v_x(:,end) = 0; v_x(1,:) = 0; v_x(end,:) = 0; 
v_y(:,1) = 0; v_y(:,end) = 0; v_y(1,:) = 0; v_y(end,:) = 0; 

make_plots(x,y,p,v_x,v_y, "Kidney_Bean",2,lambda, IM_TYPE);
end % end function kidney_bean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,v_x,v_y]=numerical(N,x,y,d_x,d_y,lambda,shp,X_mesh,Y_mesh,flag,title_str, IM_TYPE)
p = zeros(N,N); dt = 0.001;

if(flag==1) % Alpha-convex hull approximation
    for i=1:N
        xq = x(i);
        for j=1:N
            yq = y(j); 
            if (inShape(shp,xq,yq) == 1)
                p(j,i)=1;
            end
        end
    end
    p = p/(d_x*d_y*sum(p,'all'));
    
    % Initializing Block Matrices for diffusion-only time updating
    B_1 = cell(1,N); A_1 = cell(1,N-1); % Cell array of A,B,C matrices
    
    B2 = -2*((lambda/(d_x^2))+(lambda/(d_y^2)))*ones(1,N);
    B1 = (lambda/d_y^2)*ones(1,N-1);
    B_i = zeros(N,N); B_i = B_i + diag(B2,0)+ diag(B1,-1) + diag(B1,1);
    
    for i=1:N
        B_1{1,i} = B_i; 
    end
    
    A1 = (lambda/d_x^2)*ones(1,N); 
    A_i = zeros(N,N); A_i = A_i + diag(A1,0); 
    for i=1:N-1
        A_1{1,i} = A_i; 
    end
    C_1 = A_1; 
    
    M_1 = blkdiag(B_1{:}); % Main Block Diagonal
    M_A = blkdiag(A_1{:}); M_A = [[zeros(N,N*N-N);M_A] zeros(N*N,N)]; % Lower Block Diagonal
    M_C = blkdiag(C_1{:}); M_C = [zeros(N*N,N) [M_C;zeros(N,N*N-N)]]; % Upper Block Diagonal
    M_1 = sparse(M_1 + M_A + M_C); % Complete M Matrix
    p = reshape(p,[N*N,1]); 
    
    %RK4 Propagating the PDF
    for i=1:20
      K1=M_1*p;
      K2=M_1*(p+(dt/2)*K1);
      K3=M_1*(p+(dt/2)*K2);
      K4=M_1*(p+dt*K3);
      p = p+dt*((K1/6)+(K2/3)+(K3/3)+(K4/6));
    end 
    p = p/(d_x*d_y*sum(p));
    
    p = reshape(p,[N,N]); 
    [v_x,v_y] = gradient(p,d_x,d_y); 
    v_x = (lambda).*(v_x./p); v_y = (lambda).*(v_y./p);

elseif(flag==2||flag==3) % KSD Approximation
    pts = [reshape(X_mesh, [numel(X_mesh),1]) reshape(Y_mesh, [numel(Y_mesh),1])];
    [p,~] = ksdensity(shp.Points,pts,'Bandwidth',0.1); p = boundary_conditions(p,N,d_x,d_y); 
    p = reshape(p, [N,N]); [v_x,v_y] = gradient(p,d_x,d_y); 
    v_x = (lambda).*(v_x./p); v_y = (lambda).*(v_y./p);
    
    % Initializing Block Matrices for diffusion-only time updating
    B_1 = cell(1,N); A_1 = cell(1,N-1); % Cell array of A,B,C matrices
    
    B2 = -2*((lambda/(d_x^2))+(lambda/(d_y^2)))*ones(1,N);
    B1 = (lambda/d_y^2)*ones(1,N-1);
    B_i = zeros(N,N); B_i = B_i + diag(B2,0)+ diag(B1,-1) + diag(B1,1);
    
    for i=1:N
        B_1{1,i} = B_i; 
    end
    
    A1 = (lambda/d_x^2)*ones(1,N); 
    A_i = zeros(N,N); A_i = A_i + diag(A1,0); 
    for i=1:N-1
        A_1{1,i} = A_i; 
    end
    C_1 = A_1; 
    
    M_1 = blkdiag(B_1{:}); % Main Block Diagonal
    M_A = blkdiag(A_1{:}); M_A = [[zeros(N,N*N-N);M_A] zeros(N*N,N)]; % Lower Block Diagonal
    M_C = blkdiag(C_1{:}); M_C = [zeros(N*N,N) [M_C;zeros(N,N*N-N)]]; % Upper Block Diagonal
    M_1 = sparse(M_1 + M_A + M_C); % Complete M Matrix
    p = reshape(p,[N*N,1]); 
    
    %RK4 Propagating the PDF
    for i=1:5
      K1=M_1*p;
      K2=M_1*(p+(dt/2)*K1);
      K3=M_1*(p+(dt/2)*K2);
      K4=M_1*(p+dt*K3);
      p = p+dt*((K1/6)+(K2/3)+(K3/3)+(K4/6));
    end 
    p = p/(d_x*d_y*sum(p));
    p = reshape(p,[N,N]); 
    [v_x,v_y] = gradient(p,d_x,d_y); 
    v_x = (lambda).*(v_x./p); v_y = (lambda).*(v_y./p);
end

make_plots(x,y,p,v_x,v_y,title_str,2,lambda, IM_TYPE);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_plots(x,y,p,v_x,v_y,title_str,flag,lambda,IM_TYPE)
    [X_mesh,Y_mesh] = meshgrid(x,y); max_val = max(p, [], 'all');

    figure(1); clf;  hold on; view(0,90); %Top View
    set(gca, 'FontName' , 'Times','FontSize',12);

    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
    colorbar;
    clim([0 max_val]);
    if(flag==0)
        axis equal; 
    end
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    if(IM_TYPE==1)
        contour(X_mesh,Y_mesh,p,[linspace(0.01,max_val,10)], 'LineWidth',2, 'Fill', 'on');
    else
        contour(X_mesh,Y_mesh,p,[linspace(0.01,max_val,10)], 'LineWidth',2);
    end
    ax = gca;
    %exportgraphics(ax,title_str + '_relax_pdf_top.eps','Resolution',300)
    drawnow

    figure(2); clf;  hold on; view(30,30); %Top View
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
    colorbar;
    clim([0 max_val]);
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    surf(X_mesh,Y_mesh,p,'EdgeColor','none');
    ax = gca;
    %exportgraphics(ax,title_str + '_relax_pdf_iso.eps','Resolution',300)
    drawnow

    figure(3); clf;  hold on;
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    if(flag==0)
        axis equal; 
    end
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    quiver(X_mesh,Y_mesh,v_x,v_y)
    ax = gca;
    %exportgraphics(ax,title_str+'_adv.eps','Resolution',300)
    drawnow

    phi = lambda*log(p); 
    
    figure(4); clf;  hold on; view(0,90); %Top View
    set(gca, 'FontName' , 'Times','FontSize',12);
    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
    if(flag==0)
        axis equal; 
    end
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])
    cmin = min(phi,[],'all');
    cmax = max(phi,[],'all'); 
    levels = -flip(logspace(log10(abs(cmax)),log10(abs(cmin)),20),2); 
    cmp = colormap; 
    cmp = flipud(cmp);
    colormap(cmp); 
    colorbar; 
    if(IM_TYPE==1)
        contour(X_mesh,Y_mesh,phi,levels, 'LineWidth',2, 'Fill', 'on');
    else
        contour(X_mesh,Y_mesh,phi,levels, 'LineWidth',2);
    end
    ax = gca;
    %exportgraphics(ax,title_str + '_relax_phi_top.eps','Resolution',300)
    drawnow
end % end function make_plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDF = boundary_conditions(PDF,N,d_x,d_y)
    % Dirichlet Boundary Conditions, PDF = 0 at edges
    count = 1;
    for i=1:N
        for j=1:N
            if((i==1)||(i==N)||(j==1)||(j==N))
                PDF(count,1)=0;
            end
            count = count + 1; 
        end
    end
    PDF = PDF/(d_x*d_y*sum(PDF,'all'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eps_new,F1,F2] = plot_PDF(p,PDF,timestep,dt,x,y,X_mesh,Y_mesh,N,max_val,timestep_list,eps_list,flag,title_str, IM_TYPE)
 
    eps_new = sum(abs(p - PDF)); 
    PDF_plot = reshape(PDF,[N,N]);
    
    figure(5); clf;  hold on; view(30,30); %Iso View
    set(gca, 'FontName' , 'Times','FontSize',12);
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*dt), ', \Delta t = ', num2str(dt), ', ', char(949), ' = ', num2str(eps_new)]);
    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
    colorbar;
    clim([0 max_val]);
    xlim([x(1), x(end)])
    ylim([y(1), y(end)])
    surf(X_mesh,Y_mesh,PDF_plot, 'EdgeColor','none'); 
    F1 = getframe(gcf);
    if (mod(timestep,50)==0)
        title([]);
        %exportgraphics(gca,title_str+'_frame_'+num2str(timestep)+'_iso.eps','Resolution',300)
    end
    set( gcf , 'Color' , 'w' );
    drawnow

    figure(6); clf;  hold on; view(0,90); %Top View
    set(gca, 'FontName' , 'Times','FontSize',12);
    title(['iter = ', num2str(timestep), ', t = ', num2str(timestep*dt), ', \Delta t = ', num2str(dt), ', ', char(949), ' = ', num2str(eps_new)]);
    xlabel('x', 'FontSize', 16, 'FontName', 'Times')
    ylabel('y', 'FontSize', 16, 'FontName', 'Times')
    zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
    colorbar;
    clim([0 max_val]);
    if(flag==0)
        axis equal;
    end
    xlim([x(1), x(end)])
    ylim([y(1), y(end)])
    if(IM_TYPE==1)
        contour(X_mesh,Y_mesh,PDF_plot,[linspace(0.01,max_val,10)], 'LineWidth',2, 'Fill', 'on');
    else
        contour(X_mesh,Y_mesh,PDF_plot,[linspace(0.01,max_val,10)], 'LineWidth',2);
    end
    F2 = getframe(gcf);
    if (mod(timestep,50)==0)
        title([]);
        %exportgraphics(gca,title_str+'_frame_'+num2str(timestep)+'_top.eps','Resolution',300)
    end
    set( gcf , 'Color' , 'w' );
    drawnow

    figure(7); clf;  hold on;
    set(gca, 'FontName' , 'Times','FontSize',12);
    title([char(949), '(t)']);
    xlabel('timestep', 'FontSize', 16, 'FontName', 'Times')
    ylabel('\epsilon', 'FontSize', 16, 'FontName', 'Times')
    xlim([1,inf])
    ylim([0,inf])
    plot(timestep_list, eps_list);
    set( gcf , 'Color' , 'w' );
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(F,title)
    
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    for i=1:length(F)
        frame = F(i);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F] = f(d_x,d_y,v_x,v_y,N)
    F = zeros(N,N); 
    for i=1:N
        for j=1:N
            if((i~=1)&&(i~=N)&&(j~=1)&&(j~=N))
                F(i,j) = (v_x(j,i+1)-v_x(j,i-1))/(2*d_x)+(v_y(j+1,i)-v_y(j-1,i))/(2*d_y);
            elseif((i~=1)&&(i~=N)&&(j==1))
                F(i,j) = (v_x(j,i+1)-v_x(j,i-1))/(2*d_x)+(v_y(j+1,i)-v_y(j,i))/(d_y);
            elseif((i==1)&&(j~=1)&&(j~=N))
                F(i,j) = (v_x(j,i+1)-v_x(j,i))/(d_x)+(v_y(j+1,i)-v_y(j-1,i))/(2*d_y);
            elseif((i~=1)&&(i~=N)&&(j==N))
                F(i,j) = (v_x(j,i+1)-v_x(j,i-1))/(2*d_x)+(v_y(j,i)-v_y(j-1,i))/(d_y);
            elseif((i==N)&&(j~=1)&&(j~=N))
                F(i,j) = (v_x(j,i)-v_x(j,i-1))/(d_x)+(v_y(j+1,i)-v_y(j-1,i))/(2*d_y);
            elseif((i==1)&&(j==1))
                F(i,j) = (v_x(j,i+1)-v_x(j,i))/(d_x)+(v_y(j+1,i)-v_y(j,i))/(d_y);
            elseif((i==N)&&(j==N))
                F(i,j) = (v_x(j,i)-v_x(j,i-1))/(d_x)+(v_y(j,i)-v_y(j-1,i))/(d_y);
            elseif((i==1)&&(j==N))
                F(i,j) = (v_x(j,i+1)-v_x(j,i))/(d_x)+(v_y(j,i)-v_y(j-1,i))/(d_y);
            elseif((i==N)&&(j==1))
                F(i,j) = (v_x(j,i)-v_x(j,i-1))/(d_x)+(v_y(j+1,i)-v_y(j,i))/(d_y);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%