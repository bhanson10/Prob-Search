% evasive_search.m 
%
% This simulation takes a steady state PDF, generates the relaxation
% advection (either analytically or numerically depending on the
% steady-state), and includes the advection in a probablistic search
% simulation. The behavior of the target in this simulation is assumed to
% be evasive.
%
% By Benjamin L. Hanson, Feb 23 2024

clear all; close; clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%
sim.T = 4; sim.dt = 0.0025; sim.timesteps = round(sim.T/sim.dt); sim.FrameRate = round(.1/sim.dt);
theta = linspace(0,2*pi, 1000); 
N = 51; 
L_x = 7; L_y = 7;  
d_x = L_x/(N-1); d_y = L_y/(N-1);
x = linspace(-L_x/2,L_x/2,N); y = linspace(-L_y/2,L_y/2,N);
[X_mesh, Y_mesh] = meshgrid(x,y);
colors = ['g' 'b' 'r' 'k' 'm' 'c' 'y'];
%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target.lambda = 2; target.dif = target.lambda*ones(N+2); target.psi = 1e-1; 
target.stats_flag = 1; %1: Gaussian Statistics, 2: Numerical Territory
if(target.stats_flag == 1)
    target.s = 1; 
    target.xvbar=[0; 0]; target.P = [3 0; 0 3];
else
    flag = 1; % 1:alpha, 2:KSD Den, 3:KSD Lake
    if(flag==1)
        L_x = 3; L_y = 4;  
        d_x = L_x/(N-1); d_y = L_y/(N-1);
        x = linspace(-L_x/2,L_x/2,N); y = linspace(-L_y/2,L_y/2,N);
        load ./Datasets/alpha.mat
    elseif(flag==2)
        load ./Datasets/kde_den.mat
    elseif(flag==3)
        load ./Datasets/kde_lake.mat
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drone Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drones.num = 3; 
drones.ang_speed = [5 5 5]; 
drones.init_theta = [0 2*pi/3 4*pi/3];
drones.a = 1/sim.dt; drones.d = 8;

drones.orbit_flag = 4; %1: Concentric orbits, 2: Cassini ovals, 3: Lemniscates, 4: Rotating Cassini
if (drones.orbit_flag == 1)
    drones.radius = [0.5 1.5 2.5]; 
    title_str = 'circle'; 
elseif (drones.orbit_flag == 2)
    drones.focal = 1.5;
    drones.b = [1.51 1.7 1.9];
    drones.f = 1.41; 
elseif (drones.orbit_flag == 3)
    drones.phi = pi.*([1:drones.num]-1)./drones.num; 
    drones.c = cos(drones.phi); 
    drones.s = sin(drones.phi); 
    drones.f = 1/(2*drones.num);
    drones.Xm = cos(theta); drones.Ym = drones.f*sin(2.*theta); 
    drones.scale = 1.4; 
elseif (drones.orbit_flag == 4)
    theta = linspace(0,6*pi,1000); 
    drones.num = 4; 
    drones.ang_speed = [1 1 1 1]; 
    drones.init_theta = [0 pi/2 pi (3*pi)/2];
    drones.A = [2 2 2 2];
    drones.B = drones.A + [0.02 0.021 0.022 0.023]; 
    drones.v = 0.4; 
    title_str = 'rotating'; 
end

if(target.stats_flag == 1)
    drones.sigma = 0.1;  drones.P = [.25 0; 0 .25]; 
else
    drones.sigma = 0.1; drones.P = [1 0; 0 1]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drone Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Target Density/Relaxation Advection
if(target.stats_flag == 1)
    [target.p,target.v_x,target.v_y]=target_gauss(N,x,y,target.xvbar,target.P,target.s,target.lambda,d_x,d_y);
else
    [target.p,target.v_x,target.v_y]=target_numer(N,x,y,d_x,d_y,target.lambda,shp,X_mesh,Y_mesh,flag);
end

fig = figure(1); clf; hold on; fig.Position = [150 150 1200 600];
set(gca, 'FontName' , 'Times','FontSize',12);
sgtitle(['Pre-search, iter = 0, t = 0, \Delta', 't = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda)], 'FontSize', 16, 'FontName', 'Times');

subplot(1,2,1); hold on; 
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
if(target.stats_flag==1)
    axis square
end
view(0,90); colorbar;

subplot(1,2,2);  hold on; 
xlabel('x', 'FontSize', 16, 'FontName', 'Times')
ylabel('y', 'FontSize', 16, 'FontName', 'Times')
zlabel('Probability', 'FontSize', 16, 'FontName', 'Times')
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
zlim([0, inf]);
if(target.stats_flag==1)
    axis square
end
view(45,30);

subplot(1,2,1); 
plots_1(1) = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7);
subplot(1,2,2); 
plots_2(1) = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7);
for i=1:drones.num
    if (drones.orbit_flag == 1)
        x_orbit = drones.radius(i)*cos(theta);
        y_orbit = drones.radius(i)*sin(theta);
        drones.pos(i,1) = drones.radius(i)*cos(drones.init_theta(i));
        drones.pos(i,2) = drones.radius(i)*sin(drones.init_theta(i));
    elseif(drones.orbit_flag == 2)
        x_orbit = sqrt(drones.focal^2*cos(2.*theta)+sqrt(drones.focal^4*(cos(2.*theta).^2+drones.b(i)^4-drones.focal^4))).*cos(theta); 
        y_orbit = drones.f.*sqrt(drones.focal^2*cos(2.*theta)+sqrt(drones.focal^4*(cos(2.*theta).^2+drones.b(i)^4-drones.focal^4))).*sin(theta);
        drones.pos(i,1) = sqrt(drones.focal^2*cos(2.*drones.init_theta(i))+sqrt(drones.focal^4*(cos(2.*drones.init_theta(i)).^2+drones.b(i)^4-drones.focal^4))).*cos(drones.init_theta(i)); 
        drones.pos(i,2) = drones.f.*sqrt(drones.focal^2*cos(2.*drones.init_theta(i))+sqrt(drones.focal^4*(cos(2.*drones.init_theta(i)).^2+drones.b(i)^4-drones.focal^4))).*sin(drones.init_theta(i));
    elseif(drones.orbit_flag == 3)
        x_orbit = drones.scale.*(drones.c(i).*drones.Xm - drones.s(i).*drones.Ym); 
        y_orbit = drones.scale.*(drones.s(i).*drones.Xm + drones.c(i).*drones.Ym);
        Xm = cos(drones.init_theta(i)); Ym = drones.f*sin(2*drones.init_theta(i));
        drones.pos(i,1) = drones.scale.*(drones.c(i)*Xm - drones.s(i)*Ym); 
        drones.pos(i,2) = drones.scale.*(drones.s(i)*Xm + drones.c(i)*Ym); 
    elseif(drones.orbit_flag == 4)
         b=(-2*drones.A(i)^2).*cos(4.*(theta + drones.init_theta(i))); c=drones.A(i)^4-drones.B(i)^4; r=sqrt((-b+sqrt(b.^2-4*c))/2);
         x_orbit = r.*sin((theta + drones.init_theta(i))+(theta).*drones.v);
         y_orbit = r.*cos((theta + drones.init_theta(i))+(theta).*drones.v);
         drones.pos(i,1) = x_orbit(1);
         drones.pos(i,2) = y_orbit(1);
    end

    subplot(1,2,1); 
    plot(x_orbit,y_orbit,colors(i),'Linewidth',1,'HandleVisibility','off'); 
    plots_1(i+1) = scatter(drones.pos(i,1),drones.pos(i,2),75,colors(i),'filled', 'HandleVisibility', 'off');
    subplot(1,2,2); 
    plot(x_orbit,y_orbit,colors(i),'Linewidth',1,'HandleVisibility','off'); 
    plots_2(i+1) = scatter(drones.pos(i,1),drones.pos(i,2),75,colors(i),'filled', 'HandleVisibility', 'off');
end
frames(1) = getframe(gcf); 
drawnow 

pause(3); delete(plots_1); delete(plots_2);

target.relax_M = M(d_x,d_y,N,target); target.relax_v_x = target.v_x; target.relax_v_y = target.v_y; 

t = 0; 
for i=1:sim.timesteps
    drones = update_pos(drones,t);
    [target.v_x, target.v_y] = update_vel(drones,x,y,N,target);
    target.dif = target.lambda + target.psi.*update_dif(drones,x,y,N);

    p_drones_sum = zeros([N*N,1]);
    for j=1:drones.num
        p_drones_sum = p_drones_sum + search_gauss(N,x,y,[drones.pos(j,1);drones.pos(j,2)],drones.sigma,drones.a);
        subplot(1,2,1); 
        plots_1(j+1) = scatter(drones.pos(j,1),drones.pos(j,2),75,colors(j),'filled', 'HandleVisibility', 'off');
        subplot(1,2,2); 
        plots_2(j+1) = scatter(drones.pos(j,1),drones.pos(j,2),75,colors(j),'filled', 'HandleVisibility', 'off');
    end
    
    target.M = M(d_x,d_y,N,target) + drone_M(N,target,p_drones_sum);

    K1=target.M*target.p;
    K2=target.M*(target.p+(sim.dt/2)*K1);
    K3=target.M*(target.p+(sim.dt/2)*K2);
    K4=target.M*(target.p+sim.dt*K3);
    target.p = target.p+sim.dt*((K1/6)+(K2/3)+(K3/3)+(K4/6)); 
    target.p = boundary_conditions(target.p,N,d_x,d_y);

    subplot(1,2,1); 
    plots_1(1) = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7); 
    subplot(1,2,2); 
    plots_2(1) = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7); 
    
    sgtitle(['Evasive searching, iter = ',num2str(i), ', t = ', num2str(t), ', \Delta t = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda), ', \psi = ', num2str(target.psi), ', a = ', num2str(drones.a*sim.dt), '/\Delta','t, \sigma = ', num2str(drones.sigma),', d = ', num2str(drones.d)], 'FontSize', 16, 'FontName', 'Times');
    frames(i+1) = getframe(gcf);  
    %{
    if(((i==50)||(i==120))||(i==250))
        sgtitle([]);
        exportgraphics(gcf,"./Figures/Search/Evasive/evasive_search_rotating_" + num2str(i) + "_tight.eps",'Resolution',300)
    end
    %}
    drawnow;

    if(i ~= sim.timesteps)
        delete(plots_1); delete(plots_2);  
    end
    t = t + sim.dt; 
end
pause(3); delete(plots_1); delete(plots_2);

for i=sim.timesteps+1:sim.timesteps+200
    K1=target.relax_M*target.p;
    K2=target.relax_M*(target.p+(sim.dt/2)*K1);
    K3=target.relax_M*(target.p+(sim.dt/2)*K2);
    K4=target.relax_M*(target.p+sim.dt*K3);
    target.p = target.p+sim.dt*((K1/6)+(K2/3)+(K3/3)+(K4/6));
    target.p = boundary_conditions(target.p,N,d_x,d_y);
    
    subplot(1,2,1);
    plot4_1 = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7); 
    subplot(1,2,2);
    plot4_2 = surf(X_mesh,Y_mesh,reshape(target.p,[N,N]),'EdgeColor','none','FaceAlpha',0.7); 
    sgtitle(['Post-search, iter = ',num2str(i), ', t = ', num2str(t), ', \Delta', 't = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda)], 'FontSize', 16, 'FontName', 'Times');
    frames(i+1) = getframe(gcf);
    drawnow; 
    if(i ~= sim.timesteps+200)
        delete(plot4_1); delete(plot4_2);
    end
    t = t + sim.dt; 
end

%create_video(frames, sim, ['./Figures/Search/Evasive/evasive_search_',title_str,'.mp4']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,v_x,v_y]=target_gauss(N,x,y,xvbar,P,s,lambda,d_x,d_y)
    v_x = zeros(N,N); v_y =v_x; dim=2;
    
    B = gamma((dim+2)/(2*s))/(dim*gamma(dim/(2*s)));
    C = 2^(dim/2)*gamma(dim/2); 
    A = C*((s*B^(dim/2))/(gamma(dim/(2*s))))*(((2*pi)^dim)*det(P))^(-1/2);
    
    for i=1:N
        for j=1:N
            xv=[x(i); y(j)];
            p(j,i) = A*exp(-(B*(xv-xvbar)'*inv(P)*(xv-xvbar))^s);
            v = (-2*lambda*s)*B*(B*(xv-xvbar)'*inv(P)*(xv-xvbar))^(s-1)*inv(P)*(xv-xvbar);  
            v_x(j,i) = v(1,1); v_y(j,i) = v(2,1); 
        end
    end

    v_x = padarray(v_x,[1 1],'replicate','both');
    v_y = padarray(v_y,[1 1],'replicate','both');

    p = reshape(p,[N*N,1]);
    p = boundary_conditions(p,N,d_x,d_y); 
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = drone_M(N,target,p_drones_sum)
    obs = observations(target,p_drones_sum,N); 

    B = cell(1,N); 
    % Initializing block B matrices
    for i=1:N
        B_i = zeros(N,N);
        B2 = zeros(1,N); 
        for j=1:N
            B2(1,j) = -obs(i,j);
        end
        B_i = B_i + diag(B2,0);
        B{1,i} = B_i; 
    end
    
    M = sparse(blkdiag(B{:})); % Main Block Diagonal
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = M(d_x,d_y,N,target)
    B = cell(1,N); A = cell(1,N-1); C = A; % Cell array of A,B,C matrices
    
    % Initializing block B matrices
    for i=1:N
        B_i = zeros(N,N);
        B2 = zeros(1,N); 
        for j=1:N
            B2(1,j) = -((target.v_x(j+1,i+2)-target.v_x(j+1,i))/(2*d_x)+(target.v_y(j+2,i+1)-target.v_y(j,i+1))/(2*d_y))...
                -((target.dif(j+1,i+2)-4*target.dif(j+1,i+1)+target.dif(j+1,i))/(d_x^2)+(target.dif(j+2,i+1)-4*target.dif(j+1,i+1)+target.dif(j,i+1))/(d_y^2));
        end
        B_i = B_i + diag(B2,0);
        B{1,i} = B_i; 
    end
    
    for i=1:N
        B_i = B{1,i};
        B3 = zeros(1,N-1); B1 = B3;  
        for j=1:(N-1)
            B1(1,j) = (target.v_y(j+2,i+1)/(2*d_y))+((-target.dif(j+2,i+1)+target.dif(j+1,i+1)+target.dif(j,i+1))/(d_y^2));
            B3(1,j) = (-target.v_y(j+1,i+1)/(2*d_y))+((target.dif(j+2,i+1)+target.dif(j+1,i+1)-target.dif(j,i+1))/(d_y^2));
        end
        B_i = B_i + diag(B3,1)+diag(B1,-1);
        B{1,i} = B_i; 
    end
    
    for i=1:N-1
        A_i = zeros(N,N); C_i = A_i;
        A1 = zeros(1,N); C1 = A1; 
        for j=1:(N)
            A1(1,j) = (target.v_x(j+1,i+2)/(2*d_x))+((-target.dif(j+1,i+2)+target.dif(j+1,i+1)+target.dif(j+1,i))/(d_x^2));
            C1(1,j) = (-target.v_x(j+1,i+1)/(2*d_x))+((target.dif(j+1,i+2)+target.dif(j+1,i+1)-target.dif(j+1,i))/(d_x^2));
        end
        A_i = A_i + diag(A1,0); C_i = C_i + diag(C1,0);
        A{1,i} = A_i; C{1,i} = C_i; 
    end
    
    M = blkdiag(B{:}); % Main Block Diagonal
    M_A = blkdiag(A{:}); M_A = [[zeros(N,N*N-N);M_A] zeros(N*N,N)]; % Lower Block Diagonal
    M_C = blkdiag(C{:}); M_C = [zeros(N*N,N) [M_C;zeros(N,N*N-N)]]; % Upper Block Diagonal
    M = sparse(M + M_A + M_C); % Complete M Matrix
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function p=search_gauss(N,x,y,q,sigma,a)
    for i=1:N
        for j=1:N
            xv=[x(i); y(j)];
            p(j,i) = a*exp(-((xv-q)'*(xv-q))/sigma^2);
        end
    end
    p = reshape(p,[N*N,1]);
end % end function gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,v_x,v_y]=target_numer(N,x,y,d_x,d_y,lambda,shp,X_mesh,Y_mesh,flag)
p = zeros(N,N); dt = 0.0005;

if(flag==1) % lambda-convex hull approximation
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
    for i=1:200
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
    v_x = padarray(v_x,[1 1],'replicate','both');
    v_y = padarray(v_y,[1 1],'replicate','both');

    p = reshape(p,[N*N,1]);
    p = boundary_conditions(p,N,d_x,d_y); 

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
    v_x = padarray(v_x,[1 1],'replicate','both');
    v_y = padarray(v_y,[1 1],'replicate','both');

    p = reshape(p,[N*N,1]);
    p = boundary_conditions(p,N,d_x,d_y); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drones = update_pos(drones, t)
    if(drones.orbit_flag==1)
        for j=1:drones.num
            drones.pos(j,1) = drones.radius(j)*cos(drones.ang_speed(j)*t + drones.init_theta(j));
            drones.pos(j,2) = drones.radius(j)*sin(drones.ang_speed(j)*t + drones.init_theta(j));
        end
    elseif(drones.orbit_flag==2)
        for j=1:drones.num
            drones.pos(j,1) = sqrt(drones.focal^2*cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j)))+sqrt(drones.focal^4*(cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j))).^2+drones.b(j)^4-drones.focal^4))).*cos(drones.ang_speed(j)*t+drones.init_theta(j)); 
            drones.pos(j,2) = drones.f.*sqrt(drones.focal^2*cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j)))+sqrt(drones.focal^4*(cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j))).^2+drones.b(j)^4-drones.focal^4))).*sin(drones.ang_speed(j)*t+drones.init_theta(j));
        end
    elseif(drones.orbit_flag==3)
        for j=1:drones.num
            Xm = cos(drones.ang_speed(j)*t+drones.init_theta(j)); Ym = drones.f*sin(2*(drones.ang_speed(j)*t+drones.init_theta(j)));
            drones.pos(j,1) = drones.scale.*(drones.c(j)*Xm - drones.s(j)*Ym);
            drones.pos(j,2) = drones.scale.*(drones.s(j)*Xm + drones.c(j)*Ym); 
        end
    elseif(drones.orbit_flag==4)
        for j=1:drones.num
            theta = drones.init_theta(j) + drones.ang_speed(j)*t; 
            b=-2*drones.A(j)^2*cos(4*theta); c=drones.A(j)^4-drones.B(j)^4; r=sqrt((-b+sqrt(b^2-4*c))/2);
            drones.pos(j,1) = r*sin(theta+(drones.ang_speed(j)*t*drones.v));
            drones.pos(j,2) = r*cos(theta+(drones.ang_speed(j)*t*drones.v));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tot_v_x, tot_v_y] = update_vel(drones, x, y, N, target)
    tot_v_x = target.relax_v_x; tot_v_y = target.relax_v_y; dim=2; s=0.6;  
    
    B = gamma((dim+2)/(2*s))/(dim*gamma(dim/(2*s)));
    C = 2^(dim/2)*gamma(dim/2); 
    A = C*((s*B^(dim/2))/(gamma(dim/(2*s))))*(((2*pi)^dim)*det(drones.P))^(-1/2);
    
    for k=1:drones.num
        grad_p_x = zeros(N,N); grad_p_y = grad_p_x;
        for i=2:N+1
            for j=2:N+1
                xv=[x(i-1); y(j-1)]; q = [drones.pos(k,1); drones.pos(k,2)]; 
                p = A*exp(-(B*(xv-q)'*inv(drones.P)*(xv-q))^s);
                grad_p = -2*s*p*B*(B*(xv-q)'*inv(drones.P)*(xv-q))^(s-1)*inv(drones.P)*(xv-q);  
                grad_p_x(j-1,i-1) = grad_p(1,1); grad_p_y(j-1,i-1) = grad_p(2,1); 
            end
        end
        tot_v_x(2:N+1,2:N+1) = tot_v_x(2:N+1,2:N+1) - drones.d*grad_p_x; 
        tot_v_y(2:N+1,2:N+1) = tot_v_y(2:N+1,2:N+1) - drones.d*grad_p_y; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dif = update_dif(drones, x, y, N)
    dif = zeros(N,N); dim=2; s=0.6; 
    
    B = gamma((dim+2)/(2*s))/(dim*gamma(dim/(2*s)));
    C = 2^(dim/2)*gamma(dim/2); 
    A = C*((s*B^(dim/2))/(gamma(dim/(2*s))))*(((2*pi)^dim)*det(drones.P))^(-1/2);
    
    for k=1:drones.num
        p = zeros(N,N);
        for i=1:N
            for j=1:N
                xv=[x(i); y(j)]; q = [drones.pos(k,1); drones.pos(k,2)]; 
                p(j,i) = A*exp(-(B*(xv-q)'*inv(drones.P)*(xv-q))^s);
            end
        end
        dif = dif + p;  
    end

    dif = padarray(dif,[1 1],'replicate','both');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function create_video(F, sim, title)
    writerObj = VideoWriter(title, 'MPEG-4');
    writerObj.FrameRate = sim.FrameRate;
    open(writerObj);
    for i=1:length(F)
        if((i==1)||(i==sim.timesteps))
            for j=1:(sim.FrameRate*3)
                frame = F(i);
                writeVideo(writerObj, frame);
            end
        end
        frame = F(i);
        writeVideo(writerObj, frame);    
    end
    close(writerObj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obs = observations(target, p_drones_sum, N)
    obs = zeros([N,N]); p = reshape(target.p,[N,N]); p_obs = reshape(p_drones_sum,[N,N]);
    for i=1:N
        for j=1:N
            obs(i,j) = p(j,i)*p_obs(j,i);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%