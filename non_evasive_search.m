clear all; close; clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%
sim.T = 2; sim.dt = 0.001; sim.timesteps = round(sim.T/sim.dt); sim.FrameRate = round(.1/sim.dt);
theta = linspace(0,2*pi, 1000); 
N = 101; 
L_x = 3; L_y = 4;  
d_x = L_x/(N-1); d_y = L_y/(N-1);
x = linspace(-L_x/2,L_x/2,N); y = linspace(-L_y/2,L_y/2,N);
[X_mesh, Y_mesh] = meshgrid(x,y);
colors = ['g' 'b' 'r' 'y' 'm' 'c' 'k'];
%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
target.lambda = 0.25; 
target.stats_flag = 2; %1: Gaussian Statistics, 2: Numerical Territory
if(target.stats_flag == 1)
    target.s = 2; 
    target.xvbar=[0; 0]; target.P = [2 0.8; 0.8 2];
else
    flag = 1; % 1:alpha, 2:KSD Den, 3:KSD Lake
    if(flag==1)
        load ./Datasets/cougar_alpha.mat
    elseif(flag==2)
        load ./Datasets/cougar_ksd_den.mat
    elseif(flag==3)
        load ./Datasets/cougar_ksd_lake.mat
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Target Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drone Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drones.num = 3; 
drones.ang_speed = [8 8 8]; 
drones.init_theta = [0 2*pi/3 4*pi/3];

drones.orbit_flag = 3; %1: Concentric orbits, 2: Cassini ovals, 3: Lemniscates
if (drones.orbit_flag == 1)
    drones.radius = [.5 1.5 2.5]; 
    drones.sigma = 0.25; drones.a = 10/sim.dt; 
elseif (drones.orbit_flag == 2)
    drones.focal = 1.2;
    drones.b = [1.22 1.4 1.7];
    drones.f = 1.41; 
    drones.sigma = 0.25; drones.a = 10/sim.dt; 
elseif (drones.orbit_flag == 3)
    drones.phi = pi.*([1:drones.num]-1)./drones.num; 
    drones.c = cos(drones.phi); 
    drones.s = sin(drones.phi); 
    drones.f = 1/(2*drones.num);
    drones.Xm = cos(theta); drones.Ym = drones.f*sin(2.*theta); 
    drones.scale = 1.2; 
    drones.sigma = 0.1; drones.a = 8/sim.dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drone Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Target Density/Relaxation Advection
if(target.stats_flag == 1)
    [target.p,target.v_x,target.v_y]=target_gauss(N,x,y,target.xvbar,target.P,target.s,target.lambda,d_x,d_y);
else
    [target.p,target.v_x,target.v_y]=target_cougar(N,x,y,d_x,d_y,target.lambda,shp,X_mesh,Y_mesh,flag);
end

fig = figure(1); clf; hold on; fig.Position = [150 150 1200 600];
sgtitle(['Pre-search, iter = 0, t = 0, \Delta', 't = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda)]);

subplot(1,2,1); hold on; 
xlabel('x','FontSize',16,'Interpreter','latex');
ylabel('y','FontSize',16,'Interpreter','latex');
zlabel('$p$','FontSize',16,'Interpreter','latex');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
%axis square
view(0,90); colorbar;

subplot(1,2,2);  hold on; 
xlabel('x','FontSize',16,'Interpreter','latex');
ylabel('y','FontSize',16,'Interpreter','latex');
zlabel('$p$','FontSize',16,'Interpreter','latex');
xlim([x(1),x(end)]);
ylim([y(1),y(end)]);
%axis square
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
        x_unrotated = sqrt(drones.focal^2*cos(2.*theta)+sqrt(drones.focal^4*(cos(2.*theta).^2+drones.b(i)^4-drones.focal^4))).*cos(theta); 
        y_unrotated = drones.f.*sqrt(drones.focal^2*cos(2.*theta)+sqrt(drones.focal^4*(cos(2.*theta).^2+drones.b(i)^4-drones.focal^4))).*sin(theta);
        x_orbit = x_unrotated*cos(pi/4) - y_unrotated*sin(pi/4);
        y_orbit = x_unrotated*sin(pi/4) + y_unrotated*cos(pi/4); 
        x_pos = sqrt(drones.focal^2*cos(2.*drones.init_theta(i))+sqrt(drones.focal^4*(cos(2.*drones.init_theta(i)).^2+drones.b(i)^4-drones.focal^4))).*cos(drones.init_theta(i)); 
        y_pos = drones.f.*sqrt(drones.focal^2*cos(2.*drones.init_theta(i))+sqrt(drones.focal^4*(cos(2.*drones.init_theta(i)).^2+drones.b(i)^4-drones.focal^4))).*sin(drones.init_theta(i));
        drones.pos(i,1) = x_pos*cos(pi/4) - y_pos*sin(pi/4);
        drones.pos(i,2) = x_pos*sin(pi/4) + y_pos*cos(pi/4); 
    elseif(drones.orbit_flag == 3)
        x_orbit = drones.scale.*(drones.c(i).*drones.Xm - drones.s(i).*drones.Ym); 
        y_orbit = drones.scale.*(drones.s(i).*drones.Xm + drones.c(i).*drones.Ym)-.25;
        Xm = cos(drones.init_theta(i)); Ym = drones.f*sin(2*drones.init_theta(i));
        drones.pos(i,1) = drones.scale.*(drones.c(i)*Xm - drones.s(i)*Ym); 
        drones.pos(i,2) = drones.scale.*(drones.s(i)*Xm + drones.c(i)*Ym)-.25; 
    end

    subplot(1,2,1); 
    plot(x_orbit,y_orbit,colors(i),'Linewidth',2,'HandleVisibility','off'); 
    plots_1(i+1) = scatter(drones.pos(i,1),drones.pos(i,2),75,colors(i),'filled', 'HandleVisibility', 'off');
    subplot(1,2,2); 
    plot(x_orbit,y_orbit,colors(i),'Linewidth',2,'HandleVisibility','off'); 
    plots_2(i+1) = scatter(drones.pos(i,1),drones.pos(i,2),75,colors(i),'filled', 'HandleVisibility', 'off');
end
frames(1) = getframe(gcf); 
drawnow 

pause(3); delete(plots_1); delete(plots_2);
target.relax_M = M(d_x,d_y,N,target); 

t = 0; 
for i=1:sim.timesteps
    drones = update_pos(drones,t);
    
    p_drones_sum = zeros([N*N,1]);
    for j=1:drones.num
        p_drones_sum = p_drones_sum + search_gauss(N,x,y,[drones.pos(j,1);drones.pos(j,2)],drones.sigma,drones.a);
        subplot(1,2,1); 
        plots_1(j+1) = scatter(drones.pos(j,1),drones.pos(j,2),75,colors(j),'filled', 'HandleVisibility', 'off');
        subplot(1,2,2); 
        plots_2(j+1) = scatter(drones.pos(j,1),drones.pos(j,2),75,colors(j),'filled', 'HandleVisibility', 'off');
    end

    target.M = target.relax_M + drone_M(N,target,p_drones_sum); 

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
    
    sgtitle(['Non-evasive searching, iter = ',num2str(i), ', t = ', num2str(t), ', \Delta t = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda), ', a = ', num2str(drones.a*sim.dt), '/\Delta','t', ', \sigma = ', num2str(drones.sigma)]);
    frames(i+1) = getframe(gcf); 
    drawnow; 
    if(i ~= sim.timesteps)
        delete(plots_1); delete(plots_2);  
    end
    t = t + sim.dt; 
end

pause(3); delete(plots_1); delete(plots_2);

for i=sim.timesteps+1:sim.timesteps*2
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
    sgtitle(['Post-search, iter = ',num2str(i), ', t = ', num2str(t), ', \Delta', 't = ', num2str(sim.dt), ', \lambda = ', num2str(target.lambda)]);
    frames(i+1) = getframe(gcf);
    drawnow; 
    if(i ~= sim.timesteps*2)
        delete(plot4_1); delete(plot4_2);
    end
    t = t + sim.dt; 
end

%create_video(frames, sim, ['./Figures/Non-evasive/non-evasive_search_',num2str(target.lambda),'.mp4']);
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
                -(((2*target.lambda)/(d_x^2))+((2*target.lambda)/(d_y^2)));
        end
        B_i = B_i + diag(B2,0);
        B{1,i} = B_i; 
    end
    
    for i=1:N
        B_i = B{1,i};
        B3 = zeros(1,N-1); B1 = B3;  
        for j=1:(N-1)
            B1(1,j) = (target.v_y(j+2,i+1)/(2*d_y))+(target.lambda/d_y^2);
            B3(1,j) = (-target.v_y(j+1,i+1)/(2*d_y))+(target.lambda/d_y^2);
        end
        B_i = B_i + diag(B3,1)+diag(B1,-1);
        B{1,i} = B_i; 
    end
    
    for i=1:N-1
        A_i = zeros(N,N); C_i = A_i;
        A1 = zeros(1,N); C1 = A1; 
        for j=1:(N)
            A1(1,j) = (target.v_x(j+1,i+2)/(2*d_x))+(target.lambda/(d_x^2));
            C1(1,j) = (-target.v_x(j+1,i+1)/(2*d_x))+(target.lambda/(d_x^2));
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
function [p,v_x,v_y]=target_cougar(N,x,y,d_x,d_y,lambda,shp,X_mesh,Y_mesh,flag)
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
            x_pos = sqrt(drones.focal^2*cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j)))+sqrt(drones.focal^4*(cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j))).^2+drones.b(j)^4-drones.focal^4))).*cos(drones.ang_speed(j)*t+drones.init_theta(j)); 
            y_pos = drones.f.*sqrt(drones.focal^2*cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j)))+sqrt(drones.focal^4*(cos(2.*(drones.ang_speed(j)*t+drones.init_theta(j))).^2+drones.b(j)^4-drones.focal^4))).*sin(drones.ang_speed(j)*t+drones.init_theta(j));
            drones.pos(j,1) = x_pos*cos(pi/4) - y_pos*sin(pi/4);
            drones.pos(j,2) = x_pos*sin(pi/4) + y_pos*cos(pi/4); 
        end
    elseif(drones.orbit_flag==3)
        for j=1:drones.num
            Xm = cos(drones.ang_speed(j)*t+drones.init_theta(j)); Ym = drones.f*sin(2*(drones.ang_speed(j)*t+drones.init_theta(j)));
            drones.pos(j,1) = drones.scale.*(drones.c(j)*Xm - drones.s(j)*Ym);
            drones.pos(j,2) = drones.scale.*(drones.s(j)*Xm + drones.c(j)*Ym)-.25; 
        end
    end
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