% create_dataset.m 
%
% This code randomly generates the alpha-convex hull, den, and lake
% datasets to be used by the relaxation_advection, and (non)evasive search
% codes. 
%
% By Benjamin L. Hanson, 2024

n=50000;
x1=2*rand(n,1)-1; y1=2*rand(n,1)+x1.*x1-1.5;

s=0.1; shp = alphaShape(x1,y1,s);
figure(1); clf; 
plot(shp,'EdgeColor','none'); hold on; plot(x1,y1,'kx','MarkerSize',12);
axis normal; axis([-1 1 -1.5 1.5]);
save('./Datasets/alpha.mat', 'shp')

mu1 = [-0.75 -0.5]; sigma1 = [0.1 0; 0 0.1]; r1 = 1.5*sqrt(0.1); 
mu2 = [0.5 -0.3]; sigma2 = [0.1 0; 0 0.1]; r2x = 0.2; r2y = 0.6;  

count = 0;
while (count < (n/4))
    sample = mvnrnd(mu1, sigma1, 1);
    x_i = sample(1,1); y_i = sample(1,2); 
    if(inShape(shp,x_i,y_i) == 1)
        x1(end+1,1) = x_i; 
        y1(end+1,1) = y_i; 
        count = count + 1; 
    end
end
s=0.4; shp = alphaShape(x1,y1,s);
save('./Datasets/kde_den.mat', 'shp')

theta = linspace(0, 2*pi, 1000);
circ_x1 = []; circ_y1 = []; 
circ_x2 = []; circ_y2 = []; 

for i=1:length(theta)
    x = mu1(1,1) + r1 * cos(theta(i));
    y = mu1(1,2) + r1 * sin(theta(i));
    if(inShape(shp,x,y) == 1)
        circ_x1(end+1,1) = x;
        circ_y1(end+1,1) = y; 
    end
    
    x = mu2(1,1) + r2x * cos(theta(i));
    y = mu2(1,2) + r2y * sin(theta(i));
    circ_x2(end+1,1) = x;
    circ_y2(end+1,1) = y; 
end

s=0.4; circ_shp = alphaShape(circ_x2,circ_y2,s);
lake_x = []; lake_y = [];

for i=1:n
    if(inShape(circ_shp,x1(i),y1(i)) == 0)
        lake_x(end+1,1) = x1(i);
        lake_y(end+1,1) = y1(i); 
    end
end
s=0.4; shp = alphaShape(lake_x,lake_y,s);
save('./Datasets/kde_lake.mat', 'shp')

figure(2); clf; hold on;
plot(x1,y1,'kx','MarkerSize',12);
scatter(-0.75,-0.5,50,'r','filled');
scatter(circ_x1,circ_y1,10,'r', 'filled')
axis normal; axis([-1 1 -1.5 1.5]);

figure(3); clf; hold on;
plot(lake_x,lake_y,'kx','MarkerSize',12);
scatter(0.5,-0.3,50,'g','filled');
scatter(circ_x2,circ_y2,10,'g', 'filled')
axis normal; axis([-1 1 -1.5 1.5]);