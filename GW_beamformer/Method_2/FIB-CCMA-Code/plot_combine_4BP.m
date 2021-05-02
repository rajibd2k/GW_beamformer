function plot_combine_4BP(theta,f_vect,BP_matrix1,BP_matrix2,BP_matrix3,BP_matrix4)
% this is a function to plot the 2D polar beampattern, 
% 3D polar beampattern (vers different frequency), DF, 
% and WNG in the same figure.
% 
% by Gongping Huang and Jilu Jin
% Dec.17, 2018


% beampatterns -> dB plot
% theta -> rad
% plot beampatterns 

FontSize    = 12;

%% plot four figures
figure('position', [50, 20, 720, 640]);
% axes #1
axes('position', [0.07, 0.55, 0.40, 0.40]);
% # for other method, change the following code
plot_3d_bp_f(BP_matrix1, f_vect, theta)
text(-0.85, -1.3,'(a)','fontsize',FontSize,'interpreter', 'latex')

% axes #2
axes('position', [0.57, 0.55, 0.40, 0.40]);
plot_3d_bp_f(BP_matrix2, f_vect, theta)
text(-0.85, -1.3, '(b)', 'interpreter', 'latex', 'fontSize', 12);

% axes #3 % plot directive factor
axes('position', [0.07, 0.08, 0.40, 0.40]);
plot_3d_bp_f(BP_matrix3, f_vect, theta)
text(-0.85, -1.3,'(c)','fontsize',FontSize,'interpreter', 'latex')


% axes #4 % plot directive factor
axes('position', [0.57, 0.08, 0.40, 0.38]);
plot_3d_bp_f(BP_matrix4, f_vect, theta)
text(-0.85, -1.3, '(d)','fontsize',FontSize,'interpreter', 'latex')


function plot_3d_bp_f(BP_matrix, f_vect, theta)
%   plot polar beampattern vers frequency
%   
%   by: Wenxing Yang, Jilu Jin, ans Gongping Huang
%   Dec.17, 2018

% set some parameters
fontSize        = 12;
rMin            = -50;

BP_matrix = 20 * log10(BP_matrix ./ max(BP_matrix));
BMat = BP_matrix;
fBin = f_vect(:) / 1e3;
theta = theta(:);

%% set some parameters

len_f = length(fBin);
len_theta = length(theta);

% express the gain between -40~0 dB
BMat((BMat(:) < rMin)) = rMin;  
BMat = (BMat - rMin) ./ max(BMat - rMin);  

%% beampatterns polar plot

X = BMat .* cos(theta);
Y = BMat .* sin(theta);
Z = repmat(fBin, 1, length(theta)).';
R = sqrt(X.^2 + Y.^2);

blindZone = 0;
BPs = 160;
BPe = BPs + 181;

%% plot 3d beampattern

colormap(jet(64));

decimator = 3;
fBP = round(linspace(1, len_f, (len_f)/ decimator + 1));
if blindZone
    h = surf(X(BPs:BPe, fBP), Y(BPs:BPe, fBP),...
        Z(BPs:BPe, fBP), R(BPs:BPe, fBP));
else
    h = surf(X(:, fBP), Y(:, fBP), Z(:, fBP), R(:, fBP));
end

set(h, 'edgecolor', 'none')
axis([-1, 1, -1, 1, min(fBin), max(fBin)]);
set(gca, 'xTick', []);
set(gca, 'yTick', []);

zlabel('$f$ (kHz)', 'fontsize', fontSize, 'interpreter', 'latex')
set(gca, 'fontsize', fontSize, 'fontname', 'Times new roman');
view([-30, 20]);

%% plot circles

Cr = 0.25:0.25:1;
Xc = zeros(len_theta, 2);
Yc = zeros(len_theta, 2);
Zc = repmat([fBin(1), fBin(end)], len_theta, 1);

hold on;
for rIndex = 1:length(Cr)
    r = Cr(rIndex);
    Xc(:, 1) = r * cos(theta);
    Xc(:, 2) = Xc(:, 1);
    Yc(:, 1) = r * sin(theta);
    Yc(:, 2) = Yc(:, 1);

    if rIndex == length(Cr)
        plot3(Xc, Yc, Zc, 'k-', 'lineWidth', 0.2);
    else
        plot3(Xc, Yc, Zc, '-.', 'lineWidth', 0.1 , 'color', [0.7, 0.7, 0.7]);
    end
end
hold off;

%% plot spokes

 th = (1:12)'*2*pi/12;
 len_th = length(th);
 r_range = (0:1)';
 len_r_range = length(r_range);
 x = zeros(len_r_range,1);
 y = zeros(len_r_range,1);
 z_range = [1; ceil(len_f/2); ceil(len_f)];
    
hold on;
for z_index = 1:2:3
    z = ones(len_r_range,1)*fBin(z_range(z_index));   
    for th_index = 1:len_th
        x = r_range*cos(th(th_index));
        y = r_range*sin(th(th_index));   
        plot3(x, y, z, '-.', 'lineWidth', 0.1, 'color', [0.7, 0.7, 0.7]);
    end
end
% shading interp; 

%% plot others

% top and bottom
fill3(X(:, end), Y(:, end), Z(:, end), 'w');
fill3(X(:, 1), Y(:, 1), Z(:, 1), 'w');

% 2-D beampattern
plot3(X(:, end), Y(:, end), Z(:, end), 'b-', 'lineWidth', 1.2);
plot3(X(:, 1), Y(:, 1), Z(:, 1), 'b-', 'lineWidth', 1.2);
box on; 

hold off;