%%% MIE301 Lab 2
%%
close all;                                      % closes all figures
clear all;                                      % clears all variables from memory
clc;                                            % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
xmin= -200;                                     % leftmost window edge
xmax= 300;                                      % rightmost window edge
ymin= -20;                                      % bottom window edge
ymax= 400;                                      % top window edge

% xmin= -100;                                   % leftmost window edge
% xmax= 200;                                    % rightmost window edge
% ymin= 100;                                    % bottom window edge
% ymax= 200;                                    % top window edge

%% Define Angles
increments = 91; % number of theta_2 configuration steps to calculate along mechanism rotation
deg = pi/180;
min_rot_2 = 170*deg;
max_rot_2 = 200*deg;
avg_rot_2 = (max_rot_2+min_rot_2)/2;

min_rot_4 = 150*deg;
max_rot_4 = 220*deg;

min_rot_11 = 20*deg;
%max_rot_11 = 690*deg;
max_rot_11 = 1000*deg;

onez = ones(1,increments);

temp_2a= linspace(min_rot_2,max_rot_2,increments);
temp_2b= linspace(min_rot_2,avg_rot_2,increments);
temp_4 = linspace(min_rot_4,max_rot_4,increments);
temp_11= linspace(max_rot_11,min_rot_11,increments);


% theta_2 = [temp_2a flip(temp_2a) temp_2b avg_rot_2*([onez onez onez])];
% theta_4 = [min_rot_4*([onez onez onez]) temp_4 max_rot_4*([onez onez])];
% theta_11= [max_rot_11*([onez onez onez onez]) temp_11 flip(temp_11)];

% Height and Canopy Angle
% theta_2 = [temp_2a flip(temp_2a) min_rot_2*([onez onez])];
% theta_4 = [min_rot_4*([onez onez]) temp_4 flip(temp_4)];
% theta_11= [max_rot_11*([onez onez onez onez])];

% Just canopy angle
theta_2 = [180*deg*([onez])];
theta_4 = [temp_4];
theta_11= [max_rot_11*([onez])];

%just height
% theta_2 = [temp_2a flip(temp_2a)];
% theta_4 = [min_rot_4*([onez onez])];
% theta_11= [max_rot_11*([onez onez])];

% Just Opening 
% theta_2 = [min_rot_2*onez];
% theta_4 = [min_rot_4*onez];
% theta_11= [temp_11];

% t= linspace (0,4,2*increments);
t= linspace (0,4,length(theta_2));
%% Define Link Lengths and ratios
L1  = 200;
L2  = 30;
L3  = 100;
L4  = 25;
L5  = 75;
L6  = 15;                            %length of EJ
L7  = 150;
L8  = L7;
L9  = 35;
L10A= 100;
L10B= 15;
L10C= 10;

L13 = 130;                           %canopy radius

OAB_ratio = L3/L1;                  % AB / OB
FDE_ratio = L5/L3;                  % FE / DE
ADC_ratio = 0.4;
GH_ratio = 0.4;                     % GH / GF
VT_ratio = 0.15;
MN_ratio = 2.1;

L11 = VT_ratio * L13;
L12 = L11;
theta_IEJ = 60*deg;                 % Angle IEJ
lead = 5/pi;                        % Centimeters of lead screw travel per radian. 
canopy_angle = theta_4-theta_IEJ;

g=9.81;
linear_weight_density = 0.0015*g;%N/cm          %Linear density of all the links except 3 and 5
lead_linear_density = 0.005;

w1 = linear_weight_density * L1;                                                                                   
w2 = linear_weight_density * L2;                                              
w3 = 16*g;                                               
w4 = linear_weight_density * L4;                                                
w5 = 33.302147602117415;                                                
w6 = linear_weight_density * (L6 + L4*(1-GH_ratio) + sqrt(L6^2 + (L4*(1-GH_ratio))^2- 2*L6*L4*(1-GH_ratio)*cos(theta_IEJ)));                                                
w7 = linear_weight_density * L7;                                                
w8 = linear_weight_density * L8;                                               
w9 = lead_linear_density * L9;                                                
w10A=linear_weight_density * (L10A);
w10B=linear_weight_density * (L10B);  
w10C=linear_weight_density * (L10C);             
w11 =linear_weight_density * L11;                                              
w12 =linear_weight_density * L12;                                               
w13 =linear_weight_density * L13;
w14 =w13;
         
w_tot = w1 + w2 + w3 + w4 + w5 + w6 +w7 + w8 + w9 + w10A +w10B + w10C + 4*w11 + 4*w12 + 8*w13;
W_canopy = w9 + 4*w11 + 4*w12 + 8*w13;

%% set up figure
figure(1);                         %create new figure
set(1,'WindowStyle','Docked')      %dock the figure
clr = 'Color';
LW = 'LineWidth';
clr1 = [0 0 0];
clr3 = [0.69 0.1 0.1];
clr4 = [0.84 0.6 0.6];
dot_clr = 'm.';
grnd_clr = [0 0 0];

animate = true;
%% Set up zero vectors just to speed up simulation, Ignore this section
if true
    Ox = zeros(1,length(theta_2));
    Oy = zeros(1,length(theta_2));
    Ax = zeros(1,length(theta_2));
    Ay = zeros(1,length(theta_2));
    Bx = zeros(1,length(theta_2));
    By = zeros(1,length(theta_2));
    Cx = zeros(1,length(theta_2));
    Cy = zeros(1,length(theta_2));
    Dx = zeros(1,length(theta_2));
    Dy = zeros(1,length(theta_2));
    Ex = zeros(1,length(theta_2));
    Ey = zeros(1,length(theta_2));
    Fx = zeros(1,length(theta_2));
    Fy = zeros(1,length(theta_2));
    Gy = zeros(1,length(theta_2));
    Gx = zeros(1,length(theta_2));
    Hx = zeros(1,length(theta_2));
    Hy = zeros(1,length(theta_2));
    Ix = zeros(1,length(theta_2));
    Iy = zeros(1,length(theta_2));
    Jx = zeros(1,length(theta_2));
    Jy = zeros(1,length(theta_2));
    Kx = zeros(1,length(theta_2));
    Ky = zeros(1,length(theta_2));
    Lx = zeros(1,length(theta_2));
    Ly = zeros(1,length(theta_2));
    Mx = zeros(1,length(theta_2));
    My = zeros(1,length(theta_2));
    Nx = zeros(1,length(theta_2));
    Ny = zeros(1,length(theta_2));
    Px = zeros(1,length(theta_2));
    Py = zeros(1,length(theta_2));
    Qx = zeros(1,length(theta_2));
    Qy = zeros(1,length(theta_2));
    Rx = zeros(1,length(theta_2));
    Ry = zeros(1,length(theta_2));
    Sx = zeros(1,length(theta_2));
    Sy = zeros(1,length(theta_2));
    Tx = zeros(1,length(theta_2));
    Ty = zeros(1,length(theta_2));
    Ux = zeros(1,length(theta_2));
    Uy = zeros(1,length(theta_2));
    Vx = zeros(1,length(theta_2));
    Vy = zeros(1,length(theta_2));
    Wx = zeros(1,length(theta_2));
    Wy = zeros(1,length(theta_2));
end

%% calculate mechanism motion
for i=1:length(theta_2)                               % step through motion of the mechanism
    
    % Set up Coordinates:
    Ox(i) = 0;                                  
    Oy(i) = 0;                                  
    Ax(i) = 0;                                                                     
    Ay(i) = L1 - L3;                                                               
    Bx(i) = 0;                                                  
    By(i) = L1;                                                     
    Cx(i) = Ax(i)+L2*cos(theta_2(i));                                    
    Cy(i) = Ay(i)+L2*sin(theta_2(i));                                    
    Dx(i) = ADC_ratio*Ax(i)+(1-ADC_ratio)*Cx(i);                                     
    Dy(i) = ADC_ratio*Ay(i)+(1-ADC_ratio)*Cy(i);                                     
    Ex(i) = Dx(i);                                     
    Ey(i) = Dy(i)+L3;                                     
    Fx(i) = Dx(i);                                     
    Fy(i) = Dy(i)+(1-FDE_ratio)*L3;                                
    Gy(i) = Fy(i)+L4*sin(theta_4(i));                                     
    Gx(i) = Fx(i)+L4*cos(theta_4(i));                                    
    Hx(i) = GH_ratio*Fx(i)+(1-GH_ratio)*Gx(i);                                     
    Hy(i) = GH_ratio*Fy(i)+(1-GH_ratio)*Gy(i);  
    Ix(i) = Hx(i);                                                                   
    Iy(i) = Hy(i)+L5;                                                                
    Jx(i) = Ex(i)+L6*cos(canopy_angle(i));                                     
    Jy(i) = Ey(i)+L6*sin(canopy_angle(i));                                      
    Kx(i) = Ex(i)-L8*cos(theta_2(i));                                     
    Ky(i) = Ey(i)-L8*sin(theta_2(i));                                     
    Lx(i) = Kx(i)+Jx(i)-Ex(i);                                     
    Ly(i) = Ky(i)+Jy(i)-Ey(i);                                     
    Mx(i) = Lx(i)-L9*cos(canopy_angle(i));                                     
    My(i) = Ly(i)-L9*sin(canopy_angle(i));                                     
    Nx(i) = Lx(i) + MN_ratio*(Mx(i)-Lx(i));                                  
    Ny(i) = Ly(i) + MN_ratio*(My(i)-Ly(i));
    Px(i) = Nx(i);                                     
    Py(i) = Ny(i)-L10A;                                     
    Qx(i) = Px(i)+L11*cos(theta_11(i));                                     
    Qy(i) = Py(i);                                     
    Rx(i) = Qx(i);                                    
    Ry(i) = Qy(i)-L10C;
    Sx(i) = Mx(i) - (2*L13*VT_ratio-lead*theta_11(i))*cos(canopy_angle(i));%
    Sy(i) = My(i) - (2*L13*VT_ratio-lead*theta_11(i))*sin(canopy_angle(i)); %                                    
    Tx(i) = ((Sx(i) + Mx(i))/2)+sqrt((L13*VT_ratio)^2-(L13*VT_ratio-0.5*lead*theta_11(i))^2)*sin(canopy_angle(i));%
    Ty(i) = ((Sy(i) + My(i))/2)-sqrt((L13*VT_ratio)^2-(L13*VT_ratio-0.5*lead*theta_11(i))^2)*cos(canopy_angle(i));%
    Ux(i) = ((Sx(i) + Mx(i))/2)-sqrt((L13*VT_ratio)^2-(L13*VT_ratio-0.5*lead*theta_11(i))^2)*sin(canopy_angle(i));%
    Uy(i) = ((Sy(i) + My(i))/2)+sqrt((L13*VT_ratio)^2-(L13*VT_ratio-0.5*lead*theta_11(i))^2)*cos(canopy_angle(i));%%
    Vx(i) = Mx(i) + (1/VT_ratio)*(Tx(i)-Mx(i)); %                                    
    Vy(i) = My(i) + (1/VT_ratio)*(Ty(i)-My(i)); %
    Wx(i) = Mx(i) + (1/VT_ratio)*(Ux(i)-Mx(i)); %                                  
    Wy(i) = My(i) + (1/VT_ratio)*(Uy(i)-My(i)); %


end

%% Plot Mechanism

for i=1:length(theta_2)
    %///////////////////////////////////////////////////////////////////////////
    %% Draw Points:
    if animate
        plot([-30,30],[0,0], clr, grnd_clr);                                % Clear graph and draw ground
        hold on;                                                            % redraw the current figure
        axis( [xmin xmax ymin ymax] );                                      % figure axis limits
        xlabel('x (cm)', 'fontsize', 15);                                   % axis label
        ylabel('y (cm)', 'fontsize', 15);                                   % axil label
        grid off;                                                           % add a grid to the figure
        axis equal;                                                         % make sure the figure is not stretched
        title('Plot of Patio Umbrella Mechanism');                          % add a title to the figure

        plot([-30,-20],[-10 0], clr, grnd_clr);                             % draw base pivot
        plot([-20 -10],[-10 0], clr, grnd_clr);                             % draw base pivot    
        plot([-10 000],[-10 0], clr, grnd_clr);                             % draw base pivot
        plot([000 010],[-10,0], clr, grnd_clr);                             % draw base pivot    
        plot([010 020],[-10 0], clr, grnd_clr);                             % draw base pivot
        plot([020 030],[-10,0], clr, grnd_clr);                             % draw base pivot



    %     plot(Ax(i), Ay(i), dot_clr);                                      % draw a point at A
    %     plot(Bx(i), By(i), dot_clr);                                      % draw a point at B
    %     plot(Cx(i), Cy(i), dot_clr);                                      % draw a point at C
    %     plot(Dx(i), Dy(i), dot_clr);                                      % draw a point at D
    %     plot(Ex(i), Ey(i), dot_clr);                                      % draw a point at E
    %     plot(Fx(i), Fy(i), dot_clr);                                      % draw a point at F
    %     plot(Gx(i), Gy(i), dot_clr);                                      % draw a point at G
    %     plot(Hx(i), Hy(i), dot_clr);                                      % draw a point at H
    %     plot(Ix(i), Iy(i), dot_clr);                                      % draw a point at I
    %     plot(Jx(i), Jy(i), dot_clr);                                      % draw a point at J
    %     plot(Kx(i), Ky(i), dot_clr);                                      % draw a point at K
    %     plot(Lx(i), Ly(i), dot_clr);                                      % draw a point at L
    %     plot(Mx(i), My(i), dot_clr);                                      % draw a point at M
    %     plot(Nx(i), Ny(i), dot_clr);                                      % draw a point at N
    %     plot(Ox(i), Oy(i), dot_clr);                                      % draw a point at O
    %     plot(Px(i), Py(i), dot_clr);                                      % draw a point at P
    %     plot(Qx(i), Qy(i), dot_clr);                                      % draw a point at Q
    %     plot(Rx(i), Ry(i), dot_clr);                                      % draw a point at R
    %     plot(Sx(i), Sy(i), dot_clr);                                      % draw a point at S
    %     plot(Tx(i), Ty(i), dot_clr);                                      % draw a point at T
    %     plot(Ux(i), Uy(i), dot_clr);                                      % draw a point at U
    %     plot(Vx(i), Vy(i), dot_clr);                                      % draw a point at V
    %     plot(Wx(i), Wy(i), dot_clr);                                      % draw a point at W


        %% Draw Links:
        plot([Ax(i) Cx(i)], [Ay(i) Cy(i)], clr, clr1, LW, 1 );              % Link2
        plot([Ox(i) Bx(i)], [Oy(i) By(i)], clr, clr1, LW, 1 );              % Link1
        plot([Dx(i) Ex(i)], [Dy(i) Ey(i)], clr, clr1, LW, 1 );              % Link3
        plot([Gx(i) Fx(i)], [Gy(i) Fy(i)], clr, clr1, LW, 1 );              % Link4
        plot([Hx(i) Ix(i)], [Hy(i) Iy(i)], clr, clr1, LW, 1 );              % Link5
        plot([Ix(i) Ex(i)], [Iy(i) Ey(i)], clr, clr1, LW, 1 );              % Link6
        plot([Ex(i) Jx(i)], [Ey(i) Jy(i)], clr, clr1, LW, 1 );              % Link6
        plot([Ix(i) Jx(i)], [Iy(i) Jy(i)], clr, clr1, LW, 1 );              % Link6
        plot([Jx(i) Lx(i)], [Jy(i) Ly(i)], clr, clr1, LW, 1 );              % Link7
        plot([Ex(i) Kx(i)], [Ey(i) Ky(i)], clr, clr1, LW, 1 );              % Link8
        plot([Lx(i) Nx(i)], [Ly(i) Ny(i)], clr, clr1, LW, 1 );              % Link9    
        plot([Nx(i) Px(i)], [Ny(i) Py(i)], clr, clr1, LW, 1 );              % Link10
        plot([Px(i) Qx(i)], [Py(i) Qy(i)], clr, clr1, LW, 1 );              % Link10
        plot([Qx(i) Rx(i)], [Qy(i) Ry(i)], clr, clr1, LW, 1 );              % Link10
        plot([Ux(i) Sx(i)], [Uy(i) Sy(i)], clr, clr4, LW, 1 );              % Link11
        plot([Tx(i) Sx(i)], [Ty(i) Sy(i)], clr, clr4, LW, 1 );              % Link12
        plot([Mx(i) Vx(i)], [My(i) Vy(i)], clr, clr3, LW, 1 );              % Link13
        plot([Mx(i) Wx(i)], [My(i) Wy(i)], clr, clr3, LW, 1 );              % Link14
        plot([Wx(i) Vx(i)], [Wy(i) Vy(i)], clr, clr3, LW, 1 );              % 15 Umbrella canopy bottom or outer ring
        
%         lengthL = num2str(sqrt((Mx(i)-Wx(i))^2+(My(i)-Wy(i))^2));
%         text(Mx(i)+0.4, Cy(i), lengthL , clr,'k');                        %Verify that the length stays constant
    end

    axis equal;                        % make sure the figure is not stretched
    axis( [xmin xmax ymin ymax] );          % figure axis limits
    pause(0.0); 
    
%     %% Compute Averages
    avg_001x(i) = (Ox(i) + Bx(i))/2;
    avg_001y(i) = (Oy(i) + By(i))/2;              % Link1
    avg_002x(i) = (Ax(i) + Cx(i))/2;
    avg_002y(i) = (Ay(i) + Cy(i))/2;              % Link2
    avg_003x(i) = (Dx(i) + Ex(i))/2;
    avg_003y(i) = (Dy(i) + Ey(i))/2;              % Link3
    avg_004x(i) = (Gx(i) + Fx(i))/2;
    avg_004y(i) = (Gy(i) + Fy(i))/2;              % Link4
    avg_005x(i) = (Hx(i) + Ix(i))/2;
    avg_005y(i) = (Hy(i) + Iy(i))/2;              % Link5
    avg_06ax(i) = (Ix(i) + Ex(i))/2;
    avg_06ay(i) = (Iy(i) + Ey(i))/2;              % Link6
    avg_06bx(i) = (Ex(i) + Jx(i))/2;
    avg_06by(i) = (Ey(i) + Jy(i))/2;              % Link6
    avg_06cx(i) = (Ix(i) + Jx(i))/2;
    avg_06cy(i) = (Iy(i) + Jy(i))/2;              % Link6
    avg_007x(i) = (Jx(i) + Lx(i))/2;
    avg_007y(i) = (Jy(i) + Ly(i))/2;              % Link7
    avg_008x(i) = (Ex(i) + Kx(i))/2;
    avg_008y(i) = (Ey(i) + Ky(i))/2;              % Link8
    avg_009x(i) = (Lx(i) + Nx(i))/2;
    avg_009y(i) = (Ly(i) + Ny(i))/2;              % Link9    
    avg_10ax(i) = (Nx(i) + Px(i))/2;
    avg_10ay(i) = (Ny(i) + Py(i))/2;              % Link10
    avg_10bx(i) = (Px(i) + Qx(i))/2;
    avg_10by(i) = (Py(i) + Qy(i))/2;              % Link10
    avg_10cx(i) = (Qx(i) + Rx(i))/2;
    avg_10cy(i) = (Qy(i) + Ry(i))/2;              % Link10
    avg_011x(i) = (Ux(i) + Sx(i))/2;
    avg_011y(i) = (Uy(i) + Sy(i))/2;              % Link11
    avg_012x(i) = (Tx(i) + Sx(i))/2;
    avg_012y(i) = (Ty(i) + Sy(i))/2;              % Link12
    avg_013x(i) = (Mx(i) + Vx(i))/2;
    avg_013y(i) = (My(i) + Vy(i))/2;              % Link13
    avg_014x(i) = (Mx(i) + Wx(i))/2;
    avg_014y(i) = (My(i) + Wy(i))/2;              % Link14
    avg_015x(i) = (Wx(i) + Vx(i))/2;
    avg_015y(i) = (Wy(i) + Vy(i))/2;              %15   
%     %hold on;                          % draw on top of the current figure
% 
%     %% Plot averages of each link just to show that they make sense
    if false
        plot(avg_001x(i), avg_001y(i), dot_clr);                                      % draw a point at E
        plot(avg_002x(i), avg_002y(i), dot_clr);                                      % draw a point at E
        plot(avg_003x(i), avg_003y(i), dot_clr);                                      % draw a point at E
        plot(avg_004x(i), avg_004y(i), dot_clr);                                      % draw a point at E
        plot(avg_005x(i), avg_005y(i), dot_clr);                                      % draw a point at E
        plot(avg_06ax(i), avg_06ay(i), dot_clr);                                      % draw a point at E
        plot(avg_06bx(i), avg_06by(i), dot_clr);                                      % draw a point at E
        plot(avg_06cx(i), avg_06cy(i), dot_clr);                                      % draw a point at E
        plot(avg_007x(i), avg_007y(i), dot_clr);                                      % draw a point at E
        plot(avg_008x(i), avg_008y(i), dot_clr);                                      % draw a point at E
        plot(avg_009x(i), avg_009y(i), dot_clr);                                      % draw a point at E
        plot(avg_10ax(i), avg_10ay(i), dot_clr);                                      % draw a point at E
        plot(avg_10bx(i), avg_10by(i), dot_clr);                                      % draw a point at E
        plot(avg_10cx(i), avg_10cy(i), dot_clr);   
        plot(avg_011x(i), avg_011y(i), dot_clr);   
        plot(avg_012x(i), avg_012y(i), dot_clr);   
        plot(avg_013x(i), avg_013y(i), dot_clr);   
        plot(avg_014x(i), avg_014y(i), dot_clr);   
        plot(avg_015x(i), avg_015y(i), dot_clr);   
    end
    
%     %% Find h_cg
    y_avg(i) = 0;
    y_avg(i) = y_avg(i) + w1*avg_001x(i);
    y_avg(i) = y_avg(i) + w2*avg_002x(i);
    y_avg(i) = y_avg(i) + w3*avg_003x(i);
    y_avg(i) = y_avg(i) + w4*avg_004x(i);
    y_avg(i) = y_avg(i) + w5*avg_005x(i);
    y_avg(i) = y_avg(i) + w6*avg_06ax(i);
    y_avg(i) = y_avg(i) + w6*avg_06bx(i);
    y_avg(i) = y_avg(i) + w6*avg_06cx(i);
    y_avg(i) = y_avg(i) + w7*avg_007x(i);
    y_avg(i) = y_avg(i) + w7*avg_008x(i);
    y_avg(i) = y_avg(i) + w9*avg_009x(i);
    y_avg(i) = y_avg(i) + w10A*avg_10ax(i);
    y_avg(i) = y_avg(i) + w10B*avg_10bx(i);
    y_avg(i) = y_avg(i) + w10C*avg_10cx(i);
    y_avg(i) = y_avg(i) + 4*w11*avg_011x(i);
    y_avg(i) = y_avg(i) + 4*w12*avg_012x(i);
    y_avg(i) = y_avg(i) + 4*w13*avg_013x(i);
    y_avg(i) = y_avg(i) + 4*w14*avg_014x(i);
    y_avg(i) = y_avg(i) / w_tot;
    
    C = (y_avg(i)-y_avg(i))/(Sy(i)-Sy(1));
    
    if animate                                  
%         plot(-100, y_avg, dot_clr);        % draw a point at y_Avg
        axis equal;                        % make sure the figure is not stretched
        axis([xmin xmax ymin ymax]);       % figure axis limits
        pause(0.0);
        hold off;
    end
end

%% Calculate the Force required for link 5 for range of canopy angles
F1 = zeros(1,length(theta_2));
for i = 1:length(theta_2)
    
    x_avg_canopy(i) = 0;        %CG of canopy
    x_avg_canopy(i) = x_avg_canopy(i) + w9*avg_009x(i);
    x_avg_canopy(i) = x_avg_canopy(i) + w10A*avg_10ax(i);
    x_avg_canopy(i) = x_avg_canopy(i) + w10B*avg_10bx(i);
    x_avg_canopy(i) = x_avg_canopy(i) + w10C*avg_10cx(i);
    x_avg_canopy(i) = x_avg_canopy(i) + 4*w11*avg_011x(i);
    x_avg_canopy(i) = x_avg_canopy(i) + 4*w12*avg_012x(i);
    x_avg_canopy(i) = x_avg_canopy(i) + 4*w13*avg_013x(i);
    x_avg_canopy(i) = x_avg_canopy(i) + 4*w14*avg_014x(i);
    x_avg_canopy(i) = x_avg_canopy(i) / (W_canopy+w10A+w10B+w10C);

    y_avg_canopy(i) = 0;        %CG of canopy
    y_avg_canopy(i) = y_avg_canopy(i) + w9*avg_009y(i);
    y_avg_canopy(i) = y_avg_canopy(i) + w10A*avg_10ay(i);
    y_avg_canopy(i) = y_avg_canopy(i) + w10B*avg_10by(i);
    y_avg_canopy(i) = y_avg_canopy(i) + w10C*avg_10cy(i);
    y_avg_canopy(i) = y_avg_canopy(i) + 4*w11*avg_011y(i);
    y_avg_canopy(i) = y_avg_canopy(i) + 4*w12*avg_012y(i);
    y_avg_canopy(i) = y_avg_canopy(i) + 4*w13*avg_013y(i);
    y_avg_canopy(i) = y_avg_canopy(i) + 4*w14*avg_014y(i);
    y_avg_canopy(i) = y_avg_canopy(i) / (W_canopy+w10A+w10B+w10C);


    %Two force analysis calculation of vertical force input on link 5
    
    F1(i) = (W_canopy+w10A+w10B+w10C)*sqrt((x_avg_canopy(i)-Kx(i))^2+(y_avg_canopy(i)-Ky(i))^2)*cos(canopy_angle(i))/(Fx(i)-Hx(i));
end

figure(3);                         %create new figure
set(3,'WindowStyle','Docked')      %dock the figure
plot(canopy_angle/deg-90,F1+w5);
title("Force vs Canopy angle");
xlabel("Canopy angle");
ylabel("Force including link 5 weight")

