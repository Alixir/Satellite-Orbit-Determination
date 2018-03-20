% Single point positioning
% Wouter van der Wal
% Delft University of Technology

clear all
close all
clc

% Fundamental GPS constants
% gpsconst;

dt                  = 10;       % time step
number_of_epochs    = 11;       
order               = 11;       % order of interpolation of precise orbits GPS satellites
order_clk           = 3;        % order of interpolation of clock correction
t_orbit             = 60*[-78:15:73];  % time vector for precise orbits for GPS satellite provided by IGS every 15 minutes
time                = 0:dt:(number_of_epochs-1)*dt; % time vector for SWARM satellite
c                   = 299792458; % velocity of light m/sec
earth_radius        = 6371; % in km 
%% Read the IGS coordinates and clock corrections for all the GPSsatellites
[x_gps,y_gps,z_gps,clock_correction] = readsp3('igs19301.sp3');


%% Read the Swarm GPS data
fid_rinex = fopen('Data_SWARM_A.dat','r');
for ii = 1:15
    line = fgets(fid_rinex);
end
% Read the observations
for epoch = 1:number_of_epochs
    line = fgets(fid_rinex);
    ft = sscanf(line,'%c %f %f %f %f %f %f %f %f');
    no_sat(epoch) = ft(9);
    si = 0;
    for no = 1:no_sat(epoch)
        si = si + 1;
        line = fgets(fid_rinex);
        ft = sscanf(line,'%c %f %f %f %f %f %f %f %f %f');
        prns(epoch,no) = ft(2);        % PRN number of satellite
        n = prns(epoch,no);
        C1(epoch,si) = ft(5);
        % Interpolate the coordinates of the GPS satellites to the current epoch
        [x_gps_interpolated(epoch,si),x_error(epoch,si)] = polint(t_orbit,x_gps(:,n),time(epoch),order);
        [y_gps_interpolated(epoch,si),y_error(epoch,si)] = polint(t_orbit,y_gps(:,n),time(epoch),order);
        [z_gps_interpolated(epoch,si),z_error(epoch,si)] = polint(t_orbit,z_gps(:,n),time(epoch),order);
        [clock_correction_interpolated(epoch,si),t_error(epoch,si)] = polint(t_orbit,clock_correction(:,n),time(epoch),order_clk);
    end
end
fclose(fid_rinex);
ccGPS = clock_correction_interpolated*c;
x_gps_interpolated = x_gps_interpolated*1e3;                      % Convert to meters
y_gps_interpolated = y_gps_interpolated*1e3;
z_gps_interpolated = z_gps_interpolated*1e3;
r_gps_interpolated = sqrt(x_gps_interpolated(:,1).^2+y_gps_interpolated(:,1).^2+z_gps_interpolated(:,1).^2);
% Plot for checking
% figure; plot((r_gps_interpolated-6371e3)/1e3); title('altitude PRN 7 [km]')


%% Read the precise Swarm orbit
fid_obs = fopen('PreciseOrbit_SWARM_A.dat','r');
for ii = 1:22
    line = fgets(fid_obs);
end
for epoch = 1:11
        line = fgets(fid_obs);
        line = fgets(fid_obs);
        x_precise(epoch) = read_strval(line(5:18)); % These are in km from the center of the earth - Ali
        y_precise(epoch) = read_strval(line(19:32));
        z_precise(epoch) = read_strval(line(33:46));
        line = fgets(fid_obs);
end
fclose(fid_obs);
% Convert to spherical coordinates for checking
r_precise = sqrt( x_precise.^2+y_precise.^2+z_precise.^2 );
% Plot for checking
% figure; plot(r_precise-6371); title('altitude Swarm [km]')
%% ========== START OF STUDENT SCRIPT ==============
%% Author: Ali Nawaz
% Course: AE4872 Satellite Orbit Determination
% Homework Assignment 4: Extended Kalman Filter
%%===================================================
%% Data retrieval
xg = x_gps;
yg = y_gps;
zg = z_gps;
cc = clock_correction;

xgi = x_gps_interpolated;
ygi = y_gps_interpolated;
zgi = z_gps_interpolated;
rgi = r_gps_interpolated;
cci = clock_correction_interpolated;
ccf = ccGPS; % Final clock correction = cci*c

xp = x_precise;
yp = y_precise;
zp = z_precise;
rp = r_precise;

%% Task A - Data Initialisation from Assignment 2
% Final values of the epochwise distance is stored in the variable
% s_final_rcc. Please refer to previous assignment or follow the github
% link. While the clock error is stored in del_t variable. Both of them are
% used to reconstruct the initial state vector.

% Loading the initial state vector from previous script 
pos =  load('s_final_rcc','-mat'); % Position in m
clock_errors = load('del_t','-mat'); % Clock errors in s

% s_final_rcc = pos.s_final_rcc(end,:); % Final x,y,z position for last epoch with clock correction implemented
s_final_rcc = pos.s_final_rcc(1,:); % Final x,y,z position for first epoch with clock correction implemented
final_clock_error = clock_errors.del_t(1,end); % Final receiver clock correction

% Basic backward Euler scheme to evaluate velocity terms in m/s
% vel = (pos.s_final_rcc(end,:) - pos.s_final_rcc(end-1,:))/dt;
vel = (pos.s_final_rcc(2,:) - pos.s_final_rcc(1,:))/dt;
% Intializing state vector: [ x, y, z, x_dot, y_dot, z_dot] with clock
% correction implemented
s0_buf = [ s_final_rcc, vel] +  [ 0 0 0 0 0 0 ];

% Shifting initial vector slightly to visualize convergence.
% s0 = s0_buf + [ 1000 1000 1000 100 100 100];
s0 = s0_buf;

%% Task A

% Sequentially scan through all Epochs
S = s0_buf';
mu = 3.986004419e14; % [ m^3 s^-2]
% Minimum standard deviation in the state vector, using the known precise
% values for first epoch.
init_diff = s0_buf - [xp(1) yp(1) zp(1) (xp(2)-xp(1))/dt (yp(2)-yp(1))/dt (zp(2)-zp(1))/dt];
P = diag([init_diff.^2]);
% Effect of changing weights on observation vs 
% P = P + diag([ 10^20 10^20 10^20 10^10  10^10 10^10 ]);
% P = diag([ 0 0 0 0 0 0 ]);
dist = [];
phi = eye(6,6);
% Lambda parameter choices
% lambda_pos = 0.001;
% lambda_vel = 0.1;
% lambda_pos = 0;
% lambda_vel = 0;
% lambda_pos = mean([(0.4366-0.3549)^2,(0.3549-0.1093)^2]);
% lambda_vel = lambda_pos/(dt^2);
% lambda_pos = mean([(0.4366-0.3549),(0.3549-0.1093)]);
% lambda_vel = lambda_pos/(dt);
% lambda_pos = mean([((0.4366-0.3549)*10^3)^2,((0.3549-0.1093)*10^3)^2]);
% lambda_vel = lambda_pos/(dt)^2;
% lambda_pos = mean([((0.4366-0.3549)*10^3),((0.3549-0.1093)*10^3)]);
% lambda_vel = lambda_pos/(dt);
% lambda_pos = ((436.6+354.9)/2)^2;
lambda_pos = (400)^2;
lambda_vel = (436.6-354.9)^2/dt^2;
% Noise matrix added to state covariance matrix P to avoid filter
% saturation.
Q = diag([lambda_pos,lambda_pos,lambda_pos,lambda_vel,lambda_vel,lambda_vel]);
% Error co-variance parameters
lvl = [lambda_vel];
lpl = [lambda_pos];
% Collect estimated radial distances
S_collect = [norm(S)];
for e = 1:epoch
        clear H H_tilde dx4 dy4 dz4 drho4;
        r = norm(S(1:3)); % radial distance sqrt(x0, y0, z0) in meters
        % Calculating dF/dS matrix for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        dF_dS = [ 0 0 0 1 0 0; 0 0 0 0 1 0 ; 0 0 0 0 0 1 ;...
                  (mu/r^3)*( (3*S(1)^2)/(r^2) - 1),    3*mu*S(1)*S(2)/r^5,   3*mu*S(1)*S(3)/r^5, 0, 0, 0;...
                  3*mu*S(1)*S(2)/r^5, (mu/r^3)*( (3*S(2)^2)/(r^2) -1),     3*mu*S(2)*S(3)/r^5, 0, 0, 0   ;...
                  3*mu*S(1)*S(3)/r^5,  3*mu*S(2)*S(3)/r^5, (mu/r^3)*( (3*S(3)^2)/(r^2) -1), 0, 0, 0];
        
        % Calculating phi_dot
        phi_dot = dF_dS*phi;
        % Updating phi 
        phi = phi + phi_dot*dt;
        
        % Propagating state co-variance matrix P
        P = phi*P*phi' + Q ;
        
        % Calculating S_dot for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        S_dot = [ S(4); S(5); S(6); -mu*S(1)/r^3; -mu*S(2)/r^3 ; -mu*S(3)/r^3 ];
        
% Propagate Sj to Sk
% Euler integration
%         S = S + S_dot.*dt;
% State transition matrix
          S = phi*S;
% ode45
%         options = odeset('AbsTol',1e-10, 'RelTol',1e-10); 
%         [T,S] = ode45(@Sdot, [time(e)-dt, time(e)],S,options);
%         S = S(end,:)';

% Design matrix generation
        for s = 1:no_sat(e)   
            dx4(s) = S(1) - xgi(e,s); % x - x_gps for sth satellite in view
            dy4(s) = S(2) - ygi(e,s); % y - y_gps for sth satellite in view
            dz4(s) = S(3) - zgi(e,s); % z - z_gps for sth satellite in view

            drho4(s) = sqrt(dx4(s)^2 + dy4(s)^2 + dz4(s)^2); % measured range for nth satellite in view
            H(s,:) = [ dx4(s)/drho4(s) dy4(s)/drho4(s) dz4(s)/drho4(s) 0 0 0]; % Assembling epochwise information matrix.
        end
        H_tilde = H*phi;
        y_bar = transpose(C1(e,1:no_sat(e)) + (ccGPS(e,1:no_sat(e))- final_clock_error*c) - drho4);
        
        % Rk, Error covariance matrix R at kth epoch
        R = 5^2 * eye(size(y_bar));
        % Kalman Gain K for epoch k
        K = P*H'/(H*P*H' +R);
        % Filter update for x_hat
        x_hat = K*(y_bar);
        % Update S for current epoch
        S = S + x_hat;
        % Update P for current epoch
        P = (eye(size(K*H)) - K*H)*P;
        dist = [dist, norm(S(1:3))]; % Distance from the center of the ref. frame in m
        
        % Iterative Q explained in the report.
        S_collect = [S_collect,norm(S(1:3))];
        lpl = [lpl,(S_collect(e)-r_precise(e)*10^3)^2];
        lvl = [lvl, ((S_collect(end)-S_collect(end-1))^2)/dt^2];
        lambda_pos = lpl(end);
        lambda_vel = lvl(end);
        Q = diag([lambda_pos,lambda_pos,lambda_pos,lambda_vel,lambda_vel,lambda_vel]);
end
radial_dist = dist*10^(-3); % Altitude in km
dist_a2_data =  load('dist_a2_rcc','-mat'); % Position in m for batch with receiver clock correction
dist_a2 = dist_a2_data.dist_measured_rcc;
dist_a2seq_data =  load('dist_a2_seq','-mat'); % Position in km for sequential with receiver clock correction
dist_a2_seq = dist_a2seq_data.dist_measured;

figure(1)
plot(1:length(radial_dist),radial_dist-earth_radius,'r-*',1:length(r_precise),r_precise-earth_radius,'b-o',1:length(dist_a2),dist_a2,'k--s',1:length(dist_a2_seq),dist_a2_seq,'b-+');
grid on
legend('Estimated with EKF','Precise','Assignment 2 batch with dyn. paramter','Assignment 2 epochwise','Location','Best');
xlabel('Epoch');
ylabel('Altitude [km]');
title('SWARM A position comparision');

% Altitudes
diff_radial = radial_dist-r_precise; % Difference in km
diff_radial_a2 = dist_a2 - (r_precise-6371)';
diff_radial_a2seq = dist_a2_seq - (r_precise-6371)';

figure(2)
plot(1:11,zeros(1,11),'b-o',1:length(diff_radial),diff_radial,'r-*',1:length(diff_radial_a2),diff_radial_a2,'k--s',1:length(diff_radial_a2seq),diff_radial_a2seq,'b-+');
title('Difference between estimated and precise SWARM position');
legend('Precise','Estimated with EKF','Assignment 2 batch with dyn. paramter','Assignment 2 epochwise','Location','Best');
xlabel('Epoch');
ylabel('Position [km]');
grid on

% Sum of squared error residual:
sser = diff_radial*diff_radial';