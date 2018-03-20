% Single point positioning
% Wouter van der Wal
% Delft University of Technology

clear all
% close all
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
% Homework Assignment 3: Least square estimation of an initial state vector
% and a dynamic parameter
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

%% Part 1a
% Final values of the epochwise distance is stored in the variable
% s_final_rcc. Please refer to previous assignment or follow the github
% link. While the clock error is stored in del_t variable. Both of them are
% used to reconstruct the initial state vector.

% Loading the initial state vector from previous script 
pos =  load('s_final_rcc','-mat'); % Position in m
clock_errors = load('del_t','-mat'); % Clock errors in s

% s_final_rcc = pos.s_final_rcc(end,:); % Final x,y,z position for last epoch with clock correction implemented
s_final_rcc = pos.s_final_rcc(1,:); % Final x,y,z position for last epoch with clock correction implemented
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
% s0 = [ 0.1 0.1 0.1 0.1 0.1 0.1];
% Standard Gravitational parameter Earth
mu = 3.986004418e14; % [ m^3 s^-2]
% Loop initializations
diff4 = norm(s0); % Current error resolution
iter4 = 0; % Current iteration step
diff_iter = [];
Svect = [];
res = [];
while diff4>= 10^-3
% while iter4<= 10
    A4 = [];
    y_bar4 = [];
    iter4 = iter4 +1;
    n = 0;
    phi = eye(6,6); % Initializing state transition matrix phi
    S = s0; 
    Svect = [Svect; S];
    % For all epochs
    for e=1:epoch
        clear H dx4 dy4 dz4 drho4;
        r = norm(S(1,1:3)); % radial distance sqrt(x0, y0, z0) in meters
        % Calculating dF/dS matrix for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        dF_dS = [ 0 0 0 1 0 0; 0 0 0 0 1 0 ; 0 0 0 0 0 1 ;...
                  (mu/r^3)*( (3*S(1)^2)/(r^2) - 1),    3*mu*S(1)*S(2)/r^5,   3*mu*S(1)*S(3)/r^5, 0, 0, 0;...
                  3*mu*S(1)*S(2)/r^5, (mu/r^3)*( (3*S(2)^2)/(r^2) -1),     3*mu*S(2)*S(3)/r^5, 0, 0, 0   ;...
                  3*mu*S(1)*S(3)/r^5,  3*mu*S(2)*S(3)/r^5, (mu/r^3)*( (3*S(3)^2)/(r^2) -1), 0, 0, 0];
        
        % Calculating phi_dot
        phi_dot = dF_dS*phi;
        % Updating phi 
        phi = phi + phi_dot*dt;
        % Calculating S_dot for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        S_dot = [ S(4); S(5); S(6); -mu*S(1)/r^3; -mu*S(2)/r^3 ; -mu*S(3)/r^3 ]';
        % Updating S 
        S = S + S_dot.*dt;
        
        for s = 1:no_sat(e)   
            dx4(s) = S(1) - xgi(e,s); % x - x_gps for sth satellite in view
            dy4(s) = S(2) - ygi(e,s); % y - y_gps for sth satellite in view
            dz4(s) = S(3) - zgi(e,s); % z - z_gps for sth satellite in view

            drho4(s) = sqrt(dx4(s)^2 + dy4(s)^2 + dz4(s)^2); % measured range for nth satellite in view
            H(s,:) = [ dx4(s)/drho4(s) dy4(s)/drho4(s) dz4(s)/drho4(s) 0 0 0]; % Assembling epochwise information matrix.
        end
        B = H*phi; % New information matrix.
        % Accumalating primary information matrix from all the epochs into
        % one matrix
%         A4 = blkdiag(A4,B);
        A4 = vertcat(A4,B);
        % Accumulating observation for all epochs
        % y_bar4 = [y_bar4;transpose(C1(e,1:no_sat(e)) + ccGPS(e,1:no_sat(e)) - drho_rcc) ] ;
        % With inclusion of receiver clock correction from previous assignment
        y_bar4 = [y_bar4;transpose(C1(e,1:no_sat(e)) + (ccGPS(e,1:no_sat(e))- final_clock_error*c) - drho4) ] ;
        n= n+3;
    end
    
    % rank check to check for information loss in the information matrix
    rank4(iter4) = rank(A4);
    
    % LSQ with the aid of SVD
    [U4,S4,V4] = svd(A4,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
    S_inv_diag4 = [];
    for nn=1:length(S4)
        S_inv_diag4 = [S_inv_diag4,inv(S4(nn,nn))];
    end
    S_inv4 = diag(S_inv_diag4);
    
    % Estimation of small letter case, s. s = x_hat4
    x_hat4 = V4*S_inv4*U4'*y_bar4; % Parameter estimation via SVD aided LSQ
    y_hat4 = A4*x_hat4; % Measured vector 
    
    res = [res, norm(y_bar4-y_hat4)];
    
    nn = 0;
    buf = [];
    for ee = 1:epoch
        buf = [buf,nn+1,sum(no_sat(1:ee))];
        res4 = y_bar4(nn+1:sum(no_sat(1:ee))) - y_hat4(nn+1:sum(no_sat(1:ee))); % residual in meters
        res_norm4(ee,iter4) = norm(res4);% Epoch wise norm of the residual vector of all satellites               
        nn = sum(no_sat(1:ee));
    end
    % Update error and state
    diff4 = norm(x_hat4(1:3));
    diff_iter = [diff_iter;diff4];
    s0 =s0 + x_hat4';
    if iter4>=100
        disp('Iteration limit of 100 exceeds');
        break
    end
end

%% Plotting the behaviour of epochwise residual, inclusion of clock correction before batch LSQ for all epochs
for ee = 1:epoch
        figure(31)
        subplot(6,2,ee)
        plot(1:length(res_norm4(ee,:)),res_norm4(ee,:),'rx');
        hold on 
        plot(1:length(res_norm4(ee,:)),res_norm4(ee,:),'k:');
        hold off
        title(['Epoch ' num2str(ee)]);
        ylabel('Res. norm [m]');
        xlabel('Iteration no.');
        grid on
end

nn = 0;
for ee = 1:epoch
        yb =  y_bar4(nn+1:sum(no_sat(1:ee)),1)';
        yh = y_hat4(nn+1:sum(no_sat(1:ee)),1)';
        nn = sum(no_sat(1:ee));
        figure(32) 
        subplot(6,2,ee)
        plot(1:length(yb),yb,'r-x');
        hold on 
        plot(1:length(yh),yh,'b-*');
        hold off  
        title(['Epoch ' num2str(ee)]);
        ylabel('Range');
        xlabel('Sat no.');
        grid on
end
%% Part b 
% Final values of the epochwise distance is stored in the variable
% s_final_rcc. Please refer to previous assignment or follow the github
% link. While the clock error is stored in del_t variable. Both of them are
% used to reconstruct the initial state vector.

% Loading the initial state vector from previous script 
pos =  load('s_final_rcc','-mat'); % Position in m
clock_errors = load('del_t','-mat'); % Clock errors in s

% s_final_rcc = pos.s_final_rcc(end,:); % Final x,y,z position for last epoch with clock correction implemented
s_final_rcc = pos.s_final_rcc(1,:); % Final x,y,z position for last epoch with clock correction implemented
final_clock_error = clock_errors.del_t(1,end); % Final receiver clock correction

% Basic backward Euler scheme to evaluate velocity terms in m/s
% vel = (pos.s_final_rcc(end,:) - pos.s_final_rcc(end-1,:))/dt;
vel = (pos.s_final_rcc(2,:) - pos.s_final_rcc(1,:))/dt;
% Intializing state vector: [ x, y, z, x_dot, y_dot, z_dot] with clock
% correction implemented
s0_buf = [ s_final_rcc, vel, 3e14] + [ 0 0 0 0 0 0 0 ];
% Shifting initial vector slightly to visualize convergence.
s0 = s0_buf ;
% s0 = [ 0.1 0.1 0.1 0.1 0.1 0.1];
% Standard Gravitational parameter Earth
% mu = 3.986004418e14; % [ m^3 s^-2]
% Loop initializations
diff4 = norm(s0); % Current error resolution
iter4 = 0; % Current iteration step
diff_iter = [];
Svect = [];
res = [];
while diff4>= 10^-3
% while iter4<= 10
    A4 = [];
    y_bar4 = [];
    iter4 = iter4 +1;
    phi = eye(7,7); % Initializing state transition matrix phi
    S = s0; 
    Svect = [Svect; S];
    % For all epochs 
    mu = S(7);
    for e=1:epoch
        clear H dx4 dy4 dz4 drho4;
        r = norm(S(1,1:3)); % radial distance sqrt(x0, y0, z0) in meters
        % Calculating dF/dS matrix for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        dF_dS = [ 0 0 0 1 0 0 0; 0 0 0 0 1 0 0 ; 0 0 0 0 0 1 0;...
                  (mu/r^3)*( (3*S(1)^2)/(r^2) - 1),    3*mu*S(1)*S(2)/r^5,   3*mu*S(1)*S(3)/r^5, 0, 0, 0, -S(1)/r^3;...
                  3*mu*S(1)*S(2)/r^5, (mu/r^3)*( (3*S(2)^2)/(r^2) -1),     3*mu*S(2)*S(3)/r^5, 0, 0, 0, -S(2)/r^3   ;...
                  3*mu*S(1)*S(3)/r^5,  3*mu*S(2)*S(3)/r^5, (mu/r^3)*( (3*S(3)^2)/(r^2) -1), 0, 0, 0, -S(3)/r^3; ...
                  0, 0, 0, 0, 0, 0, 0];
        
        % Calculating phi_dot
        phi_dot = dF_dS*phi;
        % Updating phi 
        phi = phi + phi_dot*dt;
        % Calculating S_dot for s0 = [ x, y, z, x_dot, y_dot, z_dot ]
        S_dot = [ S(4); S(5); S(6); -mu*S(1)/r^3; -mu*S(2)/r^3 ; -mu*S(3)/r^3; 0 ]';
        % Updating S 
        S = S + S_dot.*dt;
        
        for s = 1:no_sat(e)   
            dx4(s) = S(1) - xgi(e,s); % x - x_gps for sth satellite in view
            dy4(s) = S(2) - ygi(e,s); % y - y_gps for sth satellite in view
            dz4(s) = S(3) - zgi(e,s); % z - z_gps for sth satellite in view

            drho4(s) = sqrt(dx4(s)^2 + dy4(s)^2 + dz4(s)^2); % measured range for nth satellite in view
            H(s,:) = [ dx4(s)/drho4(s) dy4(s)/drho4(s) dz4(s)/drho4(s) 0 0 0 0]; % Assembling epochwise information matrix.
        end
        B = H*phi; % New information matrix.
        % Accumalating primary information matrix from all the epochs into
        % one matrix
        A4 = vertcat(A4,B);
        % Accumulating observation for all epochs
        % y_bar4 = [y_bar4;transpose(C1(e,1:no_sat(e)) + ccGPS(e,1:no_sat(e)) - drho_rcc) ] ;
        % With inclusion of receiver clock correction from previous assignment
        y_bar4 = [y_bar4;transpose(C1(e,1:no_sat(e)) + (ccGPS(e,1:no_sat(e))- final_clock_error*c) - drho4) ] ;
    end
    
    % rank check to check for information loss in the information matrix
    rank4(iter4) = rank(A4);
    
    % LSQ with the aid of SVD
    [U4,S4,V4] = svd(A4,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
    S_inv_diag4 = [];
    for nn=1:length(S4)
        S_inv_diag4 = [S_inv_diag4,inv(S4(nn,nn))];
    end
    S_inv4 = diag(S_inv_diag4);
    
    % Estimation of small letter case, s. s = x_hat4
    x_hat4 = V4*S_inv4*U4'*y_bar4; % Parameter estimation via SVD aided LSQ
    y_hat4 = A4*x_hat4; % Measured vector 
    
    res = [res, norm(y_bar4-y_hat4)];
    
    nn = 0;
    buf = [];
    for ee = 1:epoch
        buf = [buf,nn+1,sum(no_sat(1:ee))];
        res4 = y_bar4(nn+1:sum(no_sat(1:ee))) - y_hat4(nn+1:sum(no_sat(1:ee))); % residual in meters
        res_norm4(ee,iter4) = norm(res4);% Epoch wise norm of the residual vector of all satellites               
        nn = sum(no_sat(1:ee));
    end
    % Update error and state
    diff4 = norm(x_hat4(1:3));
    diff_iter = [diff_iter;diff4];
    s0 =s0 + x_hat4';
    if iter4>=100
        disp('Iteration limit of 100 exceeds');
        break
    end
end

%% Plotting the behaviour of epochwise residual, inclusion of clock correction before batch LSQ for all epochs
for ee = 1:epoch
        figure(33) 
        subplot(6,2,ee)
        plot(1:length(res_norm4(ee,:)),res_norm4(ee,:),'rx');
        hold on 
        plot(1:length(res_norm4(ee,:)),res_norm4(ee,:),'k:');
        hold off
        title(['Epoch' num2str(ee)]);
        ylabel('Res. norm [m]');
        xlabel('Iteration no.');
        grid on
end

nn = 0;
for ee = 1:epoch
        yb =  y_bar4(nn+1:sum(no_sat(1:ee)),1)';
        yh = y_hat4(nn+1:sum(no_sat(1:ee)),1)';
        nn = sum(no_sat(1:ee));
        figure(34) 
        subplot(6,2,ee)
        plot(1:length(yb),yb,'r-x');
        hold on 
        plot(1:length(yh),yh,'b-*');
        hold off  
        title(['Epoch ' num2str(ee)]);
        ylabel('Range');
        xlabel('Sat no.');
        grid on
end

%% Part C Standard deviation estimation.
% Assumed standard deviation of observation vector
std_y = 5;
% Observation vector co-variance matrix
Pyy = (std_y^2)*eye(87,87);
% State co-variance matrix
Pxx = pinv(A4)*Pyy*pinv(A4');
% Standard deviation of gravitational parameter
std_mu = sqrt(Pxx(end,end));
% mu combined with its standard deviation
mu_with_std = s0(7)+std_mu;
% Percent difference between estimated and actual mu.
mu_percent_diff = ((3.986004419*10^(14) - mu_with_std)/(3.986004419*10^(14)) )*100;
