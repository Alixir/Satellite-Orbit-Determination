% Single point positioning
% Wouter van der Wal
% Delft University of Technology

clear all
close all
clc

% Fundamental GPS constants
gpsconst;

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
% Homework Assignment 2: GPS single point positioning
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

%% Part 1a: Epoch wise position estimation of SWARM satellite
rank_e = []; % initialising rank check per epoch for the 
%Scan through each epoch
for e = 1:epoch
    clear A dx dy dz drho; 
    s0 = [ 10^3; 10^3; 10^3 ]; % initial guess x0, y0, z0 in meters
    diff = norm(s0); % Initialise error limit
    iter =0; % initialising iteration number
    while diff>=10^(-15)
        iter = iter +1; % current iteration number
        for s = 1:no_sat(e)
            dx(s) = s0(1) - xgi(e,s); % x - x_gps for sth satellite in view
            dy(s) = s0(2) - ygi(e,s); % y - y_gps for sth satellite in view
            dz(s) = s0(3) - zgi(e,s); % z - z_gps for sth satellite in view
            
            drho(s) = sqrt(dx(s)^2 + dy(s)^2 + dz(s)^2); % measured range for nth satellite in view
            A(s,:) = [ dx(s)/drho(s) dy(s)/drho(s) dz(s)/drho(s)]; % Assembling epochwise information matrix. Final matrix for each epoch= no. of satellites in view/ epoch by 3
        end 
        % Determining observation vector y_bar
        % y_bar = rho(s) - rho(s0) + c*t_T
        % y_bar = observed pseudorange - measured range for given initial
        % guess + effect of GPS clock correction
        y_bar = transpose(C1(e,1:no_sat(e)) - drho + ccf(e,1:no_sat(e)));
        rank_e(e) = rank(A); % Rank check before application of LSQ with pinv
        % Unwieghted LSQ to obtain the state vector x_hat
        % LSQ with pseudo-inverse
%         x_hat = (A'*A)\A'*y_bar; 
        % LSQ with SVD
        [U,S,V] = svd(A,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
%         condition_number = S(1,1)/S(end,end); % Indicates the number of potential digit that can be lost in any numerical calculations when A is computed directly.
        S_inv_diag = [];
        for nn=1:length(S)
            S_inv_diag = [S_inv_diag,inv(S(nn,nn))];
        end
        S_inv = diag(S_inv_diag);
        x_hat = V*S_inv*U'*y_bar; % Parameter estimation via SVD aided LSQ
        
        % Estimated vector
        y_hat = A*x_hat; 
        
        res = y_bar - y_hat; % residual in meters
        res_norm(e,iter) = norm(res);% Epoch wise norm of the residual vector of all satellites
        
        % Updating error limit and state
        diff = abs( ( norm(x_hat + s0) - norm(s0) )/ ( norm(s0)));
        s0 = s0 + x_hat;
        % Setting a limit for iteration
        if iter>=1000
            disp(['Iteration exceeds for epoch:' num2str(e)]);
        end
        s_final(e,:) = transpose(s0); % Final value of s = [x;y;z] per epoch
        r_measured(e,:) = (norm(s_final(e,:)) * 10^(-3)); % Measured distance in km.
        dist_measured(e,:) = r_measured(e,:) - earth_radius; % Measured distance in radial direction in km.
    end    
    % Plotting epoch wise residual norms  
            figure(11)
            subplot(6,2,e);
            plot(1:length(res_norm(e,:)),res_norm(e,:),'r*');
            xlabel('Iteration Number');
            ylabel('Residual Norm [m]');
            grid on
            title(strcat('Epoch',num2str(e)));    
end

%% Part 1b: Difference between measured solution and precise orbit direction
% Plotting measured and observed radial distances in km
dist_precise = rp - earth_radius;
figure(21)
plot(1:epoch,dist_precise,'r-',1:epoch,dist_measured,'b--');
legend('Precise radial distance','Measured radial distance','Location','Best');
xlabel('Epoch');
ylabel('Altitude [km]');
grid on;
% Plotting the difference between measured and observed orbit
figure(22)
plot(1:epoch,dist_measured - dist_precise');
xlabel('Epoch');
ylabel('Difference between measured and precise radial orbit [km]');
grid on;

