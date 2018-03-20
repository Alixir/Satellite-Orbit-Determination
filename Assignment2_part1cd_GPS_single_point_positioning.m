%% Author: Ali Nawaz
% Course: AE4872 Satellite Orbit Determination
% Homework Assignment 2: GPS single point positioning
%%===================================================
%% Part 1C,D
% Initialize state vector and error limit
s0_rcc = (10^3)*ones(3*epoch +1,1);
diff_rcc = norm(s0_rcc);
iter_rcc = 0;
del_t = 0;
while diff_rcc>=10^-15
    A_rcc =[];
    y_bar_rcc = [];
    iter_rcc = iter_rcc +1;
    n= 0;
    for e=1:epoch
        clear Ac dx_rcc dy_rcc dz_rcc drho_rcc;
        for s = 1:no_sat(e)
            dx_rcc(s) = s0_rcc(n+1) - xgi(e,s); % x - x_gps for sth satellite in view
            dy_rcc(s) = s0_rcc(n+2) - ygi(e,s); % y - y_gps for sth satellite in view
            dz_rcc(s) = s0_rcc(n+3) - zgi(e,s); % z - z_gps for sth satellite in view

            drho_rcc(s) = sqrt(dx_rcc(s)^2 + dy_rcc(s)^2 + dz_rcc(s)^2); % measured range for nth satellite in view
            Ac(s,:) = [ dx_rcc(s)/drho_rcc(s) dy_rcc(s)/drho_rcc(s) dz_rcc(s)/drho_rcc(s)]; % Assembling epochwise information matrix. Final matrix for each epoch= no. of satellites in view/ epoch by 3
        end
        % Accumalating primary information matrix from all the epochs into
        % one matrix
        A_rcc = blkdiag(A_rcc,Ac);
        % Accumulating observation for all epochs
%         y_bar_rcc = [y_bar_rcc;transpose(C1(e,1:no_sat(e)) + ccGPS(e,1:no_sat(e)) - drho_rcc) ] ;
        % Inclusion of light time effect
        y_bar_rcc = [y_bar_rcc;transpose(C1(e,1:no_sat(e)) + (ccGPS(e,1:no_sat(e))- del_t(end)*c) - drho_rcc) ] ;
        n= n+3;
    end
    % Appending clock correction part to the all-epoch information matrix.
    mat_height = size(A_rcc);
    A_rcc = [A_rcc,ones(mat_height(1),1)];
    
    % rank check to check for information loss in the information matrix
    rank_rcc(iter_rcc) = rank(A_rcc);
    
    % LSQ with the aid of SVD
    [Urcc,Srcc,Vrcc] = svd(A_rcc,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
    S_inv_diag_rcc = [];
    for nn=1:length(Srcc)
        S_inv_diag_rcc = [S_inv_diag_rcc,inv(Srcc(nn,nn))];
    end
    S_inv_rcc = diag(S_inv_diag_rcc);
    
    x_hat_rcc = Vrcc*S_inv_rcc*Urcc'*y_bar_rcc; % Parameter estimation via SVD aided LSQ
    y_hat_rcc = A_rcc*x_hat_rcc; % Measured vector

    nn = 0;
    for ee = 1:epoch
        res_rcc = y_bar_rcc(nn+1:no_sat(ee)) - y_hat_rcc(nn+1:no_sat(ee)); % residual in meters
        res_norm_rcc(ee,iter_rcc) = norm(res_rcc);% Epoch wise norm of the residual vector of all satellites               
        nn = no_sat(ee);
    end
    
    % Update error and state
    diff_rcc = abs( ( norm(x_hat_rcc + s0_rcc) - norm(s0_rcc) )/ ( norm(s0_rcc + x_hat_rcc)));
    s0_rcc = s0_rcc + x_hat_rcc;
    % Clock correction per iteration
    del_t(iter_rcc) = s0_rcc(end)/c;
    if iter_rcc>=100
        disp('Iteration limit of 100 exceeds');
        break
    end
end
%% Plotting the behaviour of epochwise residual, inclusion of clock correction before batch LSQ for all epochs
for ee = 1:epoch
        figure(31) 
        subplot(6,2,ee)
        plot(1:length(res_norm_rcc(ee,:)),res_norm_rcc(ee,:),'rx');
        hold on 
        plot(1:length(res_norm_rcc(ee,:)),res_norm_rcc(ee,:),'k:');
        hold off
        title(['Epoch ' num2str(ee)]);
        ylabel('Residual norm [m]');
        xlabel('Iteration no.');
        grid on
end
%% Plotting the behaviour of clock correction per iteration
figure(33)
plot(1:length(del_t),del_t,'b:x');
ylabel('Receiver clock correction [s]');
xlabel('Iteration no.');
grid on
title('Evolution of receiver clock correction');

%% Evaluation the difference between measured and precise radial orbit
s_final_rcc = [];
n = 0;
for e = 1:epoch
    s_final_rcc(e,:) = s0_rcc(n+1:n+3,:); % final value of s = [x;y;z] per epoch
    n = n+3;
    r_measured_rcc(e,:) = (norm(s_final_rcc(e,:)) * 10^(-3)); % Measured distance in km.
    dist_measured_rcc(e,:) = r_measured_rcc(e,:) - earth_radius; % Measured distance in radial direction in km.
end

% Plotting measured and observed radial distances in km
dist_precise_rcc = rp - earth_radius;
figure(34)
plot(1:epoch,dist_precise_rcc,'r-',1:epoch,dist_measured_rcc,'b--');
legend('Precise radial distance','Measured radial distance');
xlabel('Epoch');
ylabel('Altitude [km]');
title('Batch processed difference between measured & precise radial dist.');
grid on;
% Plotting the difference between measured and observed orbit
figure(35)
plot(1:epoch,dist_measured_rcc - dist_precise_rcc');
xlabel('Epoch');
ylabel('Difference between measured and precise radial orbit [km]');
title('Batch processed difference between measured & precise orbit');
grid on;

