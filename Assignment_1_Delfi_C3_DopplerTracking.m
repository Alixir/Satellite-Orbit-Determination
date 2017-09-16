%  AE4872: Satellite Orbit Determination 
%  Homework Assignment 1: Parameter Fitting
%  Author :Ali Nawaz           
%  Facult of Aerospace Engineering, Delft University of Technology.

close all, clear all, clc;

file = importdata('Delfi-C3_32789_201602210946.rre'); % Loading data from the file
time = file.data(:,1);                                % time since recording initiated [s]
freq = file.data(:,2);                                % observed frequency [Hz]
range_rate = file.data(:,3);                          % range rate [m/s]


%% Part 1A

A = zeros(length(time),4); %initializing information matrix A
for a = 1:length(time);
    A(a,1) = 1;            % Column 1 elements of information matrix A
    A(a,2) = time(a);      % Column 2 elements of information matrix A
    A(a,3) = time(a)^3;    % Column 3 elements of information matrix A
    A(a,4) = time(a)^5;    % Column 4 elements of information matrix A
end

% Checking the rank of A matrix.
rank_A = rank(A); % Rank value of 3, A matrix loses rank.

% To prevent this rank loss we'll normalise freq and time; after LSQ
% they'll be restored to orignal states.

% Different normalisation technique tried and the best performing one is
% selected

% freq_norm = ( freq - (max(freq)+min(freq))/2 ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
% time_norm = ( time - (max(time) +min(time))/2 )/( max(time) - min(time));   % Normalised time i.e. |time_i|< 1
% 
% freq_norm = ( freq - mean(freq) ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
% time_norm = ( time - mean(time))/ ( max(time) - min(time));   % Normalised time i.e. |time_i|< 1

freq_norm = ( freq - min(freq) ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
time_norm = ( time - min(time))/ ( max(time) - min(time));   % Normalised time i.e. |time_i|< 1

% Comparing different normalizations
freq_norm1 = ( freq - (max(freq)+min(freq))/2 ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
time_norm1 = ( time - (max(time) +min(time))/2 )/( max(time) - min(time));   % Normalised time i.e. |time_i|< 1

freq_norm2 = ( freq - mean(freq) ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
time_norm2 = ( time - mean(time))/ ( max(time) - min(time));   % Normalised time i.e. |time_i|< 1

freq_norm3 = ( freq - min(freq) ) / ( max(freq) - min(freq)); % Normalised freq. i.e. |freq_i|< 1
time_norm3 = ( time - min(time))/ ( max(time) - min(time));   % Normalised time i.e. |time_i|< 1

% Plot for 3 different frequency time normalizations

figure(19)
plot(time_norm1,freq_norm1,'r-',time_norm2,freq_norm2, 'b--', time_norm3, freq_norm3,'k:');
title('Effect of Time-Frequency normalization');
legend({'mean $\mu$','Matlab mean','min()'},'Interpreter', 'latex');
xlabel('Normalized Time [s]');
ylabel('Normalized Frequency [Hz]');
grid on

% Reconstructing A matrix from the normalised dataset
A_norm = zeros(length(time),4); %initializing normalised information matrix A
for a = 1:length(time_norm);
    A_norm(a,1) = 1;            % Column 1 elements of information matrix A
    A_norm(a,2) = time_norm(a);      % Column 2 elements of information matrix A
    A_norm(a,3) = time_norm(a)^3;    % Column 3 elements of information matrix A
    A_norm(a,4) = time_norm(a)^5;    % Column 4 elements of information matrix A
end

% Rank check for A matrix 

rank_A_norm = rank(A_norm); % Rank value of 4, no rank loss

% Since there is no rank loss, LSQ with psuedo-inverse can be applied

x1 = ((transpose(A_norm)*A_norm))\transpose(A_norm)*freq_norm; % Parameter estimation via LSQ with pinv.

a0_1 = x1(1); % a0
a1_1 = x1(2); % a1
a2_1 = x1(3); % a2
a3_1 = x1(4); % a3

freq_fit_norm1 = A_norm*x1; % Normalised plot-fitted freq plot
% freq_fit1 = freq_fit_norm1* ( max(freq) - min(freq) ) + mean(freq); % De-normalising fitted frequency
% freq_fit1 = freq_fit_norm1* ( max(freq) - min(freq) ) + (max(freq) + min(freq))/2; % De-normalising fitted frequency
freq_fit1 = freq_fit_norm1* ( max(freq) - min(freq) ) + min(freq); % De-normalising fitted frequency

%Application of SVD for verification of LSQ with pseudo-inverse

[U,S,V] = svd(A_norm,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
condition_number = S(1,1)/S(4,4); %2.5508e+14 We could lose upto 15 digits in any numerical calculations when A is computed directly.
S_inv_diag = [inv(S(1,1)), inv(S(2,2)), inv(S(3,3)), inv(S(4,4))]; % Diagonal elements of inv S
S_inv = diag(S_inv_diag);
x2 = V*S_inv*U'*freq_norm; % Parameter estimation via SVD aided LSQ
a0_2 = x2(1); % a0
a1_2 = x2(2); % a1
a2_2 = x2(3); % a2
a3_2 = x2(4); % a3
freq_fit_norm2 = A_norm*x2; % Normalised plot-fitted freq plot
% freq_fit2 = freq_fit_norm2* ( max(freq) - min(freq) ) + mean(freq);
% freq_fit2 = freq_fit_norm2* ( max(freq) - min(freq) ) + (max(freq) + min(freq))/2;
freq_fit2 = freq_fit_norm2* ( max(freq) - min(freq) ) + min(freq); % De-normalised esitmated frequency


figure(11)
plot(time, freq,'r-',time,freq_fit1,'b--', time, freq_fit2, 'k-.');
legend('Observed Frequency [Hz]','Estimated Frequency Pseudoinverse[Hz]','Estimated Frequency SVD [Hz]');
title('Observed freq. vs estimated freq.');
xlabel('Time [s]');
ylabel('Frequency [Hz]');
grid on
%% Part 1B
% Residual plot between observed data and estimated function
res_pinv = freq - freq_fit1; %Residual for LSQ with pinv
res_svd = freq - freq_fit2;  %Residual for LSQ with SVD

figure(12)
plot(time, res_pinv,'r--', time, res_svd, 'b-.');
title('Residual between obs. and est. freq.');
legend( 'Residual frequency via pinv.[Hz]','Residual frequency via SVD.[Hz]');
xlabel( 'Time [s]');
ylabel( 'Frequency [Hz]');
grid on

% Histogram plot of residuals
figure(13)
subplot(2,1,1);
hist(res_pinv,100);
title('Histogram of residuals via pinv');
xlabel('Frequency [Hz]');
ylabel('Occurences');
grid on
subplot(2,1,2);
hist(res_svd,100);
title('Histogram of residuals via SVD');
xlabel('Frequency [Hz]');
ylabel('Occurences');
grid on

% Residual information for LSQ with pinv. case
res_mean_p = mean(res_pinv); % mean of residuals
res_med_p = median(res_pinv); % median of residuals
res_std_p = std(res_pinv); % standard deviation of residuals
res_rms_p = rms(res_pinv); % RMS of residuals
SSE_p = sum( res_pinv.^2);    % Sum squared error performance function
SST = sum( ( mean(freq)-freq).^2); % total sum of squares
res_r2_p = 1-SSE_p/SST;      % Co-efficient of determination R^2 of residuals
% Residual information for LSQ with SVD case 
res_mean_s = mean(res_svd); % mean of residuals
res_med_s = median(res_svd); % median of residuals
res_std_s = std(res_svd); % standard deviation of residuals
res_rms_s = rms(res_svd); % RMS of residuals
SSE_s = sum( res_svd.^2);    % Sum squared error performance function
SST = sum( ( mean(freq)-freq).^2); % total sum of squares
res_r2_s = 1-SSE_s/SST;      % Co-efficient of determination R^2 of residuals
%% Part 1C
prompt = ('Please enter the an integar value of n: ');
n = input(prompt);

An_norm = zeros(length(time_norm), n+1); % initialising normalised A matrix for nth order polynomial

An_norm(:,1) = 1;

for a = 1:length(time_norm)             % Defining normalised A matrix for nth order polynomial
    for b = 1:n
        An_norm(a,b+1) = time_norm(a)^(2*b-1);
    end
end
[Un,Sn,Vn] = svd(An_norm,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
condition_numbern = Sn(1,1)/Sn(end,end); % Indicates the number of potential digit that can be lost in any numerical calculations when A is computed directly.
Sn_inv_diag = [];
for nn=1:length(Sn)
    Sn_inv_diag = [Sn_inv_diag,inv(Sn(nn,nn))];
end
Sn_inv = diag(Sn_inv_diag);
xn = Vn*Sn_inv*Un'*freq_norm; % Parameter estimation via SVD aided LSQ

% xn = ((transpose(An_norm)*An_norm))\transpose(An_norm)*freq_norm; % Parameter estimation via LSQ with pinv.

freq_fit_norm_n = An_norm*xn; % Normalised plot-fitted freq plot
% freq_fitn = freq_fit_norm_n* ( max(freq) - min(freq) ) + mean(freq); % De-normalising fitted frequency
% freq_fitn = freq_fit_norm_n* ( max(freq) - min(freq) ) + ( max(freq) + min(freq))/2; % De-normalising fitted frequency
freq_fitn = freq_fit_norm_n* ( max(freq) - min(freq) ) + min(freq); % De-normalising fitted frequency


figure(14)
plot(time, freq,'r-',time,freq_fitn,'b--');
legend('Observed Frequency [Hz]','Estimated Frequency Pseudoinverse[Hz]');
title(['Observed freq. vs estimated freq. for n = ' num2str(n)]);
xlabel('Time [s]');
ylabel('Frequency [Hz]')
grid on
% Residual information

resn = freq - freq_fitn; % Residual between estimated and observed frequencies.

resn_mean = mean(resn); % Residual mean
resn_med = mean(resn);  % Residual median
resn_std = std(resn);   % Residual standard deviation
resn_rms = rms(resn);   % RMS of residuals
SSEn = sum( resn.^2);    % Sum squared error performance function
SST = sum( ( mean(freq)-freq).^2); % total sum of squares
res_r2n = 1-SSEn/SST;      % Co-efficient of determination R^2 of residuals

%% Part 1D
prompt = ('Please enter the max integar value of n for F-test: ');
n = input(prompt);

rss_n = [];
pol_n = [];
res_fl =[];
for aa = 1:n
    Af = zeros(length(time_norm), aa+1); % initialising normalised A matrix for nth order polynomial
    Af_norm(:,1) = 1;

    for a = 1:length(time_norm)             % Defining normalised A matrix for nth order polynomial
        for b = 1:aa
            Af_norm(a,b+1) = time_norm(a)^(2*b-1);
        end
    end
    [Uf,Sf,Vf] = svd(Af_norm,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
    condition_numbern = Sf(1,1)/Sf(end,end); % Indicates the number of potential digit that can be lost in any numerical calculations when A is computed directly.
    Sf_inv_diag = [];
    for ff=1:length(Sf)
        Sf_inv_diag = [Sf_inv_diag,inv(Sf(ff,ff))];
    end
    Sf_inv = diag(Sf_inv_diag);
    xf = Vf*Sf_inv*Uf'*freq_norm; % Parameter estimation via SVD aided LSQ
%     xf = ((transpose(Af_norm)*Af_norm))\transpose(Af_norm)*freq_norm; % Parameter estimation via LSQ with pinv.

    freq_fit_norm_f = Af_norm*xf; % Normalised plot-fitted freq plot
    % freq_fitf = freq_fit_norm_f* ( max(freq) - min(freq) ) + mean(freq); % De-normalising fitted frequency
    % freq_fitf = freq_fit_norm_f* ( max(freq) - min(freq) ) + ( max(freq) + min(freq))/2; % De-normalising fitted frequency
    freq_fitf = freq_fit_norm_f* ( max(freq) - min(freq) ) + min(freq); % De-normalising fitted frequency
    res_f = freq - freq_fitf; % Residual information
    res_fl = [res_fl,res_f]; % Appending columns of residual per n-th order polyfit
    rss = sum(res_f.^2); % Residual Sum Square
    rss_n = [rss_n;rss]; % Appending residual sum square for nth residual
    pol_n = [pol_n; 2*aa - 1]; % Appending the current order of polynomial
end
F_test = [];
for l = 2:length(rss_n)
    F_test(l,1) = ( ( rss_n(l-1)-rss_n(l) )/(pol_n(l) - pol_n(l-1)) )/ ( rss_n(l)/ (length(time_norm) - pol_n(l))) ;
end

% Plotting a CDF for F-Test
xx = 0:0.00001:5;
z = fcdf(xx,1,500);
for i = 1:length(z)
    Fcrit = xx(i);
    if z(i)>0.95
        break
    end
end

figure(15)
hold on;
plot(xx,z,'b-');
plot(xx, 0.95*ones(size(xx)),'r--');
hold off
title('Cumulative Distribution Function for F-test');
xlabel('F-value');
ylabel('1- p <()');
grid on
% VALIDATION of F_test via vartest2 function

% Vartest2(x,y)  returns a test decision for the null hypothesis that the data in vectors x and y comes from normal distributions with the same variance, using the two-sample F-test. 
%The alternative hypothesis is that they come from normal distributions with different variances. The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.

h_res = []; % Hyptothesis value 
p_res = []; % Probability of Test
ci_res = []; % Confidence interval for the true variance ratio
stats_res = [];% Structure containing info. about test statistic

size_resfl = size(res_fl); % Size of appended residuals
col_resfl = size_resfl(2);  % Total set of residual samples
for hh = 2:col_resfl
    [h,p,ci,stats] = vartest2(res_fl(:,hh),res_fl(:,hh-1),'Alpha',0.05);
    h_res =[h_res;h];
    p_res =[p_res;p];
    ci_res = [ci_res;ci];
    stats_res = [stats_res;stats];
end
% f_val = [stats_res.fstat]';

% plotting probabilities to check which n satisfies p<0.5 
figure(16)
title('F-test results');
subplot(1,2,1)
hold on
% plot(1:length(h_res), f_val,'*');
for zz=1:length(p_res)
    if p_res(zz)>0.05
        plot(zz, p_res(zz),'ro');
    else
        plot(zz,p_res(zz),'b*');
    end
end
% legend('p>0.05','p=<0.05');
title('vartest2');
ylabel('p <');
xlabel('n value');
hold off
grid on
subplot(1,2,2)
hold on 
plot(1:n,F_test,'bo');
index1 = find(z>0.95);
plot(1:0.1:n,xx(index1(1))*ones(size(1:0.1:n)),'r--');
title('CDF');
xlabel('n-value');
ylabel('F-value')
hold off
grid on
% Zoomed in view of CDF based F-Test
figure(17)
hold on 
plot(3:n,F_test(3:end),'bo');
index1 = find(z>0.95);
plot(3:0.1:n,xx(index1(1))*ones(size(3:0.1:n)),'r--');
title('CDF zoomed in range n=[3-15]');
xlabel('n-value');
ylabel('F-value')
hold off
grid on

%% Part 1E
% An optimal of n=5 is obtained. 

A5_norm = zeros(length(time_norm), 5+1); % initialising normalised A matrix for 5th order polynomial
A5_norm(:,1) = 1;

for a = 1:length(time_norm)             % Defining normalised A matrix for nth order polynomial
    for b = 1:5
        A5_norm(a,b+1) = time_norm(a)^(2*b-1);
    end
end
[U5,S5,V5] = svd(A5_norm,'econ'); % Singular value decomposition of A, to avoid singularity errors or solution manifold due to rank deficit
condition_number5 = S5(1,1)/S5(end,end); % Indicates the number of potential digit that can be lost in any numerical calculations when A is computed directly.
S5_inv_diag = [];
for nn=1:length(S5)
    S5_inv_diag = [S5_inv_diag,inv(S5(nn,nn))];
end
S5_inv = diag(S5_inv_diag);
x5 = V5*S5_inv*U5'*freq_norm; % Parameter estimation via SVD aided LSQ

% x5 = ((transpose(A5_norm)*A5_norm))\transpose(A5_norm)*freq_norm; % Parameter estimation via LSQ with pinv.

freq_fit_norm_5 = A5_norm*x5; % Normalised plot-fitted freq plot
% freq_fitn = freq_fit_norm_n* ( max(freq) - min(freq) ) + mean(freq); % De-normalising fitted frequency
% freq_fitn = freq_fit_norm_n* ( max(freq) - min(freq) ) + ( max(freq) + min(freq))/2; % De-normalising fitted frequency
freq_fit5 = freq_fit_norm_5* ( max(freq) - min(freq) ) + min(freq); % De-normalising fitted frequency

dfreq =[]; % Frequency difference between consecutive freq. elements
dtime =[]; % Time difference between consecutive time elements
for c = 2: length(time)
    dfreq = [dfreq;freq_fit5(c) - freq_fit5(c-1)];
    dtime = [dtime; time(c) - time(c-1)];
end
dfreq_dtime = dfreq./dtime;
index_FCA = find( abs(dfreq_dtime)==max(abs(dfreq_dtime))); % Finding the inflection point for closest approach

FCA = freq_fit5(index_FCA); % Frequency of closest approach
TCA = time(index_FCA);      % Time of closest approach

%Verfication of numerical results with visual results

figure(18)
plot(time(2:end),dfreq_dtime);
title('$\frac{dF}{dt}$ vs time','Interpreter','latex');
xlabel('Time [s]');
ylabel('$\frac{dF}{dt}$ [$\frac{1}{s^{2}}$]','Interpreter','latex');
grid on
