%%---------------------------------------------------------------------------------
%  Computes the correlation functions for a set of individual dipoles
%  from a file with the all the dipoles at each timestep (assumed to be in Debye) 
%  This is the newer MATLAB version, faster & better than the OCTAVE
%  version!
%  Copyright 2012-2013 Dan Elton 
%----- Load file ----------------------------------------------------
allinput = load('C:\Users\Owner\Dropbox\DIELECTRIC\DipoleGrid\dgrid2_10000_3Ang_10648.dat', 'r');

%-------------------------- User inputs -------------------------------------------
ncells =10648;          %Option for dipolegrid analysis 
ntimesteps = floor(  size(allinput,1)/ncells  );

timestep =  1000*10^(-15);%Timestep in fs => converted to seconds!

type_fit = 1; 
%supported types of fit:
%1 = one exponential
%7 = stretched exponential

%Fitting parameters to be put in by hand !
start = 3;      % timestep to start fitting
endt   = 8;   % timestep to stop fitting

plotting =  1  ;        %make plots?
printgraphs = 1;       %optional PNG file output of most of the graphs 
printdata   = 1  ;     %optional printing of the dielectric function data to ascii file
printlibrational = 1;  
printrunningeps0 = 1;  

params0 = [1, 4]; % initial trial parameters
params_str0 = [1,.9,27];

input = zeros(ntimesteps,3,ncells);

for j = 1:3
    for i = 1:ntimesteps
        starti = (i-1)*ncells + 1 ;
        endi = (i-1)*ncells+ncells;
        input(i,j,:) = allinput(starti:endi,j );
    end
end
% correlation = load('dipcorr_TIP4P2005f_1_300_.1fs.xvg', 'r');
% correlation = correlation(:,2);
% ntimesteps = size(correlation,1); 
% epsilon_0 = 61;



%----- Calculate correlation function  ------------------------------
correlationA = zeros(ntimesteps,ncells);

for i = 1:ncells 
    
    corr_components = zeros(ntimesteps,3);
    for j = 1:3
        transformed = fft(input(:,j,i),ntimesteps);
        corr_components(:,j) = ifft(abs(transformed).^2,ntimesteps);
    end
    correlationA(:,i) = sum(corr_components,2);
    
    AvgDipSq = sum(sum(input(:,:,i).^2))/ntimesteps; %divide by <M^2>
    
    correlationA(:,i) = correlationA(:,i)/AvgDipSq;
end

correlation = sum(correlationA,2)/ncells;

norm_vector = (ntimesteps - (1:(ntimesteps-1)));
norm_vector(ntimesteps) = 1; %avoid 1/0 and/or abnormally large numbers
correlation = correlation./norm_vector';

%Print out average dipole moment--------------------------------------

%AvgDip = sum(allinput(:,4))/(ntimesteps*ncells)
AvgDip = sum(sqrt(sum(allinput(:,1:3).^2,2)))/(ntimesteps*ncells)

%%

%------------- Fitting & spline variables ----------------------------
times_sec = (1:ntimesteps).*timestep -timestep;
times_ps = times_sec/(10^(-12));

times_fit_ps = times_ps(start:endt);
corr_fit = correlation(start:endt);

mid = ceil((start+endt)/2);
times_int_sec = times_sec(1:mid);

%----------- Find fit function --------------------------------------
%The fit function is of the form corr = A*exp(t/tau)
%The fit parameters are stored in the vector params = [A, tau]

%Fit is done in picoseconds for accuracy (to avoid small numbers)

fo1= fitoptions('Method','NonlinearLeastSquares','tolX',1e-12,'StartPoint',params0);
fo7= fitoptions('Method','NonlinearLeastSquares','tolX',1e-12,'StartPoint',params_str0);

ft1 = fittype('a*exp(-x/t)','options',fo1);
ft7 = fittype('a*exp(-(x/t)^b)','options',fo7);

if type_fit == 1    
    [thefit,gof2] = fit(times_fit_ps',corr_fit,ft1);
    A = thefit.a
	tau = thefit.t
    R = gof2.rsquare
	fitcurve = A*exp(-times_ps/tau);
end
if type_fit == 7 
    [thefit,gof2] = fit(times_fit_ps',corr_fit,ft7);
    A = thefit.a
	tau = thefit.t
    beta = thefit.b
    R = gof2.rsquare
	fitcurve = A*exp(-(times_ps/tau).^beta);
end

%-------------------------- Plotting ------------------------------------- 
if plotting == 1 
%--- Plotting exponential decay and fit of the corr. function ------- 
figure (2); clf;
semilogy(times_ps([1:2*endt]),correlation([1:2*endt]),times_ps([1:2*endt]),fitcurve(1:2*endt));
legend('Data', 'Data + Fit + Spline', 'Fit')
title('Dipole time autocorrelation function')
xlabel('Time (ps)')
str3 = sprintf('Fitting function: A*exp(-t/tau)\n tau = %12.6e, A = %7.4e \n R^2 = %12.6e',tau,A,R);
if type_fit == 7 
    str3 = sprintf('Fitting function: A*exp(-t/tau)\n tau = %12.6e, A = %6.4e, Beta = %6.4e \n R^2 = %12.6e',tau,A,beta,R);
end
%locate text near middle of figure
text(.6*times_ps(endt),.75*max(correlation(1:endt)),str3);
if printgraphs == 1 
    saveas(2,'correlation_shoulder_flex.png');
end 

%--- Plotting short time librational part of corr. function ----------
figure (3); clf;
plot(times_ps(1:500),correlation(1:500),times_ps(1:500),fitcurve(1:500));
legend('Data', 'Fit')
title('Dipole time autocorrelation function')
xlabel('Time (ps)')
if printgraphs == 1 
    saveas(3,'librational_flex.png');
end

%--- Plotting the entire time correlation function -------------------
% figure(4); clf;
% title('Dipole time autocorrelation function')
% xlabel('Time (ps)')
% plot(times_ps,correlation,times_ps,fitcurve);
% legend('Data', 'Fit')

end%if plotting == 1 

% %------------------ print data to files --------------------------------- 
% if printlibrational == 1 
% 	matrix2save = [times_ps', correlation(1:2000)];
% 	save -ascii librational_flex.dat matrix2save 	
% end


%----- make histgram -----------------------------------------------
[y,x] = hist(sqrt(sum(allinput.^2,2)),200);
matrix2save = [x' , y'];
save -ascii histogram.dat matrix2save 	

