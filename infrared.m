%infrared.m 
%simple infrared program 

%----- Load file ----------------------------------------------------
input = load('../../infrared/out_TTM3F128_1_300_10ps_.5fsNVE_tot_dip.dat', 'r');
%input = load('../../infrared/out_128_TTM3F_1_300_F_tot_dip.dat', 'r');
% input = load('TTM3Fdata/128_1_300td.dat', 'r');
ntimesteps = size(input,1); 


%Fitting parameters to be put in by hand !
start = 200;      % timestep to start fitting
endt   = 400;   % timestep to stop fitting
endtfit = ntimesteps; % timestep to stop integrating fit (for speed)
sizespline = 20;  % width of the spline in timesteps (multiple of 2)

ncells =  1;          %Option for analyzing multiple dipoles
latparam = 15.674789; % 31.040391 ;  %24.8 
temp = 300        ;   %Temp (K) 
timestep = .5*10^(-15);%Timestep in fs => converted to seconds!

type_fit = 1; 
%supported types of fit:
%1 = one exponential
%7 = stretched exponential

        
%----- Constants & things -------------------------------------------
kb = 1.38e-23; 
hbar = 6.626e-34; 
c = 3*10^8;
vol = latparam^3;
eps0 = 8.85418782e-12;
Debye2SI = 3.33564e-30; %Debye -> SI units (Coloumb meters)

prefac = Debye2SI^2/(6*c*vol*eps0*kb*temp);

AvgDipSq = sum(sum(input.^2));
AvgDipSq = AvgDipSq/ntimesteps;

%----- Calculate correlation function  ------------------------------
NFFT =  2^nextpow2(2*ntimesteps);
input = [input; zeros(ntimesteps,3)];%be extra safe and make sure its padded with zeros

corr_components = zeros(NFFT,3);


for j = 1:3
	transformed = fft(input(:,j),NFFT);
	corr_components(:,j) = ifft(abs(transformed).^2,NFFT);
end
correlation = sum(corr_components,2);

correlation2 = correlation(1:ntimesteps);

norm_vector = (ntimesteps - (1:(ntimesteps-1)));
norm_vector(ntimesteps) = 1; %avoid 1/0 and/or abnormally large numbers
correlation2 = correlation2./norm_vector';
% second normalization - divide by <M^2>
correlation2 = correlation2/AvgDipSq;

%------------- Fitting & Spline variables ----------------------------
times_sec = (1:ntimesteps).*timestep;
times_ps = times_sec/(10^(-12));

times_fit_ps = times_ps(start:endt);
corr_fit = correlation2(start:endt);

mid = ceil((start+endt)/2);
times_int_sec = times_sec(1:mid);

% %----------- Find fit function --------------------------------------
% %The fit function is of the form corr = A*exp(t/tau)
%The fit parameters are stored in the vector params = [A, tau]
params0 = [1, 2]; % initial trial parameters
params_str0 = [.98,10,.98];
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
    [thefit,gof2] = fit(times_fit_ps',corr_fit,ft2);
    A = thefit.a
	tau = thefit.t
    beta = thefit.b
    R = gof2.rsquare
	fitcurve = A*exp(-(times_ps/tau)^beta);
end

mid = ceil((endt+start)/2); %place where spline is joined to orginal function
correlation3 = [correlation2(1:mid) ;fitcurve(mid+1:ntimesteps)'];

%-------------- Stich together fit and data with cubic spline -------------
range = (  (mid - sizespline/2 - 100):(mid + sizespline/2) + 200  );
interpx = [(mid - sizespline/2) - 100,(mid - sizespline/2)-50, (mid - sizespline/2),(mid + sizespline/2),(mid + sizespline/2) + 100, (mid + sizespline/2) + 200];

interpy = correlation3(interpx');
spline = interp1(interpx,interpy,range);
correlation3(range) = spline; %add the spline
correlation3= correlation3*AvgDipSq; %undo normalization

%-------------- Do the FFT! -----------------------------------------------
%The scheme for finding the Fourier transform is taken directly from the 
%MATLAB website page on fft(). The fft() function automatically adds zeros. 
L = length(correlation3);
NFFT =  2^nextpow2(L);

Y = fft(correlation3,NFFT)/L;

f = (1./(timestep*2) )*linspace(0,1,NFFT/2+1); %frequencies in Hz
cms = f/(3*10^10); %inverse centimeters

f = 2*pi*f; %convert f to angular frequency omega
Y = prefac*Y(1:NFFT/2+1).*(f.^2)'; 

%plot only wavenumbers less than 8000
maxf = floor((NFFT/2+1)*8000*3*10^10*timestep*2); 

cms  = cms(1:maxf);
Y = Y(1:maxf);

%--------- Plot the IR spectra --------------------------------------------
figure(1); clf;
plot(cms,abs(Y))
xlabel('Wavenumber (cm^-1))')

%--------- Plot the correlation function + fit ----------------------------
figure (2); clf;
semilogy(times_ps([1:2*endt]),correlation2([1:2*endt]),times_ps([1:2*endt]),correlation3([1:2*endt])/AvgDipSq,times_ps([1:2*endt]),fitcurve(1:2*endt));
legend('Data', 'Data + Fit + Spline', 'Fit')
title('Dipole time autocorrelation function')
xlabel('Time (ps)')
%saveas(2,'autocorr.png');


%--------- Plot the entire correlation function ---------------------------
figure (3); clf;
semilogy(times_ps,correlation2);
title('Dipole time autocorrelation function')
xlabel('Time (ps)')




