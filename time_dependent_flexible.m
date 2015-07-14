%--------------------------------------------------------------------- 
%  Computes the frequency dependent dielectric function 
%  from a file with the total dipole at each timestep (assumed to be in Debye) 
%  This is the newer MATLAB version, faster & better than the OCTAVE
%  version!
%  2012 - 2013 Daniel C. Elton 
%  tested as working in 2015
%
%----- Load file ----------------------------------------------------
% input = load('out_128_TTM3F_1_300_F_tot_dip.dat', 'r');
% input = load('TTM3Fdata/out_128_TTM3F_1_300_tot_dip.dat', 'r');
% input = load('TTM3Fdata/seperate_files_which_have_been_added\out_128_TTM3F_1_300_tot_dip.dat', 'r');
% input = load('/home/dan/Dropbox/DIELECTRIC/Dielectric_constant/time_dependent_dielectric/TIP4P2005flex_data/dipcorr_1_340f.xvg', 'r');
% input = load('dipcorr_TIP4P2005f_1_300_.1fs.xvg','r');

%inputsingle = load('/home/delton/out_TTM3F_1_350_256_dip.dat', 'r');



% optional code to convert from single molecule input to total dipole input
ntimesteps = floor(size(inputsingle,1)/256); 
input = zeros(ntimesteps,3);
j = 1;
for i = 1,ntimesteps
	input(i,1:3) = sum(inputsingle(j:j+256,1:3));
	j = j + 256 + 1;
end


%-------------------------- User inputs -------------------------------------------
GROMACS = 0;          %Are we using a GROMACS correlation function input? (not yet tested)
calc_running_eps = 0; %Calculate eps(0) or use inputted value 
epsilon_0 = 73 ;      % 93.2 (300K)   73 (350K) 58.8;    %eps(0) to input
ncells =  1;          %Option for dipolegrid analysis 
latparam =15.674789;  %24.8 ; % 15.674789; % 31.040391; 14.72349;
temp =    350      ;  
timestep = 4*10^(-15);%Timestep in fs converted to seconds!

%supported types of fit:
%1 = one exponential
%7 = stretched exponential
type_fit = 1; 
%The fit parameters are stored in the vector params = [A, tau]
params0 = [1, 20]; % initial trial parameters
params_str0 = [.97,.99,74];

%parameters to be put in by hand !
start =  .6/.004;      % timestep to start fitting
endt  =  4/.004;       %1.5/.0005;   % timestep to stop fitting
endtfit = ntimesteps;  % timestep to stop integrating fit (option for speed)
sizespline = 100;      % width of the spline in timesteps (multiple of 2)

numfreqs = 200;      % number of points to calculate eps_omega at (spaced logrithmically)
minfreq = 1/(ntimesteps*timestep);
maxfreq = 1/(2*timestep);

loaddata    = 1;     %can be a time saver when debugging
calceps0    = 0;     %calculate dielectric constant
eps0_points = 1000;  %number of points for the running eps0 plot
plotting    = 1  ;  %make plots?
printgraphs = 1;     %optional PNG file output of most of the graphs 
printdata   = 1;     %optional printing of the dielectric function data to ascii file
printlibrational = 1;  
printrunningeps0 = 0; 




% %%% other input method: 
% correlation = load('dipcorr_TIP4P2005f_1_300_.1fs.xvg', 'r');
% correlation = correlation(:,2);
% ntimesteps = size(correlation,1); 
 
%%%-----------------------------------------------------------------
%%%----- Caculate the static dielectric constant  -----------------
%%%-----------------------------------------------------------------
if (calceps0 == 1)
  AvgDipSq = sum(sum(input(1:3,:).^2));
  AvgDip   = sum(input(1:3,:),1);
  AvgDipSq = AvgDipSq/ntimesteps;
  AvgDip   = AvgDip/ntimesteps;
  
  kb = 1.3806488e-23; 
  vac_perm = 8.854187817620e-12; 
  Debye2SI = (3.33564e-30)^2;
  vol = (latparam*10^(-10))^3;
  C = Debye2SI/(3*temp*vol*kb*vac_perm);  
  if (GROMACS ~= 1) 
      epsilon_0 = C*(AvgDipSq - sum(AvgDip.^2)) + 1;
  end
  fprintf('epsilon(0) = %f\n',epsilon_0);
end 

%%%-----------------------------------------------------------------
%%%----- Calculate running eps(0) --------------------------------- 
%%%-----------------------------------------------------------------
if (calc_running_eps == 1) 
 eps0_skip = floor(ntimesteps/eps0_points);  
 running_eps_times = eps0_skip*(1:eps0_points)*timestep/( 10^(-15) );
 running_eps_times = running_eps_times./(1000*1000*1000000); %fs -> ns
 running_eps = zeros(eps0_points,1);

 for i = 1:eps0_points 
      nps = eps0_skip*i;      
    
      AvgDipSqT = sum(sum(input(1:nps,1:2).^2));
      AvgDipT = sum(input(1:nps,1:3),1);
      
      AvgDipSqT = AvgDipSqT/nps;
      AvgDipT = AvgDipT/nps;

      running_eps(i) = C*(AvgDipSqT - sum(AvgDipT.^2)) + 1;
 end
end %if (calc_running_eps == 1)

%%%-----------------------------------------------------------------
%%%---------- get correlation function  ---------------------------
%%%-----------------------------------------------------------------
if (GROMACS == 0) 
 corr_components = zeros(ntimesteps,3);
 
 for j = 1:3
	transformed = fft(input(:,j),ntimesteps);
	corr_components(:,j) = ifft(abs(transformed).^2,ntimesteps);
 end
 correlation = sum(corr_components,2);

 norm_vector = (ntimesteps - [1:(ntimesteps-1)]);
 norm_vector(ntimesteps) = 1; %avoid 1/0 and/or abnormally large numbers
 correlation = correlation./norm_vector';
 %second normalization, to normalize to one
 correlation = correlation/correlation(1);
else
 %input (GROMACS) correlation function 
 correlation  = input(:,2);    
end;

%%%-----------------------------------------------------------------
%%%----------- Find fit function ----------------------------------
%%%-----------------------------------------------------------------
%The fit function is of the form corr = A*exp(t/tau)
%Fit is done in picoseconds for accuracy (to avoid small numbers)
%Fitting & spline variables
times_sec = (1:ntimesteps).*timestep;
times_ps = times_sec/(10^(-12));

times_fit_ps = times_ps(start:endt);
corr_fit = correlation(start:endt);

mid = ceil((start+endt)/2);
times_int_sec = times_sec(1:mid);

fo1= fitoptions('Method','NonlinearLeastSquares','tolX',1e-12,'StartPoint',params0);
fo7= fitoptions('Method','NonlinearLeastSquares','tolX',1e-13,'StartPoint',params_str0);

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

mid = ceil((endt+start)/2); %place where spline is joined to orginal function
correlation3 = [correlation(1:mid) ;fitcurve(mid+1:ntimesteps)'];

%% ------ Stich together fit and data with cubic spline 
range =  (  mid - floor(sizespline/2) ):(  mid + floor(sizespline/2)  );
interpx = [(mid - floor(sizespline/2))  ,  (mid - floor(sizespline/4)) , (mid - floor(sizespline/8))  ,(mid + floor(sizespline/8)),(mid + floor(sizespline/4)), (mid + floor(sizespline/2))];

interpy = correlation3(interpx');
spline = interp1(interpx,interpy,range);
correlation3(range) = spline; %add the spline

%%%-----------------------------------------------------------------
%%%------------calculate time-dependent dielectric function -------
%%%-----------------------------------------------------------------
freqs = logspace(log10(minfreq),log10(maxfreq),numfreqs);
%freqs = linspace(minfreq,maxfreq,numfreqs); %linear spacing option
freqs = freqs/(1e12);  %prevents underflow and makes it faster
eps_omega = zeros(1,length(freqs));
timestep_ps = timestep/(1e-12);

for w = 1:length(freqs)
	omega = 2*pi*freqs(w);
	%Do Fourier integration directly
        argument = -([diff(correlation3);0]/timestep_ps)'.*exp(-1i*omega*times_ps);
	eps_omega(w) = trapz(times_ps,argument);
	eps_omega(w) = (epsilon_0 - 1)*eps_omega(w) + 1; 
 
	%%Altnerate no derivative formula (fails at high freqs)
	%argument = correlation3'.*exp(-1i*omega*times_ps);
	%eps_omega(w) = trapz(times_ps,argument);
	%eps_omega(w) = (epsilon_0 - 1)*(1 - 1i*omega*eps_omega(w)) + 1; 

	fprintf('%i %i\n',w,length(freqs));
end
freqs = freqs*1e12/(2.99*10^10);%convert to cm^-1

%-------------- infrared spectrum ------------------------
% freqs = logspace(log10(minfreq),log10(maxfreq),numfreqs);
% freqs = freqs/(1e12);  %prevents underflow and makes it fasters
% infrared = zeros(1,length(freqs));
% 
% for w = 1:length(freqs)
% 	omega = 2*pi*freqs(w);
% 	
%    %intfun = interp1(times_sec,argument,'cubic');
%    %eps_omega(w) = integral(intfun,
%     timestep_ps = timestep/(1e-12);
%     argument = -correlation3'.*exp(-1i*omega*times_ps);
% 	infrared(w) = trapz(times_ps,argument);
% 
% 	fprintf('%i %i\n',w,length(freqs));
% end
% freqs = freqs*1e12;

%%%-----------------------------------------------------------------
%%% -------------------------- Plotting --------------------------- 
%%%-----------------------------------------------------------------
if plotting == 1 
 %%  Plotting the running average of eps(0) 
 if (calc_running_eps == 1) 
	figure(1)
	clf;
	plot(running_eps_times,running_eps);
	title('Running average of eps(0)')
	xlabel('Time (ns)');
	ylabel('Eps(0)');
	if printgraphs == 1 
		saveas(1,'running_average.png');
	end
 end

 %% Plotting exponential decay and fit of the corr. function  
 figure (2); clf;
 semilogy(times_ps([1:2*endt]),correlation([1:2*endt]),times_ps([1:2*endt]),correlation3([1:2*endt]),times_ps([1:2*endt]),fitcurve(1:2*endt));
 legend('Data', 'Data + Fit + Spline', 'Fit')
 title('Dipole time autocorrelation function')
 xlabel('Time (ps)')
 str3 = sprintf('Fitting function: A*exp(-t/tau)\n tau = %12.6e, A = %7.4e \n R^2 = %12.6e',tau,A,R);
 if type_fit == 7 
    str3 = sprintf('Fitting function: A*exp(-t/tau)\n tau = %12.6e, A = %6.4e, Beta = %6.4e \n R^2 = %12.6e',tau,A,beta,R);
 end

 % locate text near middle of figure
 text(.6*times_ps(endt),.75*max(correlation(1:endt)),str3);

 if printgraphs == 1 
    saveas(2,'correlation_shoulder_flex.png');
 end 

 %% - Plotting short time librational part of corr. function  
 figure (3); clf;
 plot(times_ps(1:1000),correlation(1:1000),times_ps(1:1000),fitcurve(1:1000));
 legend('Data', 'Fit')
 title('Dipole time autocorrelation function')
 xlabel('Time (ps)')
 if printgraphs == 1 
    saveas(3,'librational_flex.png');
 end

 %%   Plotting the entire time correlation function 
 figure(4); clf;
 title('Dipole time autocorrelation function')
 xlabel('Time (ps)')
 plot(times_ps,correlation,times_ps,fitcurve);
 legend('Data', 'Fit')

 %% --- Plotting the freq-dependent dielectric function  
 figure(5); clf;
 loglog(freqs,real(eps_omega),freqs,-imag(eps_omega));
 legend('Real part','Imaginary part');
 title('Frequency-dependent dielectric constant');
 xlabel('Frequency (Hz)');
 ylabel('Epsilon');
 if printgraphs == 1 
    saveas(5,'eps_omega_flex.png'); 
 end

 %% Plotting a Cole-Cole plot 
 figure(6); clf;
 plot(real(eps_omega),-imag(eps_omega));
 title('Cole-Cole plot');
 xlabel('real(epsilon)');
 ylabel('imag(epsilon)');
 if printgraphs == 1 
     saveas(6,'cole_cole_flex.png'); 
 end

end%if plotting == 1 

%%%-----------------------------------------------------------------
%%%------------------ write data to files -------------------------
%%%-----------------------------------------------------------------
if printdata == 1 
	matrix2save = [freqs', real(eps_omega)', -imag(eps_omega)'];
	%header = sfprintf('#Time dependent dielectric funtion - frequency - real part - imaginary part for %s',inputfilename);
	%save_header_format_string (header);
	save -ascii eps_omega_TTM3F_350.dat matrix2save 	
end
if printlibrational == 1 
	matrix2save = [times_ps(1:2000)', correlation(1:2000)];
	save -ascii librational_flex_350.dat matrix2save 	
end
if printrunningeps0 == 1 
	matrix2save = [running_eps_times', running_eps];
	save -ascii running_eps0.dat matrix2save 	
end


