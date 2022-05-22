function [general, beta, prior, samples, arrays,saving, tempering,data]=parameters(chain,data,file1,file2)
%Parameters file, most parameters are initiated here, change them here for
%different runs. You can make chain specific values of parameters by
%including an if loop as in the import_data.m file. e.g. 
% if chain>=1 && chain<=8
%       parameter_x == ....
%       ...and so on 
% end

general=[];
general.N=7e6;                       %minimum number of tempered iterations
general.N1=1e6;                      %number of iterations at T=1
general.fc=2e3;                      %check frequency (for checking acceptance rate etc.)
general.thin=100;                    %thinning of data at T=1
general.write_f=10e3;                %write data frequency (just to save time on writing out data)
general.Nt=general.fc/10;            %number of samples to estimate start temperature from
general.alpha_bar_t=1;               %average acceptance rate transformation
general.alpha_t=zeros(1,general.fc); %array of alpha vals transformation
general.alpha_bar_f=1;               %average acceptance rate fidelity
general.alpha_f=zeros(1,general.fc); %array of alpha vals fidelity
general.alpha_bar_P=1;               %average acceptance rate permutation
general.alpha_P=zeros(1,general.fc); %array of alpha vals permutation
general.acc_opt=0.234;               %optimum acceptance rate
general.acc_tol=0.1;                 %tolerance of acceptance rate
IW_mean= 0.1^2;                      %mean of hyperprior on Sigma (mean of IW dist)
IW_var=  0.2^2;                      %variance of hyperprior on Sigma (variance of IW dist)
general.nu=((2*IW_mean^2)/(IW_var)) + 6 ;         %hyperparameter 
general.Psi=(general.nu - 4)* (IW_mean).*eye(3);  %hyperparameter 

general.writing_count=0;             %keep a check of how many times we write data during T=1, this is to assist in data unpacking later on
general.adaptive_on=0;               %make sure adaptive MCMC is off for tempered sampling (this is changed later on)
general.d=3;                         %number of dimensions of problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=[];       %step size related parameters involved in RW of parameters
beta.t=1;      %step-size term for transformation sampling %%%% later we scale as 2.38^2/d where d=12 or 12+3n
beta.f=1;      %step-size for fidelity term sampling  %%%% later we scale as 2.38^2/d where d=n1
beta.adj=0.05; %adjustment factor for beta when acceptance rate is not within limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prior=[];                                          %structure ot contain all of the prior data on params
Xtheta=unifrnd(-pi,pi,1,10000);                    %Samples on uniform dist (prior on angles) bounded params
Xtautheta=log( ((2*pi)./(Xtheta+pi)) -1 );         %transform them to transformed params
prior.ang_propstd = std(Xtautheta) ;               %calc. approximate variance for proposals on bounded terms
prior.sigb=0.1;                                    %sigma for b in transformation
prior.sigs=0.1;                                    %sigma for scaling in transformation
prior.sigmom=1                                     %sigma for momentum prior
prior.kernel=1;                                    %kernel in deformation calculation to do with the repulsion of cells (~distance before cels repel each other)
prior.mom_C=diag(ones(1,data.n2)*prior.sigmom^2);  %intial covariance matrix for momentum proposals
prior.fid_alpha= 2;                                %fidelity prior alpha term (for beta dist)
prior.fid_beta=2;                                  %fidelity prior beta term (for beta dist)
Xf=betarnd(prior.fid_alpha,prior.fid_beta,1,10000);%Samples on beta dist (prior on fidelity) bounded params
Xf(Xf==1)=1-1e-10;
Xf(Xf==0)=1e-10;
Xtauf= log((1./Xf)-1);                             %transform them to transformed params
prior.fid_propstd = std(Xtauf);                    %calc. approximate variance for fidelity proposals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samples=[];                         %structure to store current values of parameters
samples.alpha_t=1;                  %current values of acceptance ratio for transformation
samples.alpha_f=1;                  %current values of acceptance ratio for fidelity
samples.alpha_P=1;                  %current values of acceptance ratio for permutation
samples.nlog_normconst=0;           %normalisation constant of nlog posterior (dep. on fidelity)
samples.avg_nlog_posterior=0;       %average nlog of posterior
samples.avg_A=zeros(1,9);           %avg. Affine matrix
samples.avg_b=zeros(1,3);           %avg. b vector
samples.avg_p=zeros(3,data.n2);     %avg. momenta
samples.avg_fid=zeros(data.n1,1);   %avg. fidelity
samples.avg_nlog_normconst=0;       %avg. normalisation constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arrays=[];                                                    %array structure to save out data before writing to files
arrays.temp=zeros(general.write_f/general.thin,1,'single');   %temperatures
arrays.accr_t=zeros(general.write_f/general.thin,1,'single'); %acceptance ratio transformation
arrays.beta_t=zeros(general.write_f/general.thin,1,'single'); %step-size vals transformation
arrays.accr_f=zeros(general.write_f/general.thin,1,'single'); %acceptance ratio fidelity
arrays.beta_f=zeros(general.write_f/general.thin,1,'single'); %step-size vals fidelity
arrays.accr_p=zeros(general.write_f/general.thin,1,'single'); %acceptance ratio permutation
arrays.theta=zeros(general.write_f/general.thin,12,'single'); %store transformation parameters
arrays.A=zeros(general.write_f/general.thin,9,'single');      %store affine matrix vals
arrays.f=zeros(general.write_f/general.thin,data.n1,'single');%store fidelity params
arrays.matches=zeros(data.n1,data.n2);                        %n1 x n2 matrix to store counts of matching during T=1
arrays.px=zeros(general.write_f/general.thin,data.n2,'single');%store momentum x, y, z components
arrays.py=zeros(general.write_f/general.thin,data.n2,'single');
arrays.pz=zeros(general.write_f/general.thin,data.n2,'single');
arrays.cov_update=[];                                         %store iterations where covariance is updated
arrays.cell2matchdist=zeros(general.write_f/general.thin,data.n1,'single');%store the cell to match distances during T=1

saving=[];                                       %structure to save all saving related information
saving.files=[];                                 %empty array to store all ecessary .bin file names for creation
ftime=clock;                                     %time of starting to help name files
saving.directory=strcat('./Results/',sprintf('%d%d%d_%d_%s,%s/',ftime(3),ftime(2),ftime(1),chain,file1(1:end-4),file2(1:end-4)));
saving.matches= repmat('%d ',1,data.n2);         %saving format for matches array
saving.matches=[saving.matches '\n'];
saving.matches=repmat(saving.matches,1,data.n1);

tempering=[];           %structure to contain all tempering related info
tempering.tol=0.01;     %user-set tolerance to define start temperature and cooling rate
end
