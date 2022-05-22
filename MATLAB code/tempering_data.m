function [tempering]=tempering_data(Nt,tempering, data, prior, general,samp_trans,samp_def,samp_perm,samp_fid)
%This function caluclates an appropriate start temperature for the
%tempering, as well as an appropriate cooling rate based on the initial
%number of iterations described.

nlog_likelihood_array=zeros(Nt,1); %empty array for storing the nlogs of the posterior predictive

for i=1:Nt   %number of random samples to perform
    %sample random values for Theta from priors
    if samp_trans==1
        phis1=(rand(1,3)-0.5)*(2*pi());  %sample from priors on affine params
        phis2=(rand(1,3)-0.5)*(2*pi());
        s=normrnd(0,prior.sigs,1,3);
        b=normrnd(0,prior.sigb,1,3);
        theta=[phis1,s,phis2,b];
        [theta_A,theta_b]=matrix_gen(theta); %generate affine A matrix and b vector        
        if samp_def==1            
            p=mvnrnd(zeros(3,data.n2),prior.mom_C); %sample from prior on momentum 
            
            y0=[single(reshape(data.Y2',3*data.n2,1));reshape(p',3*data.n2,1)]; %assemble Y2 and momentum into column vector for input into ode solver
            [~,yout]=ode45(@(t,y) p_deform_mex(t,y,single(prior.kernel),single(data.n)),single([0,1]),y0);  %perform deformation using compiled p_deform function, tspan=[0,1], we use singles to further improve performance here;
            yfinal=yout(end,:)'; %deformed Y2
            Y2d=double(reshape(yfinal(1:(3*data.n2)),data.n2,3)'); 
            fY2=(theta_A*Y2d)+theta_b; %apply affine transformation to the deformed Y2
        else
            fY2=(theta_A*data.Y2)+theta_b; %if not applying non-linear deformation, then just apply the affine transformation
        end     
    else
        fY2=data.Y2; %if doing neither affine or non-linear, just set equal to Y2
    end
    
    if samp_perm==1
        P=[randperm(data.n2-data.n_refpoints),linspace(data.n2-data.n_refpoints,data.n2,data.n_refpoints)]; %generate random permutation P
    else
        P=linspace(1,data.n2,data.n2); %or just use linear permutation vector
    end
    
    
    if samp_fid==1        
        f=betarnd(prior.fid_alpha,prior.fid_beta,1,data.n1); %sample random fidelity terms from prior     
        f(f==1)=1-1e-10;                                     %impose these assignments so that we don't have exact 0s and 1s
        f(f==0)=1e-10;
        D=diag(f).^2;                                        %construct fidelity matrix 
        
        [nlog_normconst]=fid_norm_const(general,f);          %construct fidelity normalisation constant dependent on D
    else
        f=ones(1,data.n1);
        D=diag(f);
        [nlog_normconst]=fid_norm_const(general,f); %
    end
    
    
   
    nlog_likelihood_array(i,1)=neg_log_PP(data.Y1,data.n1,fY2,D,general.Psi,general.nu,P) + nlog_normconst; %now calculate the negative log of the posterior predictive and append
    %store in array so we can calc start temp and cooling rate. 
    
end

P95=prctile(nlog_likelihood_array,95); %calculate the 95th and 5th percentiles of the nlogPP array
P5=prctile(nlog_likelihood_array,5);

tempering.T0=(P95-P5)/log(1+tempering.tol);                %calculate start temperature
tempering.tc=nthroot(1/tempering.T0,general.N/general.fc); %calculate cooling rate

tempering.betac= (tempering.tc)^0.5;                       %beta adustment due to cooling
tempering.T=tempering.T0;                                  %initialise start temperature

end
    
