%Bayesian Landmark Matching (BALM) pipeline
%J.E.Forsyth May 2022
%Joint work with Prof. S.Cotter 

%this is a pipeline to register two sets of landmarks via affine, (and
%non-linear) transformation, infer on the permutation vector. Inclusion of
%Bayesian data selection to account for points with non-corresponding
%matches. 

%This pipeline is set up to run on x=32 cores, but can easily be changed to
%suit the system being used. 

function run_multiple_chains()
close all 
clear
clc

nc=32; %number of chains - choose to fit your system e.g. number of cores
parpool(nc) %intiialise matlab workers for the number of chains 
fprintf('\n')

parfor chain=1:nc
    run_matching(chain) %now run this individual chain (each chain is allocated to individual worker)
end
fprintf('\n')

delete(gcp('nocreate')) %close all the workers down after you've finished

end