function [general]=trans_t1covariance_gen(samples, general, it)
%calculate the rolling average of the empirical sample covariance for the transformation (affine and non-linear deformation) parameters. 
%see section S1.10

    deln=(samples.transformation-samples.avg_transformation); %difference between mean and it-1 vals.
    %this avg transformation array is from it-1 as it has NOT yet been updated
    general.C_transformation_temp = (((it-2)/(it-1)) .* general.C_transformation_temp) + ((1/it).*(deln*deln'));
   
end 
