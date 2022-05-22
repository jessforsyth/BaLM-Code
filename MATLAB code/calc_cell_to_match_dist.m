function [samples,arrays]=calc_cell_to_match_dist(it,data,samples,saving,general,beta,tempering,samp_def,arrays)
%Calculate the cell to match distances for the current values of transformation and permutation vector. 


%first of all calculate the cell to match distances
P=samples.P(1:data.n1);
samples.cell2matchdist=sqrt(sum((data.Y1-samples.fY2(:,P)).^2)); %this is a vector of cell to match distances for each cell in Y1.

%now store the cell2match distances and write data to .bin file
if mod(it,general.write_f)==0
    write_your_data('d',saving, general, beta, tempering,samp_def,arrays,samples,data);
elseif  mod(it,general.thin)==0 && mod(it,general.write_f)~=0
    itap=(it/general.thin)-(floor(it/general.write_f)*(general.write_f/general.thin));
    arrays.cell2matchdist(itap,:)=single(samples.cell2matchdist); 
end

end