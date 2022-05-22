function [data]=data_import(chain,filename1,filename2)
%This function imports the relevant data files for the current chain running. 

data=[];                                        %initialise data structure to store all of your data in
path='Data';                                    %path to the data files 

%get file 1 data
file1=filename1;
fname=fullfile(path,file1);
fileID = fopen(fname,'r');
Y1=fscanf(fileID,'%f');
Y1=transpose(reshape(Y1,[length(Y1)/3,3]));     %coordinates for Y1 imported as 3xn1 array
fclose(fileID);
data.Y1=Y1;                                     %assign data.Y1=Y1

%get file 2 data
file2=filename2;
fname=fullfile(path,file2);
fileID = fopen(fname,'r');
Y2=fscanf(fileID,'%f');
Y2=transpose(reshape(Y2,[length(Y2)/3,3]));     %coordinates for Y1 imported as 3xn1 array
fclose(fileID);
data.Y2=Y2;                                     %assign data.Y1=Y1

data.n1=length(data.Y1);                        %calculate how many cells in Y1
data.n2=length(data.Y2);                        %calculate how many cells in Y2
%check that n1<n2 , if not ....swap the data round
if data.n1>data.n2  
    temp=data.Y2;
    tempn=data.n2;
    data.Y2=data.Y1;
    data.n2=data.n1;
    data.Y1=temp;
    data.n1=tempn;
    fprintf('Data swapped \n');
end
data.n=data.n2;   %make n==n2 (i.e. max cell number)                         


%Re-scale Y1 and Y2 (this tries to account for any gross shifts/transformation of the data.
%centre points about 0
data.Y1=data.Y1-mean(data.Y1,2);
data.Y2=data.Y2-mean(data.Y2,2);

%change data to have a median distance between neighbours =1
nbr_distances_1=zeros(1,data.n1);
for i=1:data.n1
    [~,D] = knnsearch(data.Y1', data.Y1(:,i)','K',2);
    nbr_distances_1(i)=D(2);
end
nbr_distances_2=zeros(1,data.n2);
for i=1:data.n2
    [~,D] = knnsearch(data.Y2', data.Y2(:,i)','K',2);
    nbr_distances_2(i)=D(2);
end

cell_dist_scale1=min(nbr_distances_1); %minimum cell-cell distance
cell_dist_scale2=min(nbr_distances_2); %minimum cell-cell distance

data.Y1=data.Y1./cell_dist_scale1; %re-scale Y1
data.Y2=data.Y2./cell_dist_scale2; %re-scale Y2

end