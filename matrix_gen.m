function [A,b]=matrix_gen(w)
%Assembles sampled parameter values into matrix A and vector b for
%the affine transformation of data. 

%scaling and shear matrix
sh1=[[cos(w(1)),sin(w(1)),0];[-sin(w(1)),cos(w(1)),0];[0,0,1]]; %to rotate the frame out of standard xyz, and induce a shear when stretched
sh2=[[cos(w(2)),0,-sin(w(2))];[0,1,0];[sin(w(2)),0,cos(w(2))]];
sh3=[[1,0,0];[0,cos(w(3)),sin(w(3))];[0,-sin(w(3)),cos(w(3))]];
sh=sh1*sh2*sh3;


D=[[w(4)+1,0,0];[0,w(5)+1,0];[0,0,w(6)+1]]; %remember to add 1 as we have centered our priors about 0.
S=sh*D;%this matrix now will scale and shear coordinates

%rotation matrix
r1=[[cos(w(7)),sin(w(7)),0];[-sin(w(7)),cos(w(7)),0];[0,0,1]]; %to rotate the frame out of standard xy, and induce a shear when stretched
r2=[[cos(w(8)),0,-sin(w(8))];[0,1,0];[sin(w(8)),0,cos(w(8))]];
r3=[[1,0,0];[0,cos(w(9)),sin(w(9))];[0,-sin(w(9)),cos(w(9))]];
R=r1*r2*r3;  

%complete rotation, scale, shear matrix -A
A=S*R;

%translatoin vector-b
b=[w(10);w(11);w(12)];

end
