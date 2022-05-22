function [data,general, beta, prior, samples, arrays,saving, tempering]=import_files(chain)
%Define what tests you want to run here for each of the chains. We have
%included all of the examples for our runs for the tests included in the
%paper. Just comment out the tests you do/don't want, or include your own. 

%Section 6.1 tests, direct matches , known match , all cells have corresponding match
if chain>=1 && chain<=8
    file1='Y1_8c.txt';
    file2='Y2_8c.txt';
elseif chain>=9 && chain<=16
    file1='Y1_15c.txt';
    file2='Y2_15c.txt';
elseif chain>=17 && chain<=24
    file1='Y1_33c.txt';
    file2='Y2_33c.txt';
elseif chain>=25 && chain<=32
    file1='Y1_62c.txt';
    file2='Y2_62c.txt';
end

%Section 6.2 Data selection tests
%{
if chain>=1 && chain<=8
    file1='Y1_33c_f.txt';
    file2='Y2_33c_f.txt';
elseif chain>=9 && chain<=16
    file1='Y1_62c_f.txt';
    file2='Y2_62c_f.txt';
elseif chain>=17 && chain<=24
    file1='Y1_33c_f_20perc.txt';
    file2='Y2_33c_f_20perc.txt';
elseif chain>=25 && chain<=32
    file1='Y1_62c_f_20perc.txt';
    file2='Y2_62c_f_20perc.txt';
end
%}

%Section 6.3: Non-linear deformation tests
%{
if chain>=1 && chain<=8
    file1='Y1_deftest_pscale=1.txt';
    file2='Y2_33c.txt';
elseif chain>=9 && chain<=16
    file1='Y1_deftest_pscale=1_1.txt';
    file2='Y2_33c.txt';
elseif chain>=17 && chain<=24
    file1='Y1_deftest_pscale=1.txt';
    file2='Y2_33c.txt';
elseif chain>=25 && chain<=32
    file1='Y1_deftest_pscale=1_1.txt';
    file2='Y2_33c.txt';
end
%}


%Section 6.4 : Biological tests with reference markers
%{
if chain>=1 && chain<=8
    file1='Y1_ali_injectiondata.txt';
    file2='Y2_ali_injectiondata.txt';
end
%}

%Section 6.5 : Real movie to stained image tests, example for embryo 1 to
%A,B,C,D. Change to embryo 2,3,4 as well! 
%{
if chain>=1 && chain<=8
    file1='em1.txt';
    file2='emA.txt';
elseif chain>=9 && chain<=16
    file1='em1.txt';
    file2='emB.txt';
elseif chain>=17 && chain<=24
    file1='em1.txt';
    file2='emC.txt';
elseif chain>=25 && chain<=32
    file1='em1.txt';
    file2='emD.txt';
end
%}

%import the data files
[data]=data_import(chain,file1,file2); 

%import the parameters to initiate runs
[general, beta, prior, samples, arrays,saving, tempering,data]=parameters(chain,data,file1,file2);

end
