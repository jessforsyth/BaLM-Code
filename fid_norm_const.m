function [nlog_normconst]=fid_norm_const(general,f)

    nlog_normconst= -general.d * sum(log(f)); %as in section S1.7
    
end