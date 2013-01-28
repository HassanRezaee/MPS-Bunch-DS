%---------------------------------------------------------------------------------------------
% Author: Hassan Rezaee, Email: h.rezaee@ut.ac.ir
% Paper: "Multiple-Point Geostatistical Simulation Using Bunch-Pasting Direct Sampling Method"
%---------------------------------------------------------------------------------------------

% This function carries out the post-processing calculations on the
% realizations produced by Bunch_DS.m. The post-maps includes E-type,
% Interquartile Range and Quantile (probability 0.5) maps.

%---------------------------------------------------------------------------------------------
% Input parameters:
%   reals: the realizations; each column represents one realization
%   [realz-1 realz-2 realz-3,... realz-n]

%---------------------------------------------------------------------------------------------

% Outputs:
%   post: each post-map is saved in a separate column in the output.

function [post] = sim_post_proc(reals)

    post(:,1)=sum(reals,2)/size(reals,2);
    post(:,2)=iqr(reals,2);
    post(:,3)=quantile(reals,0.5,2);
    
end