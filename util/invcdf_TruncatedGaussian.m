function [val] = invcdf_TruncatedGaussian(cdf,x_bar,sigma_bar,lowerBound,upperBound)
%stats_TruncatedGaussian - stats for a truncated gaussian distribution
% INPUT
%   - cdf: evaluated at the values at cdf
%   - x_bar,sigma_bar,lowerBound,upperBound: suppose X~N(mu,sigma^2) has a normal distribution and lies within
%   the interval lowerBound<X<upperBound
%   *The size of cdfTG and pdfTG is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
% OUTPUT
%   - val: the x corresponding to the cdf TG value
% Written by: X. Du 05/20/2015

term2 = (normcdf_my(upperBound,x_bar,sigma_bar) - normcdf_my(lowerBound,x_bar,sigma_bar));
const2 =  cdf*term2 + normcdf_my(lowerBound,x_bar,sigma_bar);

const3 = const2*2-1;
inner_temp = erfinv(const3);
val = inner_temp*sqrt(2)*sigma_bar + x_bar;


end


function [p] = normcdf_my(x,mu,sigma)
% INOPUT
%   - x: the x to compute the cdf value
%   - mu: mean of the Gaussian distribution
%   - sigma: sigma of the Gaussian distribution
% OUTPUT
%   - val: the x corresponding to the cdf value
% based on MATLAB normcdf() function.

z = (x-mu) ./ sigma;
p = 0.5 * erfc(-z ./ sqrt(2));


end