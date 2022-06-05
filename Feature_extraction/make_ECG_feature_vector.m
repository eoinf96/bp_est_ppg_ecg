function [feature_vector, feature_names] = make_ECG_feature_vector(ts, peaks)

    
scale_list = 2:2:8;
     
     
if ~isempty(peaks)
    [mobility, complexity ] = HjorthParameters(ts);
    fractal = Higuchi_FD(ts', 17);
    
    p = histcounts(ts, 'Normalization', 'Probability');
    shannon_entropy = -sum(p(p>0).*log2(p(p>0)) );

    
    % Complexity features
    m = 2;
    r_factor = 0.2;
    approxEnt = approximateEntropy(ts, 'Dimension', m, 'Radius', r_factor * std(ts));
    sampEnt = sampleEntropy(ts, m, r_factor );
   
    for scale_idx = 1:length(scale_list)
        mseEnt(scale_idx) = multiscaleSampleEntropy( ts, m, r_factor, scale_list(scale_idx));
    end
    
    feature_vector = [mobility; complexity; fractal; shannon_entropy; approxEnt; sampEnt; mseEnt(:)];
else
    feature_vector = nan(6+length(scale_list), 1);
end

% if do_HRV

% RR.ts = t(peaks(2:end)) - t(peaks(1:end-1));
% RR.t = t(peaks(2:end));
% 
% HRV_feats = func.HRV.get_HRV_features(RR.ts, RR.t);
% HRV_feats = rmfield(HRV_feats, intersect(fieldnames(HRV_feats), {'t'; 'VLF_power';'nLF_power';'nHF_power';'T_power';'samp_entropy'; 'LF_div_HF'}));
% feature_vector = [feature_vector; struct2array(HRV_feats)'];

% feature_names = vertcat({'Mobility'; 'Complexity'; 'Fractal'; 'Shannon_Entropy'; 'approxEnt'; 'sampEnt' }, strcat('MSE_scale_', cellstr(string(scale_list)))');
feature_names = vertcat({'Hjorth_mobility'; 'Hjorth_complexity'; 'Fractal_dimension'; 'SE'; 'approxEnt'; 'sampEnt' }, strcat('MSE_scale_', cellstr(string(scale_list)))');

% feature_names = vertcat(feature_names, fieldnames(HRV_feats));
end



function value = sampleEntropy(signal, m, r, dist_type)
    % Error detection and defaults
    if nargin < 3, error('Not enough parameters.'); end
    if nargin < 4
        dist_type = 'chebychev';
%         fprintf('[WARNING] Using default distance method: chebychev.\n');
    end
    if ~isvector(signal)
        error('The signal parameter must be a vector.');
    end
    if ~ischar(dist_type)
        error('Distance must be a string.');
    end
    if m > length(signal)
        error('Embedding dimension must be smaller than the signal length (m<N).');
    end
    
    % Useful parameters
    signal = signal(:)';
    N = length(signal);     % Signal length
    sigma = std(signal);    % Standard deviation
    
    % Create the matrix of matches
    matches = NaN(m+1,N);
    for i = 1:1:m+1
        matches(i,1:N+1-i) = signal(i:end);
    end
    matches = matches';
    % Check the matches for m
    d_m = pdist(matches(:,1:m), dist_type);
    if isempty(d_m)
        % If B = 0, SampEn is not defined: no regularity detected
        %   Note: Upper bound is returned
        value = Inf;
    else
        % Check the matches for m+1
        d_m1 = pdist(matches(:,1:m+1), dist_type);
        
        % Compute A and B
        %   Note: logical operations over NaN values are always 0
        B = sum(d_m  <= r*sigma);
        A = sum(d_m1 <= r*sigma);
        % Sample entropy value
        %   Note: norm. comes from [nchoosek(N-m+1,2)/nchoosek(N-m,2)]
        value = -log((A/B)*((N-m+1)/(N-m-1))); 
    end
    
    % If A=0 or B=0, SampEn would return an infinite value. However, the
    % lowest non-zero conditional probability that SampEn should
    % report is A/B = 2/[(N-m-1)(N-m)]
    if isinf(value)
        % Note: SampEn has the following limits:
        %       - Lower bound: 0
        %       - Upper bound: log(N-m)+log(N-m-1)-log(2)
        value = -log(2/((N-m-1)*(N-m)));
    end
end



function [ e, A, B ] = multiscaleSampleEntropy( x, m, r, tau )
%MULTISCALESAMPLEENTROPY
%
% Based on "Multiscale entropy analysis of biological signals"
% By Madalena Costa, Ary L. Goldberger, and C.-K. Peng
% Published on 18 February 2005 in Phys. Rev. E 71, 021906.
%
% This code was implemented by John Malik on 26 April 2017.
% Contact: john.malik@duke.edu
switch nargin
    case 1
        m = 2;
        r = 0.15;
        tau = 1;
    case 2
        r = 0.15;
        tau = 1;
    case 3
        tau = 1;
end
% coarse signal
y = mean(buffer(x(:), tau), 1);
% (m+1)-element sequences
X = buffer(y, m + 1, m, 'nodelay')';
% matching (m+1)-element sequences
A = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% matching m-element sequences
X = X(:, 1:m);
B = sum(pdist(X, 'chebychev') < r * nanstd(x, 1));
% take log
if A == 0 || B == 0
    e = NaN;
    return
end
e = log(B / A);
end


function [mobility,complexity] = HjorthParameters(xV)
% [mobility,complexity] = HjorthParameters(xV)
% HJORTHPARAMETERS computes the Hjorth parameters mobility and complexity.
% INPUTS:
% - xV          : The given scalar time series (vector of size n x 1).
% OUTPUTS
% - mobility
%========================================================================
%     <HjorthParameters.m>, v 1.0 2010/02/11 22:09:14  Kugiumtzis & Tsimpiris
%     This is part of the MATS-Toolkit http://eeganalysis.web.auth.gr/

%========================================================================
% Copyright (C) 2010 by Dimitris Kugiumtzis and Alkiviadis Tsimpiris 
%                       <dkugiu@gen.auth.gr>

%========================================================================
% Version: 1.0

% LICENSE:
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program. If not, see http://www.gnu.org/licenses/>.

%=========================================================================
% Reference : D. Kugiumtzis and A. Tsimpiris, "Measures of Analysis of Time Series (MATS): 
% 	          A Matlab  Toolkit for Computation of Multiple Measures on Time Series Data Bases",
%             Journal of Statistical Software, in press, 2010

% Link      : http://eeganalysis.web.auth.gr/
%========================================================================= 
n = length(xV);
dxV = diff([0;xV]);
ddxV = diff([0;dxV]);
mx2 = mean(xV.^2);
mdx2 = mean(dxV.^2);
mddx2 = mean(ddxV.^2);

mob = mdx2 / mx2;
complexity = sqrt(mddx2 / mdx2 - mob);
mobility = sqrt(mob);

function [HFD] = Higuchi_FD(serie, Kmax) 
narginchk(1, inf)
if nargin < 2
    Kmax = 282;
end
%{
Script for computing the Higuchi Fractal Dimension (HDF) of a signal.
INPUT:
    serie: is the temporal series that one wants to analyze by HDF. 
    It must be a row vector.
    Kmax: maximum number of sub-series composed from the original. To
    determine its values, we have followed the recommendation of Doyle et
    al at "Discriminating between elderly and young using a fractal 
    dimension analysis of centre of pressure". 
OUTPUT:
    HFD: the HFD of the temporal series.
PROJECT: Research Master in signal theory and bioengineering - University of Valladolid
DATE: 02/03/2014
AUTHOR: Jes�s Monge �lvarez
%}
%% Checking the ipunt parameters:
control = ~isempty(serie);
assert(control,'The user must introduce a series (first inpunt).');
control = ~isempty(Kmax);
assert(control,'The user must introduce the Kmax parameter (second inpunt).');
%% Processing:
% Composing of the sub-series:
N = length(serie); 
X = NaN(Kmax,Kmax,N);
for k = 1:Kmax
    for m = 1:k
        limit = floor((N-m)/k);
        j = 1;
        for i = m:k:(m + (limit*k))
            X(k,m,j) = serie(i);
            j = j + 1;
        end  
    end
end
% Computing the length of each sub-serie:
L = NaN(1, Kmax);
for k = 1:Kmax
    L_m = zeros(1,k);
    for m = 1:k
        R = (N - 1)/(floor((N - m)/k) * k);
        aux = squeeze(X(k,m,logical(~isnan(X(k,m,:))))); %We get the sub-serie without the NaNs.
        for i = 1:(length(aux) - 1)
            L_m(m) = L_m(m) + abs(aux(i+1) - aux(i));
        end
        L_m(m) = (L_m(m) * R)/k;
    end
    L(k) = sum(L_m)/k;
end
% Finally, we compute the HFD:
x = 1./(1:Kmax);
aux = polyfit(log(x),log(L),1);
HFD = aux(1); %We only want the slope, not the independent term. 