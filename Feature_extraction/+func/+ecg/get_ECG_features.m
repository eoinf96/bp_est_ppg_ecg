function [ECG_feats, feature_names] = get_ECG_features(ts, configs)
% This function returns the complexity features computes on an ECG segment ts
% Many methods used in this function have been taken from MathWorks
% FileExchage -- the authors of these methods have been highlighted
%
% INPUT: ts: ECG segment time series
%       configs: configs struct, see below for details
%
% OUTPUT: ECG_feats: Vector of feature values
%         feature_names: String cell of feature names
% ---
% Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure.
%
% Released under the GNU General Public License
%
% Copyright (C) 2022  Eoin Finnegan
% University of Oxford, Insitute of Biomedical Engineering, CIBIM Lab
% eoin.finnegan@eng.ox.ac.uk
%
% Referencing this work
%
% Finnegan, E., Davidson, S., Harford, M., Jorge, J., Watkinson, P., Tarassenko, L. and Villarroel, M., 2022. Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure. Submitted to Scientific reports
%
%
% Relevant literature:
% - Monika Simjanoska et al. “Non-invasive blood pressure estimation from ECG using machine learning techniques”. In: Sensors (Switzerland) 18.4 (2018), p. 1160.
% - Sen Yang et al. “Blood pressure estimation with complexity features from electrocardiogram and photoplethysmogram signals”. In: Optical and Quantum Electronics 52.3 (2020), p. 135.
% - 

narginchk(1, inf)
if nargin < 2
    configs = struct();
end
default_configs.MSE_scale_list = 2:2:8;
configs = func.aux_functions.update_with_default_opts(configs, default_configs);
%% Get complexity features
% Hjorth params
[mobility, complexity ] = HjorthParameters(ts);

% Fractal dimension
fractal = Higuchi_FD(ts', 17);

% Shannon entropy
p = histcounts(ts, 'Normalization', 'Probability');
shannon_entropy = -sum(p(p>0).*log2(p(p>0)) );


% Entropy features
m = 2;
r_factor = 0.2;%* std(ts);
approxEnt = approxEntropy(ts, m, r_factor );
sampEnt = sampleEntropy(ts, m, r_factor );

mseEnt = nan(length(configs.MSE_scale_list), 1);
for scale_idx = 1:length(configs.MSE_scale_list)
    mseEnt(scale_idx) = multiscaleSampleEntropy( ts, m, r_factor, configs.MSE_scale_list(scale_idx));
end

%% Combine features
ECG_feats = [mobility; complexity; fractal; shannon_entropy; approxEnt; sampEnt; mseEnt(:)];
feature_names = vertcat({'Hjorth_mobility'; 'Hjorth_complexity'; 'Fractal_dimension'; 'SE'; 'approxEnt'; 'sampEnt' }, strcat('MSE_scale_', cellstr(string(configs.MSE_scale_list)))');

end

%% Local functions -- all from MathWorks


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
function apen = approxEntropy(data,  dim, r, tau )
%ApEn
%   dim : embedded dimension
%   r : tolerance (typically 0.2 * std)
%   data : time-series data
%   tau : delay time for downsampling
%   Changes in version 1
%       Ver 0 had a minor error in the final step of calculating ApEn
%       because it took logarithm after summation of phi's.
%       In Ver 1, I restored the definition according to original paper's
%       definition, to be consistent with most of the work in the
%       literature. Note that this definition won't work for Sample
%       Entropy which doesn't count self-matching case, because the count 
%       can be zero and logarithm can fail.
%
%       A new parameter tau is added in the input argument list, so the users
%       can apply ApEn on downsampled data by skipping by tau. 
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Ver 0 : Aug 4th, 2011
% Ver 1 : Mar 21st, 2012
%---------------------------------------------------------------------
if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end
    
N = length(data);
result = zeros(1,2);
for j = 1:2
    m = dim+j-1;
    phi = zeros(1,N-m+1);
    dataMat = zeros(m,N-m+1);
    
    % setting up data matrix
    for i = 1:m
        dataMat(i,:) = data(i:N-m+i);
    end
    
    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        tempMat = abs(dataMat - repmat(dataMat(:,i),1,N-m+1));
        boolMat = any( (tempMat > r),1);
        phi(i) = sum(~boolMat)/(N-m+1);
    end
    
    % summing over the counts
    result(j) = sum(log(phi))/(N-m+1);
end
apen = result(1)-result(2);
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
end
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
end