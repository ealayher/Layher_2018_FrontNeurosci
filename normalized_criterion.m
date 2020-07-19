function [nc, resid, slope, int, r, r2] = normalized_criterion(cVec, dVec)
% normalized_criterion.m: Residualize c against d' to get normalized c values
%
% function [nc, resid, slope, int, r, r2] = normalized_criterion(cVec, dVec)
%
% Created: 04/10/2018, Evan Layher
% Revised: 05/02/2018, Evan Layher % Added residuals output
% Revised: 02/16/2020, Evan Layher % Fixed cSq calculation & minor updates
%
% --- LICENSE INFORMATION --- %
% Modified BSD-2 License - for Non-Commercial Use Only
%
% Copyright (c) 2018-20, The Regents of the University of California
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted for non-commercial use only provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list
%    of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice, this list
%    of conditions and the following disclaimer in the documentation and/or other
%    materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
% WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% For permission to use for commercial purposes, please contact UCSB's Office of
% Technology & Industry Alliances at 805-893-5180 or info@tia.ucsb.edu.
% --------------------------- %
%
% Signal detection theory: Adjust criterion (c) values to account for relationship with d'
% Use linear regression to find relationship between c and d'
% The regression residuals are added to mean c to get normalized c values
%
% Inputs: 2 vectors of equal length
% (1) cVec: vector of c values
% (2) dVec: vector of d' values
%
% Outputs:
% (1) nc: normalized c values (residuals + c mean)
% (2) resid: residuals
% (3) slope: regression slope
% (4) int: regression intercept
% (5) r: Pearson's r
% (6) r2: Pearson's r squared

if length(cVec) ~= length(dVec) % Crash if vector lengths are unequal
    error('VECTOR LENGTHS MUST BE EQUAL: "cVec" (%d) AND "dVec" (%d)\n', ...
        length(cVec), length(dVec))
end

mCs = mean(cVec); % mean c (y)
mDs = mean(dVec); % mean d' (x)

cSq = sum((cVec - mCs) .^ 2); % sum of squares: c (y) ** FIXED: 02/16/2020 **
dSq = sum((dVec - mDs) .^ 2); % sum of squares: d' (x)
cdSq = sum((dVec - mDs) .* (cVec - mCs)); % sum of squares: c x d' (y)

slope = cdSq / dSq; % regression slope
int = mCs - (slope * mDs); % regression intercept
r = cdSq / sqrt((dSq * cSq)); % Pearson's r
r2 = r ^ 2; % explained variance (%)
ln = slope .* dVec + int; % y = mx + b
resid = (cVec - ln); % residuals
nc = resid + mCs; % normalized c (residuals + mean)

end