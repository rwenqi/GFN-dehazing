function [peaksnr, snr] = psnr(A, ref, peakval)
%PSNR Peak Signal-To-Noise Ratio.
%   PEAKSNR = PSNR(A, REF) calculates the peak signal-to-noise ratio for
%   the image in array A, with the image in array REF as the reference. A
%   and REF can be N-D arrays, and must be of the same size and class.
% 
%   PEAKSNR = PSNR(A, REF, PEAKVAL) uses PEAKVAL as the peak signal value
%   for calculating the peak signal-to-noise ratio.
% 
%   [PEAKSNR, SNR] = PSNR(A, REF, __) also returns the simple
%   signal-to-noise in SNR, in addition to the peak signal-to-noise ratio.
%
%   Example
%   ---------
%   This example shows how to compute PSNR for noisy image given the
%   original reference image.
% 
%   ref = imread('pout.tif');
%   A = imnoise(ref,'salt & pepper', 0.02);
% 
%   [peaksnr, snr] = psnr(A, ref);
% 
%   fprintf('\n The Peak-SNR value is %0.4f', peaksnr);
%   fprintf('\n The SNR value is %0.4f \n', snr);
% 
%   Class Support
%   -------------
%   Input arrays A and REF must be one of the following classes: uint8,
%   int16, uint16, single, or double. Both A and REF must be of the same
%   class. They must be nonsparse. PEAKVAL is a scalar of any numeric
%   class. PEAKSNR and SNR are scalars of class double, unless A and REF
%   are of class single in which case PEAKSNR and SNR are scalars of class
%   single.
%
%   See also IMMSE, MEAN, MEDIAN, SSIM, SUM, VAR.

%   Copyright 2013-2015 The MathWorks, Inc. 

checkImages(A,ref);

if nargin < 3
    peakval = diff(getrangefromclass(A));
else
    checkPeakval(peakval, A);
    peakval = double(peakval);
end

if isempty(A) % If A is empty, ref must also be empty
    peaksnr = zeros(0, 'like', A);
    snr     = zeros(0, 'like', A);
    return;
end

err = immse(A,ref);

peaksnr = 10*log10(peakval.^2/err);

if nargout > 1
    if isinteger(ref)  
        ref = double(ref);
    end                      
    snr = 10*log10(mean(ref(:).^2)/err);
end

end

function checkImages(A, ref)

validImageTypes = {'uint8','uint16','int16','single','double'};

validateattributes(A,validImageTypes,{'nonsparse'},mfilename,'A',1);
validateattributes(ref,validImageTypes,{'nonsparse'},mfilename,'REF',2);

if ~isa(A,class(ref))
    error(message('images:validate:differentClassMatrices','A','REF'));
end   
if ~isequal(size(A),size(ref))
    error(message('images:validate:unequalSizeMatrices','A','REF'));
end

end

function checkPeakval(peakval, A)

validateattributes(peakval,{'numeric'},{'nonnan', 'real', ...
    'nonnegative','nonsparse','nonempty','scalar'}, mfilename, ...
    'PEAKVAL',3);
if isinteger(A) && (peakval > diff(getrangefromclass(A)))
    warning(message('images:psnr:peakvalTooLarge', 'A', 'REF'));
end

end
