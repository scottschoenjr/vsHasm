%**************************************************************************
%
% 3D Heterogeneous ASA Reconstruction for Verasonics
%
%   Function computes volumetric PAM from received signals on a matrix
%   array with a known SoS field. It should work with a 1D linear probe, 
%   but I have not tested this.
%
%   If used, please cite:
%   [1] Arvanitis et al. "Passive Acoustic Mapping with the Angular 
%       Spectrum Method" IEEE T. Med. Imag. 36(4) (2016)
%       DOI: 10.1109/TMI.2016.2643565
%   [2] Schoen Jr & Arvanitis "Heterogeneous Angular Spectrum Method for 
%       Trans-skull Imaging and Focusing" IEEE T. Med. Imag. 39(5) (2019)
%       DOI: 10.1109/TMI.2019.2953872
%   [3] Schoen Jr & Arvanitis "Acoustic Source Localization with the 
%       Angular Spectrum Approach in Continuously Stratified Media" 
%       J. Acoust. Soc. Am. 148 (EL333–EL339) (2020)
%       DOI: https://doi.org/10.1121/10.0002095
%   [4] Schoen Jr et al. "Experimental Demonstration of Trans-skull 
%       Volumetric Passive Acoustic Mapping with the Heterogeneous 
%       Angular Spectrum Approach" IEEE T. Ultrason. Ferro. 69(2) (2021)
%       DOI: 10.1109/TUFFC.2021.3125670
%
% Inputs
%   rfData  - Matrix of receive data from probe
%   Trans   - Verasonics Trans structure
%   Receive - Verasonics Receive Structure
%   ASAOptions - Struct with fields:
%     .fVector - Vector of frequencies [Hz]
%     .binWidth - Number of frequency bins to use
%     .zVector - Vector of depths [m]
%     .numPadChannels - Number of pad channels to add in each direcion
%     .tukeyWindowFraction - Fraction of physical aperture
%   props.
%     .cField - Sound speed at each point
%     .xVec - Vector of x positions
%     .yVec - Vector of y positions
%     .zVec - Vector of y positions
%
% Outputs
%   asamap  - Volumetric PAM
%   xVec   |
%   yVec   |- Pixel position arrays [m]
%   zVec   |
%   cmpTime - Total computation time
%
%
%        Scott Schoen Jr | scottschoenjr@gatech.edu | 2020-2021
%
%**************************************************************************

function [asamap, xVec, yVec, zVec, cmptime] = ...
    getVsHASA3D(rfData, Trans, Receive, ASAOptions )

% Get total number of sensors and number of time points
[Nt, numSensors] = size( rfData );

% Get the number of elements in each direction
Nx = length( unique( Trans.ElementPos( :, 1 ) ) );
Ny = length( unique( Trans.ElementPos( :, 2 ) ) );

% Store number of actual channels
Nx_txdr = Nx;
Ny_txdr = Ny;

% Reshape the data
rfData = double( rfData );
rfDataNew = zeros( Nx, Ny, Nt );

channelCount = 0;
for xCount = 1 : Nx
    for yCount = 1 : Ny
        
        % Increment
        channelCount = channelCount + 1;
        
        % Store data to grid
        rfDataNew( xCount, yCount, : ) = rfData( :, channelCount );
        
    end
end

% Window spatial data
try
    R = abs( ASAOptions.tukeyWindowFraction );
    if R > 1
        R = 0.25; % Default setting
    end
catch
    R = 0.25;
end

% Define 2D window
windowCols = tukeywin( Nx, R );
windowRows = tukeywin( Ny, R );
[wc, wr] = meshgrid( windowCols, windowRows );
window2D = wc.*wr;

% Get the element spacings
dx = 1E-3.*Trans.spacingMm;
dy = 1E-3.*Trans.spacingMm;

% If a number of pad channels hasn't been set, use number of elements in
% each direction
try
    numPadChannelsX = ASAOptions.numPadChannels(1);
    numPadChannelsY = ASAOptions.numPadChannels(2);
catch
    numPadChannelsX = Nx;
    numPadChannelsY = Ny;
end

% Pad RF data
rfDataPadded = padarray( rfDataNew, ...
    [numPadChannelsX, numPadChannelsY], 0, 'both' );
window2DPadded = padarray( window2D, ...
    [numPadChannelsX, numPadChannelsY], 0, 'both' );

% Update number of elements
[Nx, Ny, ~] = size( rfDataPadded );

% Get data into time frequency domain
pTilde = fft( rfDataPadded, [], 3 );

% Apply to spatial data
window2D_volume = repmat( window2DPadded, [1, 1, Nt] );
pTilde = pTilde.*window2D_volume;

% Define sensor positions in space
xSensors = ( 0 : dx : dx.*(Nx - 1) );
ySensors = ( 0 : dy : dy.*(Ny - 1) );

% Shift so centered at 0, 0
xShift = mean( xSensors );
yShift = mean( ySensors );

xSensors = xSensors - xShift;
ySensors = ySensors - yShift;

% [xPlot, yPlot] = meshgrid( xSensors, ySensors );

% Get z vector
zVector = ASAOptions.zVector;
Nz = length( zVector );
dz = zVector(2) - zVector(1);

% Get time and frequency vectors
Fs = 1E6.*Receive.decimSampleRate;
dt = 1./Fs;

tVec = 0 : dt : ( dt.*(Nt-1) ); % Time vector
fVec = linspace(0, Fs, length(tVec)); % Frequency vector

% Find indices of frequency vector to use for reconstruction
fIndices = [];
freqsToUse = ASAOptions.fVector;
fbinWidth = ASAOptions.fbinWidth;
for fCount = 1 : length( freqsToUse )
    
    % Find closest index
    [~, fIndex] = min( abs( fVec - freqsToUse(fCount) ) );
    
    % Get frequencies in this bin
    halfBinWidth = round( (fbinWidth - 1) / 2 );
    startIndex = fIndex - halfBinWidth;
    endIndex = fIndex + halfBinWidth;
    
    % Make sure in bounds
    startIndex = max( startIndex, 1 );
    endIndex = min( endIndex, length( fVec ) );
    
    % Append to list of indices to use in reconstruction (and make sure we
    % don't use any bins twice).
    fIndices = unique( [ fIndices, (startIndex : endIndex) ] );
    
end


% Interpolate pTilde at the required frequencies
try
    interpFactor = ASAOptions.interpFactor;
catch
    interpFactor = 1;
end

if interpFactor > 1
    
    % Store orinal grid points
    x0 = xSensors;
    y0 = ySensors;
    [x0, y0] = meshgrid( xSensors, ySensors );
    
    % Get new
    Nx = round( interpFactor.*Nx );
    Ny = round( interpFactor.*Ny );
    
    dx = dx./interpFactor;
    dy = dy./interpFactor;
    
    % Define sensor positions in space
    xSensors = ( 0 : dx : dx.*(Nx - 1) );
    ySensors = ( 0 : dy : dy.*(Ny - 1) );
    
    % Shift so centered at 0, 0
    xShift = mean( xSensors );
    yShift = mean( ySensors );
    
    xSensors = xSensors - xShift;
    ySensors = ySensors - yShift;
    
    % Get interpolation points
    [xi, yi] = meshgrid( xSensors, ySensors );
    
    % Interpolate each slice and store to final array
    for fIndCount = 1 : length(fIndices)
        
        pTildeLoop = squeeze( pTilde( :, :, fIndices(fIndCount) ) );
        pTildeLoopInterp = interp2( x0, y0, pTildeLoop, xi, yi );
        
        % Make sure no NaNs
        pTildeLoopInterp(isnan( pTildeLoopInterp) ) = 0;
        
        % Store to full array
        pTildeSlices( :, :, fIndCount ) = pTildeLoopInterp;
        
    end
    
else
    
    % Get pTilde just at the necessary frequencies
    pTildeSlices = pTilde( :, :, fIndices );
    
end

% Interpolate sound speed to coputational grid

% First check if we need to pad
xMax = max(abs(xSensors(:)));
yMax = max(abs(ySensors(:)));
zMax = max(abs(zVector(:)));

% Get material properties
props = ASAOptions.props;
c = props.cField;

% Pad if needed
% x
if xMax > max( props.xVec )
    
    % Get number of pixels to add in x direction
    NxProps = length( props.xVec );
    dxProps = props.xVec(2) - props.xVec(1);
    extraPixelsX = ceil( (max(xMax) - max(props.xVec))./dxProps );
    
    % Add on either end (centered at 0 still)
    NxPadded = NxProps + 2.*extraPixelsX;
    xProps = ( 0 : dxProps : dxProps.*(NxPadded - 1) );
    xProps = xProps - mean(xProps);
    
    % Pad sound speed array
    c = padarray( c, [0, extraPixelsX], 'replicate', 'both' );
    
else
    xProps = props.xVec;
end
% y
if yMax > max( props.yVec )
    
    % Get number of pixels to add in y direction
    NyProps = length( props.yVec );
    dyProps = props.yVec(2) - props.yVec(1);
    extraPixelsY = ceil( (max(yMax) - max(props.yVec))./dyProps );
    
    % Add on either end (centered at 0 still)
    NyPadded = NyProps + 2.*extraPixelsY;
    yProps = ( 0 : dyProps : dyProps.*(NyPadded - 1) );
    yProps = yProps - mean(yProps);
    
    c = padarray( c, [extraPixelsY, 0], 'replicate', 'both' );
    
else
    yProps = props.yVec;
    
end
% z
if zMax > max( props.zVec )
    
    % Get number of pixels to add in y direction
    NzProps = length( props.zVec );
    dzProps = props.zVec(2) - props.zVec(1);
    extraPixelsZ = ceil( (max(zMax) - max(props.zVec))./dzProps );
    
    % Add to end in z
    zAppend = (dzProps : dzProps : dzProps.*extraPixelsZ) + max(props.zVec);
    zProps = [ props.zVec, zAppend ];
    
    c = padarray( c, [0, 0, extraPixelsZ], 'replicate', 'post' );
    
else
    zProps = props.zVec;
end


% Interpolate
[xr, yr, zr] = meshgrid( xSensors, ySensors, zVector ); % Reconstruction
[xc, yc, zc] = meshgrid( xProps, yProps, zProps ); % Data
c = interp3( xc, yc, zc, c, xr, yr, zr );

% Get sound speed to use
c0 = mean( c(:) );

% ASA Reconstruction ------------------------------------------------------
tic;

% Iniatialize ASA map
asamap = zeros(Nx, Ny, Nz);

% Compute x and y wavenumbers
dkx = 2.*pi./dx;
dky = 2.*pi./dy;

kxVec = (( (0 : Nx - 1) - Nx./2 )./Nx).*dkx;
kyVec = (( (0 : Ny - 1) - Ny./2 )./Ny).*dky;

% Arrange into transverse wavenumber grid
[kx, ky] = meshgrid( kxVec, kyVec );

% Repeat z for matrix calculation
% TODO: Probably can be done with repmat for speed, but just for now...
z_vol = zeros( Nx, Ny, Nz );
for xCount = 1 : Nx
    for yCount = 1 : Ny
        z_vol( xCount, yCount, : ) = zVector;
    end
end


% Now compute for each frequency
for fCount = 1 : length( fIndices )
    
    % Get current index
    fInd = fIndices( fCount );
    
    % Get the rfData at that frequency
    pTilde0 = squeeze( pTildeSlices( :, :, fCount ) ); % [Nx-by-Ny]
    
    % Perform 2D FFT to get angular spectrum at receiver for this
    % (time) frequency
    P0 = fftshift( fft2( pTilde0 ) );
    
    % Get the propagating wavenumber for this frequency
    omega0 = 2.*pi.*fVec( fInd );
    k0 = omega0./c0;
    kz = sqrt( (omega0./c0).^(2) - kx.^(2) - ky.^(2) );
    
    % Initialize P
    Pn = P0;
    P = repmat( P0, 1, 1, Nz );
    
    % Get cutoff frequency
    cutoffThreshold = 0.0;
    realInds = find( real( kz./k0 ) > cutoffThreshold );
    
    % Define window function
    xLocal = real( kx./k0 );
    yLocal = real( ky./k0 );
    rLocal = sqrt( xLocal.^(2) + yLocal.^(2) );
    R = 0.25;
    window = ...
        (rLocal < (1 - R) ).*( ...
        1 ...
        ) ...
        + ...
        ( rLocal >= (1 - R) & rLocal < 1).*( ...
        0.5.*( 1 + cos( (2.*pi./(2.*R)).*( rLocal - R ) ) ) ...
        ) ...
        + ...
        ( rLocal >= 1).*( ...
        0 ...
        );
    
    % Define weighting window
weighting = window;
    
    % Define lambda for this frequency
    mu = (c0./c).^(2);
    lambda = (omega0.^(2)./c0.^(2)).*( 1 - mu );
    
    % Perform marching algorithm
    for zCount = 2 : length(zVector)
        
        % Perform convolution
        pn = ifft2( ifftshift( Pn ) ); % Spatial domain pressure
        LambdaStarPn = fftshift( fft2( lambda(:, :, zCount).*pn ) );
        
        % Compute propagation step
        Pn1 = Pn.*exp(1i.*kz.*dz) + ...
            exp(1i.*kz.*dz)./(2i.*kz).*LambdaStarPn.*dz;
               
        % Update previous value and store
        P( :, :, zCount ) = Pn1;
        Pn = weighting.*Pn1;
        
    end
    
    % KLUDGE There's probably a more efficient way to do this, but just
    % gonna loop through z for now
    pTilde_vol = zeros( Nx, Ny, Nz );
    for zCount = 1 : Nz
        % Get 2D inverse transform of each z-slice.
        pTilde_vol( :, :, zCount ) = ...
            ifft2( ifftshift( P(:, :, zCount ) ) );
    end
    
    
    % Add into contribution
    asamap = asamap + abs( pTilde_vol ).^(2);   
    
end

% Return variables
cmptime = toc;
[xVec, yVec, zVec] = meshgrid( xSensors, ySensors, zVector );


end