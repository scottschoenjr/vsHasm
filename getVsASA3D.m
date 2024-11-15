%**************************************************************************
%
% 3D Homogeneous ASA Reconstruction for Verasonics
%
%   Function computes volumetric PAM from received signals on a matrix
%   array. It should work with a 1D linear probe, but I have not tested
%   this.
%
%   If used, please cite:
%   [1] Arvanitis et al. "Passive Acoustic Mapping with the Angular 
%       Spectrum Method" IEEE T. Med. Imaging 36(4) (2016)
%       DOI: 10.1109/TMI.2016.2643565
%   [2] Schoen Jr et al. "Experimental Demonstration of Trans-skull 
%       Volumetric Passive Acoustic Mapping with the Heterogeneous 
%       Angular Spectrum Approach" IEEE T. Ultrason. Ferro. 69(2) (2021)
%       DOI: 10.1109/TUFFC.2021.3125670
%
% Inputs
%   rfData  - Matrix of receive data from probe
%   Trans   - Verasonics Trans structure
%   Receive - Verasonics Receive Structure
%   ASAOptions - Struct with fields:
%      .fVector        - Vector of frequencies [Hz]
%      .binWidth       - Number of frequency bins to use
%      .zVector        - Vector of depths [m]
%      .c0             - Sound Speed [m/s]
%      .numPadChannels - Number of spatial padding channels
%      .interpFactor   - Spatial interpolation factor
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
%--------------------------------------------------------------------------

function [asamap, xVec, yVec, zVec, cmptime] = ...
    getVsASA3D(rfData, Trans, Receive, ASAOptions )

% Get total number of sensors and number of time points
[Nt, numSensors] = size( rfData );

% Get the number of elements in each direction
Nx = length( unique( Trans.ElementPos( :, 1 ) ) );
Ny = length( unique( Trans.ElementPos( :, 2 ) ) );

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

% Update number of elements
[Nx, Ny, ~] = size( rfDataPadded );

% Get data into time frequency domain
pTilde = fft( rfDataPadded, [], 3 );

% Define sensor positions in space
xSensors = ( 0 : dx : dx.*(Nx - 1) );
ySensors = ( 0 : dy : dy.*(Ny - 1) );

% Shift so centered at 0, 0
xShift = mean( xSensors );
yShift = mean( ySensors );

xSensors = xSensors - xShift;
ySensors = ySensors - yShift;

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


% Get sound speed to use
c0 = ASAOptions.soundSpeed;

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
kx = repmat( kxVec, Ny, 1 );
ky = repmat( kyVec, Nx, 1 );

[kx, ky] = meshgrid( kxVec, kyVec );

% Repeat z for matrix calculation
% TODO: Probably can be done with repmat, but just want to get working
%       for now...
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
    kz = sqrt( (omega0./c0).^(2) - kx.^(2) - ky.^(2) );
    
    % Repeat matrices over length of z so we can do the propagation in a
    % single step
    P0_vol = repmat( P0, 1, 1, Nz );
    kz_vol = repmat( kz, 1, 1, Nz );
    
    % Perform reconstrcution
    P_vol = P0_vol.*exp( 1i.*kz_vol.*z_vol );
    
    % Transform back to pressure
    %     pTilde_vol = ifft2( ifftshift( P_vol ) );
    
    % KLUDGE There's probably a more efficient way to do this, but just 
    % gonna loop through z for now
    pTilde_vol = zeros( Nx, Ny, Nz );
    for zCount = 1 : Nz
        % Get 2D inverse transform of each z-slice.
        pTilde_vol( :, :, zCount ) = ...
            ifft2( ifftshift( P_vol(:, :, zCount ) ) );
    end  
    
    % Add into contribution
    asamap = asamap + abs( pTilde_vol ).^(2);  
    
end

% Return variables
cmptime = toc;
[xVec, yVec, zVec] = meshgrid( xSensors, ySensors, zVector );

end