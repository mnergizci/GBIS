function [] = GBISrun(inputFileName, insarDataCode, gpsDataFlag, modelCode, nRuns, skipSimulatedAnnealing)

%%  Geodetic Bayesian Inversion Software (GBIS)
%   Software for the Bayesian inversion of geodetic data.
%   Copyright: Marco Bagnardi, 2018
%
%%  =========================================================================
%   Usage: GBISrun(inputFileName, insarDataCode, gpsDataFlag, modelCode, nRuns, skipSimulatedAnnealing)
%
%   inputFileName:  name and extension of input file (e.g., 'VolcanoExercise.inp')
%
%   insarDataCode:  select data to use (e.g., [1,3] to use InSAR data
%                   with insarID = 1 and 3; insarID specified in input file *.inp).
%                   Leave empty (e.g.,[]) if no InSAR data.
%
%   gpsDataFlag:    'y' to use GPS data, 'n' to not use GPS data
%
%   modelCode:      select forward models to use;
%                   'M' for Mogi source (point source, [Mogi, 1958])
%                   'T' for McTigue source (finite sphere, [McTigue, 1987])
%                   'Y' for Yang source (prolate spheroid, [Yang et al., 1988])
%                   'S' for rectangular horizontal sill with uniform opening [Okada, 1985]
%                   'P' for penny-shaped crack [Fialko et al., 2001]
%                   'D' for rectangular dipping dike with uniform opening [Okada, 1985]
%                   'F' for rectangular dipping fault with uniform slip [Okada, 1985]
%
%                   Custom made sources:
%                   'H' for two hinged rectangular dikes (hinged at the bottom edge of the upper dike)
%
%   nRuns:          number of iterations (samples) to be performed (e.g., 1000000)
%
%   skipSimulatedAnnealing: 'y' to skip initial Simulated Annealing phase
%
%   Example: >> GBISrun('VolcanoExercise.inp',[1,2],'y','MD',1e06, 'n')
%   =========================================================================
%
%   Email: gbis.software@gmail.com
%
%   Reference: 
%   Bagnardi M. & Hooper A, (2018). 
%   Inversion of surface deformation data for rapid estimates of source 
%   parameters and uncertainties: A Bayesian approach. Geochemistry, 
%   Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585
%
%   Last update: 8 August, 2018

%% Check number of input arguments and return error if not sufficient

if nargin == 0
    help GBISrun;
    return;
end

if nargin < 6
    disp('#############################################')
    disp('####### Not enough input arguments. #########')
    disp('#############################################')
    disp(' ')
    disp('Type "help GBISrun" for more information');
    disp(' ')
    return;
end

%% Start timer and initialise global variables
clc % Clean screen
tic % Start timer

clear  outputDir  % Clear global variables
global outputDir  % Set global variables

%% Diplay header
disp('**********************************************')
disp('Geodetic Bayesian Inversion Software (GBIS)')
disp('Software for the Bayesian inversion of geodetic data.')
disp('Markov chain Monte Carlo algorithm incorporating the Metropolis algorithm')
disp(' ')
disp('by Marco Bagnardi and Andrew Hooper (COMET, University of Leeds)')
disp('Emal: M.Bagnardi@leeds.ac.uk')
disp('Last update: 12/05/2017')
disp('**********************************************')
disp(' ')
disp(' ')

%% Read input data from input text file *.inp

inputFileID = fopen(inputFileName, 'r');
textLine = fgetl(inputFileID); 

while ischar(textLine)
    eval(textLine)
    textLine = fgetl(inputFileID);
end

fclose(inputFileID);

%% Create output directories and output file name

[inputFile.path, inputFile.name, inputFile.ext] = fileparts(inputFileName); % Extract name
outputDir = inputFile.name; % Name output directory as input file name

   % InSAR + GPS
if gpsDataFlag == 'y' && ~isempty(insarDataCode)
    
    % Add InSAR dataset ID
    for i = 1:length(insarDataCode)
        insarDataNames{i} = num2str(insarDataCode(i));
    end
    
    % Add GPS
    saveName = ['invert_',strjoin(insarDataNames, '_'),'_GPS'];
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
    
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    
    % Add .mat extension
    saveName = [saveName,'.mat'];

    % InSAR only
elseif gpsDataFlag == 'n' && ~isempty(insarDataCode)
    
    % Add InSAR dataset ID
    for i = 1:length(insarDataCode)
        insarDataNames{i} = num2str(insarDataCode(i));
    end

    saveName = ['invert_',strjoin(insarDataNames, '_')];
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
    
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    
    % Add .mat extension
    saveName = [saveName,'.mat'];
    
    % GPS only
elseif gpsDataFlag == 'y' && isempty(insarDataCode)
    
    % Add GPS
    saveName = ['invert_GPS'];
    
    % Add models
    for i = 1:length(modelCode)
        saveName = [saveName,'_',modelCode(i)];
    end
            
    % Create output directories
    outputDir = [outputDir,'/',saveName];
    disp(['Output directory: ', outputDir])
    mkdir(outputDir)
    mkdir([outputDir,'/Figures']) % Create directory for Figures
    
    % Add .mat extension
    saveName = [saveName,'.mat'];
end

%% Initialise variables

nObs = 0;   % Initialise number of observations variable
obs  = [];  % Initialise observation points (x,y,z) matrix

%% Ingest InSAR data

% Create colormaps for plotting InSAR data (call third party colormap_cpt.m function and GMT *.cpt files)
cmap.seismo = colormap_cpt('GMT_seis.cpt', 100);    % GMT 'Seismo' colormap for wrapped data
cmap.redToBlue = colormap_cpt('polar.cpt', 100);    % Red to Blue colormap for unwrapped data

% Select InSAR datasets to use from list in input file
if ~isempty(insarDataCode)
    disp('InSAR datasets used in this inversion:')
    
    for i = 1:length(insarDataCode)
        selectedInsarData{i} = insar{insarDataCode(i)};
        disp(insar{insarDataCode(i)}.dataPath)    % display filename of InSAR dataset
    end
    
    % Create subsampled InSAR datasets for inversion
    disp(' ')
    disp('Ingesting InSAR data and performing Quadtree subsampling ...')
    [insar, obsInsar, nObsInsar] = loadInsarData(selectedInsarData, geo, cmap); % Load and subsample InSAR data
    nObs = nObsInsar;   % Add number of InSAR data points to total number of points
    obs  = obsInsar;    % Add InSAR observation points to observation points
else
    disp(' ')
    disp 'No InSAR datasets will be used in this inversion.'
    insar = [];
end

%% Ingest GPS data

if gpsDataFlag == 'y'
    disp(' ')
    disp('GPS dataset used in this inversion:')
    disp(gps.dataPath) % display filename of GPS data file
    disp(' ')
    disp('Ingesting GPS data ...')
    [gps, obsGps, nObsGps] = loadGpsData(gps, geo);
    disp([num2str(nObsGps), ' GPS sites will be used in this inversion'])
    gps.ix = [nObs+1:nObs+nObsGps]; % Index of GPS data in observation vector
    nObs = [nObs + nObsGps];   % Add number of GPS data points to total number of points
    obs = [obs, obsGps];       % Add GPS observation points to observation points
else
    disp(' ')
    disp 'No GPS datasets will be used in this inversion.'
    gps = [];
end

%% Plug-in here any further type of data to ingest (i.e., differential DEMs)

%% Setup inversion parameters

% Pause and press key to continue
disp ' '
disp '#################   Press any key to continue ...'
pause

% Define inversion parameters
disp 'Preparing for inversion ...'

invpar.nSave = 1000;    % Save output to file every 1000 iterations (every 10,000 after 20,000 iterations)
invpar.sensitivitySchedule = [1:100:10000,11000:1000:30000,40000:10000:nRuns]; % sensitivity schedule (when to change step sizes)

if skipSimulatedAnnealing == 'y'
    invpar.TSchedule = 1; % No temperature schedule if Simulated Annealing is not performed
else
    invpar.TSchedule = 10.^(3:-0.2:0); % Cooling schedule for Simulated Annealing
end
invpar.TRuns = 1000; % Number of runs at each temperature (Simulated Annealing only)

invpar.nModels = length(modelCode); % number of models used (e.g., 'MD' = 2x models)

invpar.nRuns = nRuns;

% Switch model code to full model name
% ADD HERE ANY NEW CUSTOMISED MODEL

for i = 1:invpar.nModels
    modelName = modelCode(i);
    switch modelName
        case 'M'
            invpar.model{i}='MOGI';    % Mogi source
        case 'T'
            invpar.model{i}='MCTG';    % McTigue source
        case 'P'
            invpar.model{i}='PENN';   % Penny-shaped crack (Fialko 2001)
        case 'Y'
            invpar.model{i}='YANG';    % Yang source
        case 'S'
            invpar.model{i}='SILL';    % Sill simulated as horizontal dislocation (Okada)
        case 'D'
            invpar.model{i}='DIKE';    % Dipping dike dislocation (Okada)
        case 'F'
            invpar.model{i}='FAUL';    % Dipping fault dislocation (Okada)
            
            % Custum made models
        case 'H'
            invpar.model{i}='HING';    % Two dikes hinged along L at depth (2x Okada)
            
            % End of custum made models
            
        otherwise
            error 'Invalid model code.'
    end
end

model = prepareModel(modelInput, invpar, insar, gps);

%% Run inversion

invResults = runInversion(geo, gps, insar, invpar, model, modelInput, obs, nObs);

%% Create *.mat file with final results

if gpsDataFlag == 'y' && ~isempty(insar)

    cd(outputDir)
    save(saveName, 'insarDataCode', 'geo', 'inputFile', 'invpar', 'gps', 'insar', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName')
    delete temporary.mat
    cd ../..
    
elseif gpsDataFlag == 'n' && ~isempty(insar)

    cd(outputDir)
    save(saveName, 'insarDataCode', 'geo', 'inputFile', 'invpar', 'insar', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName')
    delete temporary.mat
    cd ../..
    
elseif gpsDataFlag == 'y' && isempty(insar)

    cd(outputDir)
    save(saveName, 'geo', 'inputFile', 'invpar', 'gps', 'model', 'modelInput', 'invResults', 'obs', 'nObs', 'saveName')
    delete temporary.mat
    cd ../..
end

%% Display inversion duration
disp('=========================================================')
disp(['Time since start (HH:MM:SS):  ',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')])