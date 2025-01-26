% Runs a simulation of a user-specified pulse sequence for each 
% orientation found in the user-specified grid and returns the 
% corresponding angles, integration weights, tessellation trian-
% gles and sequence results (as a cell array of whatever it is 
% that the pulse sequence returns). 
%
% The function supports parallel processing via Matlab's Distri-
% buted Computing Toolbox - different system orientations are eva-
% luated on different labs. Syntax:
%
% [alphas,betas,gammas,weights,triangles,fid]=...
%        roadmap(spin_system,pulse_sequence,parameters,assumptions)
%
% Input arguments:
%
%    pulse_sequence  -  pulse sequence function handle. See the
%                       /experiments directory for the list of
%                       pulse sequences that are supplied with
%                       Spinach.
% 
%    parameters      -  a structure containing sequence-specific
%                       parameters. See the pulse sequence header
%                       for the list of parameters each sequence
%                       requires. The parameters that this inter-
%                       face itself requires are:
%                      
%                       parameters.spins - a cell array giving the
%                       spins that the pulse sequence works on, in
%                       the order of channels, e.g. {'1H','13C'}
%
%                       parameters.offsets - a cell array giving
%                       transmitter offsets on each of the spins
%                       listed in parameters.spins.
%
%                       parameters.grid - name of the spherical 
%                       averaging grid file (see the grids direc-
%                       tory in the kernel).
%
%    assumptions     -  context-specific assumptions ('nmr', 'epr',
%                       'labframe', etc.) - see the pulse sequence
%                       header and assumptions.m file for further
%                       information on this setting.
%
% Output arguments:
%
%            alphas  -  an array of alpha Euler angles from the 
%                       spherical grid
%
%            betas   -  an array of beta Euler angles from the 
%                       spherical grid
%
%            gammas  -  an array of gamma Euler angles from the 
%                       spherical grid
%
%           weights  -  an array of summation weights appropri-
%                       ate for the corresponding points of the
%                       spherical grid
%
%         triangles  -  an array of triangle specification from
%                       the tesselation information supplied in
%                       the grid file
%
%               fid  -  a cell array of whatever it is that the
%                       pulse sequence function returns at each
%                       point inthe spherical grid
%
% i.kuprov@soton.ac.uk

function [alphas,betas,gammas,weights,triangles,answer]=roadmap(spin_system,pulse_sequence,parameters,assumptions)

% Show the banner
banner(spin_system,'sequence_banner'); 

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,pulse_sequence,parameters,assumptions);

% Report to the user
report(spin_system,'building the Liouvillian...');

% Set the secularity assumptions
spin_system=assume(spin_system,assumptions);

% Get the rotational expansion of the Hamiltonian
[I,Q]=hamiltonian(spin_system);

% Get kinetics
K=kinetics(spin_system);

% Apply the offsets
I=frqoffset(spin_system,I,parameters);

% Get carrier operators
C=cell(size(parameters.rframes));
for n=1:numel(parameters.rframes)
    C{n}=carrier(spin_system,parameters.rframes{n}{1});
end

% Get problem dimensions
parameters.spc_dim=1; parameters.spn_dim=size(I,1);

% Get the spherical averaging grid
spherical_grid=load([spin_system.sys.root_dir '/kernel/grids/' parameters.grid],...
                    'alphas','betas','gammas','weights','triangles'); 

% Assign local variables
alphas=spherical_grid.alphas; betas=spherical_grid.betas; 
gammas=spherical_grid.gammas; weights=spherical_grid.weights;
triangles=spherical_grid.triangles;
                           
% Preallocate the answer
answer=cell(numel(weights),1);

% Shut up and inform the user
report(spin_system,['powder average being computed over ' num2str(numel(weights)) ' orientations.']);
if ~isfield(parameters,'verbose')||(parameters.verbose==0)
    report(spin_system,'pulse sequence has been silenced to avoid excessive output.')
    spin_system.sys.output='hush';
end

% Crunch the orientations in parallel
parfor n=1:numel(weights)
    
    % Get the full Liouvillian at the current orientation
    H=I+orientation(Q,[alphas(n) betas(n) gammas(n)]); H=(H+H')/2;
    
    % Apply rotating frames
    for k=1:numel(parameters.rframes) %#ok<PFBNS>
        H=rotframe(spin_system,C{k},H,parameters.rframes{k}{1},parameters.rframes{k}{2}); %#ok<PFBNS>
    end
    
    % Get the relaxation superoperator at the current orientation
    R=relaxation(spin_system,[alphas(n) betas(n) gammas(n)]);

    % Report to the user
    report(spin_system,'running the pulse sequence...');
    
    % Run the simulation
    answer{n}=pulse_sequence(spin_system,parameters,H,R,K); %#ok<PFBNS>
    
end
    
end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'decouple')
    report(spin_system,'parameters.decouple field not set, assuming no decoupling.');
    parameters.decouple={};
end
if ~isfield(parameters,'rframes')
    report(spin_system,'parameters.rframes field not set, assuming no additional rotating frame transformations.');
    parameters.rframes={};
end
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'verbose')
    report(spin_system,'parameters.verbose field not set, silencing array operations.');
    parameters.verbose=0;
end
end

% Consistency checking
function grumble(spin_system,pulse_sequence,parameters,assumptions)

% Spherical grid
if ~isfield(parameters,'grid')
    error('spherical averaging grid must be specified in parameters.grid variable.');
elseif isempty(parameters.grid)
    error('parameters.grid variable cannot be empty.');
elseif ~ischar(parameters.grid)
    error('parameters.grid variable must be a character string.');
end

% Pulse sequence
if ~isa(pulse_sequence,'function_handle')
    error('pulse_sequence argument must be a function handle.');
end

% Assumptions
if ~ischar(assumptions)
    error('assumptions argument must be a character string.');
end

% Active spins
if isempty(parameters.spins)
    error('parameters.spins variable cannot be empty.');
elseif ~iscell(parameters.spins)
    error('parameters.spins variable must be a cell array.');
elseif ~all(cellfun(@ischar,parameters.spins))
    error('all elements of parameters.spins cell array must be strings.');
elseif any(~ismember(parameters.spins,spin_system.comp.isotopes))
    error('parameters.spins refers to a spin that is not present in the system.');
end

% Offsets
if isempty(parameters.offset)
    error('parameters.offset variable cannot be empty.');
elseif ~isnumeric(parameters.offset)
    error('parameters.offset variable must be an array of real numbers.');
elseif ~isfield(parameters,'spins')
    error('parameters.spins variable must be specified alongside parameters.offset variable.');
elseif numel(parameters.offset)~=numel(parameters.spins)
    error('parameters.offset variable must have the same number of elements as parameters.spins.');
end

% Rotating frame transformations
if ~isfield(parameters,'rframes')
    error('parameters.rframes variable must be specified.');
elseif ~iscell(parameters.rframes)
    error('parameters.rframes must be a cell array.');
end
for n=1:numel(parameters.rframes)
    if ~iscell(parameters.rframes{n})
        error('elements of parameters.rframes must be cell arrays.');
    end
    if numel(parameters.rframes{n})~=2
        error('elements of parameters.rframes must have exactly two sub-elements each.');
    end
    if ~ischar(parameters.rframes{n}{1})
        error('the first part of each element of parameters.rframes must be a character string.');
    end
    if ~ismember(parameters.rframes{n}{1},spin_system.comp.isotopes)
        error('parameters.rframes refers to a spin that is not present in the system.');
    end
    if ~isnumeric(parameters.rframes{n}{2})
        error('the second part of each element of parameters.rframes must be a number.');
    end
end

end

% The fact that we live at the bottom of a deep gravity well, on the
% surface of a gas covered planet going around a nuclear fireball 90
% million miles away and think this to be normal is obviously some
% indication of how skewed our perspective tends to be.
%
% Douglas Adams

