% Coherence selection function.
%
% <http://spindynamics.org/wiki/index.php?title=Coherence.m>

function rho=coherence(spin_system,rho,spec)

% Check consistency
grumble(spin_system,rho,spec);

% Store dimension statistics
spn_dim=size(spin_system.bas.basis,1);
spc_dim=numel(rho)/spn_dim;
problem_dims=size(rho);

% Fold indirect dimensions
rho=reshape(rho,[spn_dim spc_dim]);

% Compute projection quantum numbers of basis states
[~,M]=lin2lm(spin_system.bas.basis);

% Preallocate state mask array
state_mask=false(spn_dim,numel(spec));

% Loop over specifications
for n=1:numel(spec)

    % Parse spin specification
    if ischar(spec{n}{1})
        
        % Symbolic specification
        if strcmp(spec{n}{1},'all')
            spins=1:numel(spin_system.comp.isotopes);
        else
            spins=find(strcmp(spec{n}{1},spin_system.comp.isotopes));
        end
        
    else
        
        % Specification by number
        spins=spec{n}{1};
        
    end
    
    % Determine coherence order of each basis state
    coherence_orders_present=sum(M(:,spins),2);
  
    % Wipe all coherence orders except those specified by the user
    state_mask(:,n)=ismember(coherence_orders_present,spec{n}{2});
        
end

% Intersect state masks
state_mask=all(state_mask,2);

% Apply the state mask
rho(~state_mask,:)=0;

% Unfold indirect dimensions
rho=reshape(rho,problem_dims);

% Report overly destructive calls
if norm(rho,1)<1e-10
    report(spin_system,'WARNING - all magnetization appears to have been destroyed by this call.');
end

end

% Consistency enforcement
function grumble(spin_system,rho,spec) %#ok<INUSD>
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('analytical coherence order selection is only available for sphten-liouv formalism.');
end
if ~isnumeric(rho)
    error('the state vector(s) must be numeric.');
end
if mod(numel(rho),size(spin_system.bas.basis))~=0
    error('the number of elements in rho must be a multiple of the dimension of the spin state space.');
end
end

% There was a time when men were afraid that somebody would reveal
% some secret of theirs that was unknown to their fellows. Nowadays,
% they're afraid that somebody will name what everybody knows. Have
% you practical people ever thought that that's all it would take
% to blast your whole, big, complex structure, with all your laws
% and guns -- just somebody naming the exact nature of what it is
% you're doing?
%
% Ayn Rand, "Atlas Shrugged"

