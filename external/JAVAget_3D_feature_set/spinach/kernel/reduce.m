% Symmetry and trajectory-level state space reduction for a user-
% specified Liouvillian L and initial state rho. Applies all reduc-
% tion methods (unless disabled during the call to create.m) and
% returns a cell array of projectors into a set of independent re-
% duced subspaces. Syntax:
%
%              projectors=reduce(spin_system,L,rho)
%
% where L is the Liouvillian and rho is the initial state (in the
% case of source state screening) or the detection state (if des-
% tination state screening is used). The output is a cell array of
% projectors into independently evolving reduced subspaces. Those
% projectors are to be used as follows:
%
%             L_reduced=P'*L*P;    (for matrices)
%             rho_reduced=P'*rho;  (for state vectors)
%
% Further information is available here:
%
%         http://dx.doi.org/10.1016/j.jmr.2008.08.008
%         http://link.aip.org/link/doi/10.1063/1.3398146
%         http://dx.doi.org/10.1016/j.jmr.2011.03.010
%
% i.kuprov@soton.ac.uk

function projectors=reduce(spin_system,L,rho)

% Check the input
grumble(spin_system,L,rho);

% Check for blanket ban
if ismember('trajlevel',spin_system.sys.disable)
    report(spin_system,'WARNING - trajectory level state space restriction disabled by user.');
    projectors{1}=speye(size(L)); return
end

% Decide how to proceed
switch spin_system.bas.formalism
    
    case 'zeeman-hilb'
        
        % If a cell array is supplied, make a representative density matrix
        if iscell(rho)
            rho_rep=abs(rho{1});
            for n=2:numel(rho)
                rho_rep=rho_rep+abs(rho{n});
            end
            rho=rho_rep/numel(rho);
        end
        
        % Run symmetry factorization
        if ismember('symmetry',spin_system.sys.disable)
            
            % Issue a reminder to the user
            report(spin_system,'WARNING - permutation symmetry treatment disabled by user.');
            
            % Return unit matrix
            projectors{1}=speye(size(L));
            
        elseif ~isfield(spin_system.bas,'irrep')
            
            % Inform the user
            report(spin_system,'no permutation symmetry information has been supplied.');
            
            % Return unit matrix
            projectors{1}=speye(size(L));
            
        else
            
            % Decide which irreps to keep
            n_irreps=numel(spin_system.bas.irrep); irrep_keep_index=true(n_irreps,1);
            for n=1:n_irreps
                
                % Check the irrep contribution to the total norm
                if spin_system.bas.irrep(n).dimension==0
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' is has dimension zero - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                elseif norm(spin_system.bas.irrep(n).projector'*rho*spin_system.bas.irrep(n).projector,1)<spin_system.tols.irrep_drop
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(spin_system.bas.irrep(n).dimension) ' has less than '...
                                        num2str(spin_system.tols.irrep_drop) ' of the state norm - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                else
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(spin_system.bas.irrep(n).dimension) ', is active - kept.']);
                    
                end
                
            end
            
            % Compile the projector array
            projectors={spin_system.bas.irrep(irrep_keep_index).projector};
            
            % Flatten out the cell array
            projectors=[projectors{:}];
        
        end
        
    case {'zeeman-liouv','sphten-liouv'}

        % If a stack is supplied, choose a representative state vector
        if size(rho,2)>1, rho=mean(abs(rho),2); end
        
        % Run symmetry factorization and ZTE
        if ismember('symmetry',spin_system.sys.disable)
            
            % Issue a reminder to the user
            report(spin_system,'WARNING - permutation symmetry treatment disabled by user.');
            
            % Run zero track elimination
            report(spin_system,'attempting zero track elimination...');
            projectors{1}=zte(spin_system,L,rho);
            
        elseif ~isfield(spin_system.bas,'irrep')
            
            % Inform the user
            report(spin_system,'no permutation symmetry information has been supplied.');
            
            % Run zero track elimination
            report(spin_system,'attempting zero track elimination...');
            projectors{1}=zte(spin_system,L,rho);
            
        else
            
            % Decide which irreps to keep
            n_irreps=numel(spin_system.bas.irrep); irrep_keep_index=true(n_irreps,1);
            for n=1:n_irreps
                
                % Check the irrep contribution to the total norm
                if spin_system.bas.irrep(n).dimension==0
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ' is has dimension zero - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                elseif norm(spin_system.bas.irrep(n).projector'*rho,1)<spin_system.tols.irrep_drop
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(spin_system.bas.irrep(n).dimension) ' has less than '...
                                        num2str(spin_system.tols.irrep_drop) ' of the state norm - dropped.']);
                    
                    % Flag the irrep for dropping
                    irrep_keep_index(n)=0;
                    
                else
                    
                    % Update the user
                    report(spin_system,['irrep #' num2str(n) ', dimension '...
                                        num2str(spin_system.bas.irrep(n).dimension) ' is active - kept.']);
                    
                end
                
            end
            
            % Compile the projector array
            projectors={spin_system.bas.irrep(irrep_keep_index).projector};
            
            % Loop over permutation group irreps
            for n=1:numel(projectors)
                
                % Report to the user
                report(spin_system,['irrep #' num2str(n) ', attempting zero track elimination...']);
                
                % Run zero track elimination
                zte_projector=zte(spin_system,projectors{n}'*L*projectors{n},projectors{n}'*rho);
                
                % Project the projectors
                projectors{n}=projectors{n}*zte_projector;
                
            end
            
        end
        
        % Run path tracing
        for n=1:numel(projectors)
            
            % Inform the user
            report(spin_system,['path-tracing subspace #' num2str(n) '...']);
            
            % Run the path tracing
            pt_projectors=path_trace(spin_system,projectors{n}'*L*projectors{n},projectors{n}'*rho);
            
            % Project the projectors
            for k=1:numel(pt_projectors)
                pt_projectors{k}=projectors{n}*pt_projectors{k};
            end
            projectors{n}=pt_projectors; %#ok<AGROW>
            
        end
        
        % Flatten out the cell array
        projectors=[projectors{:}];
        
    otherwise
        
        % Complain and bomb out
        error('unknown formalism specification.');
        
end

end

% Consistency enforcement
function grumble(spin_system,L,rho) %#ok<INUSD,INUSL>
if ~isnumeric(L)
    error('L must be numeric.');
end
if size(L,1)~=size(L,2)
    error('L must be a square matrix.');
end
end

% For centuries, the battle of morality was fought between those who claimed
% that your life belongs to God and those who claimed that it belongs to your
% neighbours � between those who preached that the good is self-sacrifice for
% the sake of ghosts in heaven and those who preached that the good is self-
% sacrifice for the sake of incompetents on earth. And no one came to say that
% your life belongs to you and the good is to live it.
%
% Ayn Rand, "Atlas Shrugged"

