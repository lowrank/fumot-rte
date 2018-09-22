classdef FUMOT < handle
    %FUMOT 
    
    properties (Access = public)
        RTE
        sigma
        gamma
        boundary
    end
    
    methods (Access = public)
        function obj = FUMOT(opt)
            obj.RTE = rte(opt);
            
            if ~(isfield(opt, 'sigma') )
                error('Reference to non-existing field:FUMOT');
            end
            
            if ~(isfield(opt, 'gamma') )
                error('Reference to non-existing field:FUMOT');
            end
            
            if ~(isfield(opt, 'boundary') )
                error('Reference to non-existing field:FUMOT');
            end            
           
            
            % for the first stage only.
            % function handles in opt. OUTPUT rows vector.
            obj.sigma.xs = opt.sigma.xs(obj.RTE.nodes);
            obj.sigma.xa = opt.sigma.xa(obj.RTE.nodes);
            obj.sigma.xf = opt.sigma.xf(obj.RTE.nodes);

            obj.boundary = opt.boundary; 
            
            % constants
            obj.gamma.x  = opt.gamma.x;
            obj.gamma.f  = opt.gamma.f;
            
        end
        
        function [H, psi] = ExciteForwardOp(obj, XF)
            obj.RTE.setBoundaryCondition(obj.boundary);
            obj.RTE.sigmaT = obj.sigma.xa + XF + obj.sigma.xs; ...
            obj.RTE.sigmaS = obj.sigma.xs;
            
            nAngle = obj.RTE.nAngle;
            
            % The rte routine gives the adjoint solution.
            % U_{-} = u(x, -v);
            Un = obj.RTE.ForwardSolve(); 
            
            % Reverse the directions by flipping the bottom half to top.
            % U_{+} = u(x, v) as the exact solution.
            Up = [Un(nAngle/2+1:end, :) ; Un(1:nAngle/2, :)]; 
            psi = sum(Un .* Up , 1) / nAngle;
            
            H = obj.gamma.f * XF .* psi;
            
        end
       
        
        % do a stochatic gradient?
        
    end
    
end

