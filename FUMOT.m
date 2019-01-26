classdef FUMOT < handle
    %FUMOT 
    
    properties (Access = public)
        RTE
        sigma
        gamma
        boundary
        H
        model
        
        beta
        
        edge % in case boundary values are excluded
        
        curXF
        h
    end
    
    methods (Access = public)
        function obj = FUMOT(opt)
            obj.RTE = rte(opt);
            
            if ~(isfield(opt, 'sigma') )
                error('Reference to non-existing field:FUMOT');
            else
                % for the first stage only.
                % function handles in opt. OUTPUT rows vector.
                obj.sigma.xs = opt.sigma.xs(obj.RTE.nodes);
                obj.sigma.xa = opt.sigma.xa(obj.RTE.nodes);
                obj.sigma.xf = opt.sigma.xf(obj.RTE.nodes);
                obj.sigma.ms = opt.sigma.ms(obj.RTE.nodes);
                obj.sigma.ma = opt.sigma.ma(obj.RTE.nodes);
                obj.sigma.eta = opt.sigma.eta(obj.RTE.nodes);
            end
            
            if ~(isfield(opt, 'gamma'))
                obj.gamma.x = 1;
                obj.gamma.f = 1;
            else
                % constants, default values are 1. This situation is not
                % supported anymore.
                obj.gamma.x  = opt.gamma.x;
                obj.gamma.f  = opt.gamma.f;  
            end
           
            
            if ~(isfield(opt, 'boundary') )
                error('Reference to non-existing field:FUMOT');
            else
                obj.boundary = opt.boundary; 
                obj.h        = opt.h; % for emission.
            end            
            
            obj.beta = 1e-3;
            tau = opt.tau;
            
            [obj.H, ~] = obj.ExciteForwardOp(obj.sigma.xf);
            
            % polution
            obj.H =  obj.H .* (1 + tau * (2 * rand(size(obj.H)) - 1));
            
            % edges idx
            obj.edge = unique(obj.RTE.segms);
        end
        
        function S = ScatterOp(obj, X)
            rX = reshape(X, obj.RTE.nAngle, obj.RTE.nPoint);
            S  = ifft(bsxfun(@times, obj.RTE.pFHG', fft(rX))) * obj.RTE.dtheta;  
        end
        
        function resetBoundaryCondition(obj)
            obj.RTE.boundarySource = zeros(obj.RTE.nAngle, obj.RTE.nPoint);
        end
        
        function [H, Up, psi] = ExciteForwardOp(obj, XF)    
            obj.RTE.setBoundaryCondition(obj.boundary);
            sigmaXTF = obj.sigma.xa + XF + obj.sigma.xs; ...
            sigmaS   = obj.sigma.xs;

            obj.RTE.sigmaT = sigmaXTF;
            obj.RTE.sigmaS = sigmaS;

            nAngle = obj.RTE.nAngle;

            % The rte routine gives the adjoint solution.
            % U_{-} = u(x, -v);
            % Reverse the directions by flipping the bottom half to top.
            % U_{+} = u(x, v) as the exact solution.

            Un = obj.RTE.ForwardSolve(); 
            Up = [Un(nAngle/2+1:end, :) ; Un(1:nAngle/2, :)]; 
            psi = sum(Un .* Up               , 1)/ nAngle;           
            phi = sum(Un .* obj.ScatterOp(Up), 1)/nAngle;

            H = -obj.gamma.f * sigmaXTF .* psi + obj.gamma.x * sigmaS .* phi;
            
            
        end
       
        function [f, g] = ExciteBackwardOp(obj, XF)
            
            % LOAD INTO CLASS, then it is OK to stop at anytime.
            obj.curXF = XF;
            
            [H0, Up, psi] = ExciteForwardOp(obj, XF');
            mismatch = H0 - obj.H; % row vec
            f = 0.5 * mismatch * obj.RTE.M * mismatch' + 0.5 * obj.beta * XF' * obj.RTE.S * XF;
            g = obj.beta * obj.RTE.S * XF;
            
            % first part in gradient
            g = g - obj.RTE.M * (mismatch .* psi)';
            
            % solves adjoint equation with internal source
            sigmaXTF = obj.sigma.xa + XF' + obj.sigma.xs; ...
            sigmaS   = obj.sigma.xs;
        
            % the transpose is to make sure output as angle x space.
            Q  = - bsxfun(@times, 2 * sigmaXTF', Up')' + ...
                bsxfun(@times, 2 * sigmaS', (obj.ScatterOp(Up))' )';
            
            Z = -bsxfun(@times, mismatch', Q')';
            Y = reshape(  obj.RTE.rays.interior_transport(...
                obj.RTE.nodes, obj.RTE.elems, obj.RTE.sigmaT,  Z), obj.RTE.nAngle * obj.RTE.nPoint, 1);
            % solve the transport equation
            forwardMap = @(X) (X - obj.RTE.mapping(reshape(X, obj.RTE.nAngle, obj.RTE.nPoint))  );
            [Xp, ~] = gmres(forwardMap, Y, 10, 1e-12, 400); 
            
            Xp = reshape(Xp, obj.RTE.nAngle, obj.RTE.nPoint);
            Xn = [Xp(obj.RTE.nAngle/2+1:end, :) ; Xp(1:obj.RTE.nAngle/2, :)]; 
            
            T = sum(Up .* Xn               , 1)/ obj.RTE.nAngle; 
            g = g + obj.RTE.M * T';
            
            
            %g(obj.edge) = 0;
  
        end
        
        function [xf,hist]  = ExciteSolve(obj, XF)
            
%             options = optimset('Diagnostics','on','DerivativeCheck','off',...
%                 'FinDiffType','central','LargeScale','off',...
%                 'GradObj','on','Display','iter-detailed',...
%                 'TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',800,....
%                 'MaxIter',300,'HessUpdate','bfgs');
%           
%             [xf,fval,exitflag,output] = fminunc(@obj.ExciteBackwardOp, XF,options);


            options    = struct( 'factr', 1e0, 'pgtol', 1e-7, 'm', 20, ...
                'x0',XF, 'maxIts', 300, 'maxTotalIts', 800);
            
            options.printEvery     = 1;            
            
            [xf, ~, hist] =...
                lbfgsb_c(@obj.ExciteBackwardOp, -inf * ones(size(XF)), inf * ones(size(XF)), options);  
%             
        end
        
        
        function [S, wp, phi] = EmissionForwardOp(obj, ETA, H)
            % solve interior source problem instead.
             obj.RTE.setBoundaryCondition(obj.h);
             sigmaMT = obj.sigma.ma + obj.sigma.ms; ...
             sigmaS   = obj.sigma.ms;

             obj.RTE.sigmaT = sigmaMT;
             obj.RTE.sigmaS = sigmaS;

             nAngle = obj.RTE.nAngle;
             
             Wn = obj.RTE.ForwardSolve(); 
             Wp = [Un(nAngle/2+1:end, :) ; Un(1:nAngle/2, :)]; % aux function.
             
             
             
             
            
        end

        
    end
    
end

