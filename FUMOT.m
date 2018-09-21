classdef FUMOT < handle
    %FUMOT 
    
    properties (Access = public)
        RTE
        sigma
    end
    
    methods (Access = public)
        function obj = FUMOT(opt)
            obj.RTE = rte(opt);
            
            if ~(isfield(opt, 'FuncXS') && ...
                    isfield(opt, 'FuncXA') && isfield(opt, 'FuncXF'))
                error('Reference to non-existing field:FUMOT');
            end
            
            % for the first stage only.
            obj.sigma.xs = opt.FuncXS(obj.RTE.nodes)';
            obj.sigma.xa = opt.FuncXA(obj.RTE.nodes)';
            obj.sigma.xf = opt.FuncXF(obj.RTE.nodes)';
            
            
            
        end

        
        
        
    end
    
end

