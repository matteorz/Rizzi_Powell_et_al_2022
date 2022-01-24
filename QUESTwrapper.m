classdef QUESTwrapper < handle
    % Palamedes wrapper
    
    properties (GetAccess = 'public', SetAccess = 'private')
        pObj
    end
    
    %% ====================================================================
    %  -----PUBLIC METHODS-----
    %$ ====================================================================
    
    methods (Access = 'public')
        
        %% == CONSTRUCTOR =================================================
        
        function obj=QUESTwrapper(qparams)
            qparams.priorAlphaRange = qparams.deltaMin:qparams.deltaStepsize:qparams.deltaMax;
            qparams.PF = eval(qparams.PF);
            obj.pObj = PAL_AMRF_setupRF(PAL_AMRF_setupRF(), 'priorAlphaRange',qparams.priorAlphaRange, 'beta',qparams.beta, 'lambda',qparams.lambda, 'gamma',qparams.gamma, 'PF',qparams.PF, 'stopCriterion',qparams.stopCriterion, 'stopRule',qparams.stopRule);
            % test (don't save result)
            PAL_AMRF_updateRF(obj.pObj, obj.pObj.xCurrent, false);
        end
        % Destructor
        function obj = delete(obj)
            clear obj;
        end
        
        %% == METHODS =====================================================
        
        function [delta] = getDelta(obj)
            delta = obj.pObj.xCurrent;
        end
        
        function fin = isFinished(obj)
            fin = obj.pObj.stop();
        end

        function [] = update(obj, wasCorrect)
             obj.pObj = PAL_AMRF_updateRF(obj.pObj, obj.pObj.xCurrent, wasCorrect);
        end
        
        function [est] = estThresh(obj)
            est = obj.pObj.xCurrent;
        end
        
    end
    
    %% ====================================================================
    %  -----STATIC METHODS-----
    %$ ====================================================================
    
    methods(Static)
        % useful when debugging
        function params = getDummyParams()
            params = [];
            params.startVal = 100;
            params.stepSize = [10 5 2.5]; % 1./[1.1 1.05 1.001];
            params.downMod = [1 .5 .5];
            params.nReversals = [4 8 8];
            params.nUp = 1; % [1 1 1];
            params.nDown = 1; % [1 1 1];
            params.isAbsolute = true;
            params.minVal = 50;
            params.maxVal = 200;
            params.minNTrials = 0;
            params.maxNTrials = 50;
            params.verbosity = 2;
        end
    end
    
    %% ====================================================================
    %  -----PRIVATE METHODS-----
    %$ ====================================================================
    
    methods(Access = 'private')
        
        function figHandles = createPlot(obj, maxN,initialVal,minVal,maxVal)
            figHandles.hFig=figure(length(findobj('Type','figure'))+1); % create on top of any existing
            set(figHandles.hFig, 'Position', [300 100 600 800]); % [x y width height]
            hold on
            figHandles.hPerf = plot(-1,-1);
            figHandles.hRight = plot(-1,-1 ...
                ,'o' ...
                ,'LineWidth',2 ...
                ,'MarkerFaceColor','g' ...
                ,'MarkerEdgeColor','k' ...
                ,'MarkerSize',10);
            figHandles.hWrong = plot(-1,-1 ...
                ,'o' ...
                ,'LineWidth',2 ...
                ,'MarkerFaceColor','r' ...
                ,'MarkerEdgeColor','k' ...
                ,'MarkerSize',10);
            figHandles.hNextTarget = plot(1,initialVal ...
                ,'x' ...
                ,'MarkerEdgeColor','k' ...
                ,'MarkerSize',10);
            figHandles.vlines = cell(1, length(obj.stepSize));
            hold off
            if maxN>1; xlim([1 maxN]); end
            ylim([minVal, maxVal]);
            xlabel('Trial Number','FontSize',16);
            ylabel('\Delta','FontSize',16);
            set(figHandles.hFig,'Name','Adaptive Track','NumberTitle','off'); % set window title
            set(gcf, 'Renderer', 'Painters'); % Painters? OpenGL? OpenGL conflicts with psychtoolbox  
        end
        
        function [] = updatePlot(obj)
            % update markers
            trialNums = 1:length(obj.deltaHistory);
            set(obj.figHandles.hPerf,'XData',trialNums,'YData',obj.deltaHistory);
            set(obj.figHandles.hRight,'XData',trialNums(obj.wasCorrectHistory==1),'YData',obj.deltaHistory(obj.wasCorrectHistory==1));
            set(obj.figHandles.hWrong,'XData',trialNums(obj.wasCorrectHistory==0),'YData',obj.deltaHistory(obj.wasCorrectHistory==0));
            
            % highlight reversals (each stage in different colour)
            nStages = length(obj.stepSize);
            
           	axes(gca);   
            for i = 1:nStages
                if ishandle(obj.figHandles.vlines{i})
                    delete(obj.figHandles.vlines{i});
                end
                idx = (obj.stageHistory==i) & (obj.reversalsHistory~=0);
                if sum(idx)>0
                    obj.figHandles.vlines{i} = vline(find(idx), obj.stageColours{i});
                end
            end

            % update next target value [cross]
            set(obj.figHandles.hNextTarget,'XData',length(obj.deltaHistory)+1,'YData',obj.getDelta());
               
            % refresh graphics
            drawnow();
        end
        
    end
end