classdef VS_mrGratings < VStim
    properties (SetAccess=public)
        
        grating_width       = 53;       %pixel/cycle, choose sp=225, tp=5 for 30deg/s at 60x
        tmporalFreq         = 2;
        txtStmpFrq          = 2;        %cycles/sec
        txtSmaskRadius      = 500;      %radius of circular mask
        txtSnumDirs         = 8;
        txtSnumTrials       = 5;
        txtSduration        = 4;        %txtSduration in secs
        txtSpreStimWait     = 5;
        txtSdelay           = 2;
        txtSinterTrialWait  = 2;
        save_stimulus       = true;
        chkSmakeGrating     = true;
        txtSdirList         = nan;
        txtSbrtIntensity    = 255;
        popSgrtColor        = [1 1 1];
        txtSdrkIntensity    = 0;
        chkScorrectGrating  = true;
        x_shift             = 0;
        y_shift             = 0;
        chkSmaskRect        = false;
        txtSrectWidth       = 100;
        txtSrectHeight      = 400;
        txtSscrIntensity    = 139;
        popSscrColor        = [1 1 1];
        chkSsaveImage       = false;
        txtSsaveImageTime   = 2;
    end
    properties (Hidden,Constant)
        defaultTrialsPerCategory=50; %number of gratings to present
        defaultBackground=0;
        grating_widthTxt        = "scalar,width of black+white grating, pixels/cycle";
        txtStmpFrqTxt           = "scalar,cycles/sec";
        txtSmaskRadiusTxt       = "scalar,half size of square texture";
        txtSnumDirsTxt          = "scalar,number of directions";
        txtSnumTrialsTxt        = "scalar,number of repetitions";
        txtSdurationTxt         = "scalar,txtSduration of each repetition in secs";
        txtSpreStimWaitTxt      = "scalar,seconds of black screen before stimulus start";
        txtSdelayTxt            = "scalar,seconds of black screen before each trial start";
        txtSinterTrialWaitTxt   ="scalar,seconds of waiting between trials";
        save_stimulusTxt        = "0 or 1,save stimuli?";
        chkSmakeGratingTxt      = "0 or 1,1 creates grating instead of sinusoid";
        txtSdirListTxt          = "vector, list of stimulus directions in degrees if nan (default), generates evenly spaced directions based on txtSnumDirs argument.";
        txtSbrtIntensityTxt     = "scalar, 0 to 255. Maximum intensity of grating or sinusoid";
        txtSdrkIntensityTxt     = "scalar, 0 to 255. minimum intensity of grating or sinusoid";
        chkScorrectGratingTxt   = "0 or 1, if 0, makes half rectified sinusoids instead of squarewaves like old method.";
        shift_xTxt              = "scalar, shifts the center of stimuli in the x axis";
        shift_yTxt              = "scalar, shifts the center of stimuli in the y axis";
        chkSmaskRectTxt         = "0 or 1. if 1, the masking of the grating is rectangular";
        txtSrectWidthTxt        = "in case of a rectangular mask";
        txtSrectHeightTxt       = "in case of a rectangular mask";
        txtSscrIntensityTxt     = "color of baground. default is gray";
        remarks                 = {'Categories in stimuli are: speed, offset'};
    end
    properties (SetAccess=protected)
        speeds
        %         directions
        offsets
        barTrajectories1X
        barTrajectories2X
        barTrajectories1Y
        barTrajectories2Y
        nFrames
    end
    properties (Hidden, SetAccess=protected)
        flip
        stim
        flipEnd
        miss
    end
    methods
        
        function obj=run(obj)
            screen_full_color = obj.txtSscrIntensity*obj.popSscrColor;
            
            % configure directions
            
            if isnan(obj.txtSdirList) || isempty(obj.txtSdirList)
                dirs = 0:(360/obj.txtSnumDirs):(360-(360/obj.txtSnumDirs));
            else
                dirs = obj.txtSdirList;
                dims = size(dirs);
                if dims(1)>dims(2)
                    dirs = dirs';
                end
            end
            
            directions = [];
            for r=1:obj.txtSnumTrials
                dirs_shuffled = Shuffle(dirs);
                directions=cat(2,directions,dirs_shuffled);
            end
            directions=directions';
            
            % configure correct mask
            
            if obj.chkSmaskRect
                obj.txtSmaskRadius = 0;
            end
            
            Screen('BlendFunction', obj.PTB_win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('FillRect', obj.PTB_win, screen_full_color, []);
            Screen('Flip',obj.PTB_win);
            
            [width, height]=Screen('WindowSize', obj.PTB_win);
            
            if ~obj.chkSmaskRect
                mask = makeCircularMaskForGUI(obj.txtSmaskRadius,width, height,'color',screen_full_color);
            else
                obj.txtSmaskRadius = max(ceil(obj.txtSrectWidth/2),ceil(obj.txtSrectHeight/2));
                mask = makeRectangularMaskForGUI(obj.txtSrectWidth,obj.txtSrectHeight);
            end
            
            mask = mask(:,:,[2,4]);
            masktex=Screen('MakeTexture', obj.PTB_win, mask);
            
            % Calculate parameters of the Sinusoid
            
            gray=round((obj.txtSbrtIntensity+obj.txtSdrkIntensity)/2);
            inc=obj.txtSbrtIntensity-gray;
            
            p=ceil(obj.grating_width);  % pixels/cycle
            fr=(1/obj.grating_width)*2*pi;
            visiblesize=2*obj.txtSmaskRadius+1;
            
            [x,~]=meshgrid(-obj.txtSmaskRadius:obj.txtSmaskRadius + abs(p), -obj.txtSmaskRadius:obj.txtSmaskRadius);
            sinusoid = gray + inc*cos(fr*x); %change x to fix mightex polygon
            
            if obj.chkSmakeGrating
                grating_mask = sinusoid>gray; %returns a 0 and 1 mask
                if obj.chkScorrectGrating
                    sinusoid(grating_mask) = obj.txtSbrtIntensity; %makes actual squarewave
                    sinusoid(~grating_mask) = obj.txtSdrkIntensity;
                else
                    sinusoid = sinusoid.*grating_mask; %old method generates half rectified sinusoids
                end
            end
            
            if any(obj.popSgrtColor ~= [1 1 1])
                color_sinusoid = repmat(sinusoid,1,1,3);
                color_sinusoid = color_sinusoid.*reshape(obj.popSgrtColor,1,1,3);
            else
                color_sinusoid = sinusoid;
            end
            
            sinusoidtex=Screen('MakeTexture', obj.PTB_win, color_sinusoid);
            
            % configure display window
            
            dstRect=[0 0 visiblesize visiblesize];
            dstRect=CenterRect(dstRect, obj.rect);
            [x_center, y_center] = RectCenterd(dstRect);
            dstRect = CenterRectOnPointd(dstRect,x_center+obj.x_shift,y_center+obj.y_shift);
            shiftperframe= obj.txtStmpFrq * p * obj.ifi;
            
            obj.sendTTL(1,true);
            WaitSecs(obj.txtSpreStimWait);
            
            % stimulation start
            for trial=1:length(directions)
                direction=directions(trial);
                direction=mod((direction+180),360);               
                
                i=0;
                WaitSecs(obj.txtSdelay);
                vbl=Screen('Flip', obj.PTB_win);
                vblendtime = vbl + obj.txtSduration;
                saveImageTime=obj.txtSsaveImageTime+vbl;
                %%%%%%
                while (vbl < vblendtime)
                    xoffset = mod(i*shiftperframe,abs(p));
                    i=i+1;
                    srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
                    Screen('DrawTexture', obj.PTB_win, sinusoidtex, srcRect, dstRect, direction);
                    Screen('DrawTexture',  obj.PTB_win, masktex, [0 0 visiblesize visiblesize], dstRect, direction);
                    obj.applyBackgound;
                    obj.sendTTL(3,true);
                    vbl=Screen('Flip',  obj.PTB_win);
                    obj.sendTTL(3,false);
                    
                    if obj.chkSsaveImage
                        if vbl>=saveImageTime || round(vbl)>=saveImageTime
                            obj.chkSsaveImage=0;
                            imageArray=Screen('GetImage', window, screenRect);
                            imwrite(imageArray,[datestr(now,'yymmdd'),'_sinusoidsScr_vbl',num2str(vbl),...
                                '_time',num2str(obj.txtSsaveImageTime),'.bmp'],'bmp');
                        end
                    end
                    
                end
                Screen('Flip', obj.PTB_win);
                WaitSecs(obj.txtSinterTrialWait);%pause to allow recording device to prepare for new trial
                obj.sendTTL(2,false);
                disp(['Direction ' num2str(trial) '/' num2str(obj.txtSnumTrials*obj.txtSnumDirs)]);
            end
            obj.applyBackgound;
            Screen('Flip', obj.PTB_win);
            obj.sendTTL(1,false);
            if obj.save_stimulus
                filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrGratings_%s.mat', datestr(now,'mm_dd_yyyy_HHMM'));
                save(filename, 'directions', 'obj', '-v7.3');
            end
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats.props=obj.getProperties;
        end
        %class constractor
        function obj=VS_mrGratings(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
        end
    end
end %EOF