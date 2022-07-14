classdef app1_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        GridLayout              matlab.ui.container.GridLayout
        LeftPanel               matlab.ui.container.Panel
        UIAxes2                 matlab.ui.control.UIAxes
        UIAxes                  matlab.ui.control.UIAxes
        RightPanel              matlab.ui.container.Panel
        DropDown                matlab.ui.control.DropDown
        DropDownLabel           matlab.ui.control.Label
        BeatsDropDown           matlab.ui.control.DropDown
        BeatsDropDownLabel      matlab.ui.control.Label
        FilenamesTextArea       matlab.ui.control.TextArea
        FilenamesTextAreaLabel  matlab.ui.control.Label
        ECGDropDown             matlab.ui.control.DropDown
        ECGDropDownLabel        matlab.ui.control.Label
        FILTERButton            matlab.ui.control.Button
        RESETButton             matlab.ui.control.Button
        LoadBeatButton          matlab.ui.control.Button
        SPLITButton             matlab.ui.control.Button
        LOADSIGNALButton        matlab.ui.control.Button
        UIAxes_2                matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        signal % signal property
        smoothed % smoothed signal
    end
    
    methods (Access = private)
        
        
        function baseline_removed = baseline_wandering(app)
            x=app.signal;
            n=2;
            span=floor(length(x)/n);
            ini=0;
            window_length=span;
            initial=window_length;
                
            for i=1:1:(floor(length(x))/window_length)
                prev_value_all=0;
                
                m=1;
                for k=(1+ini):1:(window_length+ini)
                    
                    data_all(m,i)=x(k);
                    m=m+1;
                end
                
                %%% data detetrending portion %%%%%
                
                order_pol=10;
                t=(1:length(data_all(:,i)));
                
                [p,s,mu]=polyfit(t', data_all(:,i), order_pol);
                f_y=polyval(p,t,[],mu);
                dt_ecgnl=data_all(:,i)-f_y';
                %%%%%%%%%%%%% fresh Data %%%%%%%%%%%  
                
                fresh_data(:,i)=dt_ecgnl;
                 
                
            ini=ini+initial;
            end
            for x=1:1:n   
              re_data(x,:)=[zeros(1,(span*(x-1))) fresh_data(:,x)' zeros(1,(span*(n-x)))];
            end
            final_data=zeros(1,n*span);
            for xx=1:1:n
                final_data=final_data+re_data(xx,:);
            end
            baseline_removed = [final_data];
            end    
        end
    
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LOADSIGNALButton
        function LOADSIGNALButtonPushed(app, event)
            data = app.ECGDropDown.Value;
            mat=load(data);
            app.signal=mat.val(1,:);
            N= length(app.signal); 
            fs=360 ;
            T=1/fs; 
            seconds=linspace(0,(N*T),N); 
            
            plot(app.UIAxes,seconds,app.signal);
        end

        % Button pushed function: FILTERButton
        function FILTERButtonPushed(app, event)
            %app.smoothed = smoothdata(app.signal,"movmedian","SmoothingFactor",0.5);
            
            filtered_signal=baseline_wandering(app);
           
            app.smoothed = smoothdata(filtered_signal,"movmedian","SmoothingFactor",0.5);
            N= length(app.smoothed); 
            fs=360 ;
            T=1/fs; 
            seconds=linspace(0,(N*T),N); 
            plot(app.UIAxes2,seconds,app.smoothed);
            
            
        end

        % Button pushed function: RESETButton
        function RESETButtonPushed(app, event)
            cla(app.UIAxes)
            cla(app.UIAxes2)
            cla(app.UIAxes_2);
            app.BeatsDropDown.Items=" "
            app.FilenamesTextArea.Value=" "

        end

        % Button pushed function: SPLITButton
        function SPLITButtonPushed(app, event)
            y=app.smoothed;
            Fs=360;
            peak_threshold=0.3;
            Q_th = round((20/360)*Fs);
            Q_th_1st = round((60/360)*Fs);
            Q_th_2nd = round((5/360)*Fs);
            P_th_1st = round((40/360)*Fs);
            P_th_2nd = round((1/360)*Fs);
            S_th_1st = round((5/360)*Fs);
            S_th_2nd = round((60/360)*Fs);
            T_th_1st = round((40/360)*Fs);
            T_th_2nd = round((120/360)*Fs);
            % ------------------QRS detection by FS1 and FS2 Algorithm ---------%
            y0       = diff(y);
            for i = 2:length(y)-2
                y1(i) = (y0(i-1) + 2*y0(i) + y0(i+1))/4;
            end
            y1(1) = 0;
            y2 = diff(y0);
            for i = 1:length(y)-2
                y3(i) = y1(i) + y2(i);
            end
             
            th = peak_threshold*max(y3);
            QRS = [];
            for i = 2:length(y3)-6
                if (y3(i) > th)  &&  y3(i) > y3(i+1) && y3(i) > y3(i-1)    
                    QRS = [QRS i];
                end
            end
            QRS_Actual = [];
            for i = 2:length(QRS)-1
                if QRS(i+1) -QRS(i) > (Fs*0.6)
                   QRS_Actual = [ QRS_Actual QRS(i) ];
                end
            end
            % -------------------------------------------------------------%
            %  ---------- R peak detection in original signal ------------------%
            R_L = [];
            for i = 1:length(QRS_Actual)-1
                a = QRS_Actual(i)-Q_th: QRS_Actual(i)+S_th_2nd;
                m = max(y(a));
                b = find(y(a)== m);
                if numel(b) ~= 0
                    b = b(1);
                    b = a(b);
                    R_L = [R_L b];
                end
            end
            % -------------------------------------------------------------------%
            % ----------- Q peak detection signal -------------------------------%
            Q_L = [];
            for i = 1:length(R_L)
                a = R_L(i)-P_th_1st: R_L(i)-P_th_2nd;
                m = min(y(a));
                b = find(y(a)== m);
                if numel(b) ~= 0
                    b = b(1);
                    b = a(b);
                    Q_L = [Q_L b];
                end
            end
            % ------------------------------------------------------------------ %
            % -------------- S peak detection -----------------------------------%
            S_L = [];
            for i = 1:length(R_L)
                x = R_L(i)+S_th_1st:R_L(i)+S_th_2nd;
                m = min(y(x));
                p = find(y(x) == m);
                if numel(p) ~= 0
                    p = p(1);
                    p = x(p);
                    S_L = [S_L p];
                end
            end
            % -----------------------------------------------------------------%
            % ------------- P Peak detection ----------------------------------------%
            P_L = [];
            for i = 1:length(Q_L)
                a = Q_L(i)-Q_th_1st: Q_L(i)-Q_th_2nd;
                m = max(y(a));
                b = find(y(a)==m);
                b = b(1);
                b = a(b);
                P_L = [P_L b];
            end
            % ----------------------------------------------------------------------- %
            % ----------------- T peak detection ------------------------------------%
            T_L = [];
            for i = 1:length(S_L)
                a = S_L(i)+T_th_1st: S_L(i)+T_th_2nd;
                m = max(y(a));
                b = find(y(a)==m);
                b = b(1);
                b = a(b);
                T_L = [T_L b];
            end
            %-----------------------------------------------------------------------%
            % -----------------------------------------------------------------%
            abc_out = [P_L; T_L];
            
            for i =1:length(P_L)
                beat = app.smoothed(P_L(i):T_L(i));
                filename=strcat(num2str(i),'_beat.mat')
                %fn=str(filename);
                 N= length(beat); 
                 fs=360 ;
                 T=1/fs; 
                 seconds=linspace(0,(N*T),N); 
                plot(app.UIAxes_2,seconds,beat);
                hold(app.UIAxes_2,"on");
                app.BeatsDropDown.Items{i}=filename;
                

                app.FilenamesTextArea.Value{i+1} = filename ;
                save(filename,'beat');
            end
            hold(app.UIAxes_2,"off");
           %plot(app.UIAxes_2,beat)
            
        end

        % Button pushed function: LoadBeatButton
        function LoadBeatButtonPushed(app, event)
            data = app.BeatsDropDown.Value;
            mat=load(data);
            
            plot(app.UIAxes_2,mat.beat);
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {634, 634};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {413, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.4745 0.8588 0.8275];
            app.UIFigure.Position = [100 100 1092 634];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {413, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create UIAxes
            app.UIAxes = uiaxes(app.LeftPanel);
            title(app.UIAxes, 'ECG SIGNAL')
            xlabel(app.UIAxes, 'Seconds')
            ylabel(app.UIAxes, 'Amplitude(mV)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [1 316 395 305];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.LeftPanel);
            title(app.UIAxes2, 'FILTERED ECG SIGNAL')
            xlabel(app.UIAxes2, 'Seconds')
            ylabel(app.UIAxes2, 'Amplitude(mV)')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [1 8 389 295];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.RightPanel);
            title(app.UIAxes_2, 'HEART BEAT')
            xlabel(app.UIAxes_2, 'Seconds')
            ylabel(app.UIAxes_2, 'Amplitude(mV)')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.Position = [45 316 447 305];

            % Create LOADSIGNALButton
            app.LOADSIGNALButton = uibutton(app.RightPanel, 'push');
            app.LOADSIGNALButton.ButtonPushedFcn = createCallbackFcn(app, @LOADSIGNALButtonPushed, true);
            app.LOADSIGNALButton.Position = [384 223 100 22];
            app.LOADSIGNALButton.Text = 'LOAD SIGNAL';

            % Create SPLITButton
            app.SPLITButton = uibutton(app.RightPanel, 'push');
            app.SPLITButton.ButtonPushedFcn = createCallbackFcn(app, @SPLITButtonPushed, true);
            app.SPLITButton.Position = [384 136 100 22];
            app.SPLITButton.Text = 'SPLIT';

            % Create LoadBeatButton
            app.LoadBeatButton = uibutton(app.RightPanel, 'push');
            app.LoadBeatButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBeatButtonPushed, true);
            app.LoadBeatButton.Position = [384 50 100 22];
            app.LoadBeatButton.Text = 'LOAD BEAT';

            % Create RESETButton
            app.RESETButton = uibutton(app.RightPanel, 'push');
            app.RESETButton.ButtonPushedFcn = createCallbackFcn(app, @RESETButtonPushed, true);
            app.RESETButton.Position = [384 7 100 22];
            app.RESETButton.Text = 'RESET';

            % Create FILTERButton
            app.FILTERButton = uibutton(app.RightPanel, 'push');
            app.FILTERButton.ButtonPushedFcn = createCallbackFcn(app, @FILTERButtonPushed, true);
            app.FILTERButton.Position = [384 179 100 22];
            app.FILTERButton.Text = 'FILTER';

            % Create ECGDropDownLabel
            app.ECGDropDownLabel = uilabel(app.RightPanel);
            app.ECGDropDownLabel.HorizontalAlignment = 'right';
            app.ECGDropDownLabel.Position = [337 267 32 22];
            app.ECGDropDownLabel.Text = 'ECG';

            % Create ECGDropDown
            app.ECGDropDown = uidropdown(app.RightPanel);
            app.ECGDropDown.Items = {'300m001.mat', '300m002.mat', '300m003.mat', '300m004.mat', '300m005.mat', '300m010.mat', '300m011.mat', '300m012.mat', '300m013.mat', '300m014.mat', '300m015.mat', '300m020.mat'};
            app.ECGDropDown.Placeholder = 'Pick Signal';
            app.ECGDropDown.Position = [384 267 100 22];
            app.ECGDropDown.Value = '300m001.mat';

            % Create FilenamesTextAreaLabel
            app.FilenamesTextAreaLabel = uilabel(app.RightPanel);
            app.FilenamesTextAreaLabel.HorizontalAlignment = 'right';
            app.FilenamesTextAreaLabel.Position = [44 262 61 22];
            app.FilenamesTextAreaLabel.Text = 'Filenames';

            % Create FilenamesTextArea
            app.FilenamesTextArea = uitextarea(app.RightPanel);
            app.FilenamesTextArea.Position = [120 49 208 237];

            % Create BeatsDropDownLabel
            app.BeatsDropDownLabel = uilabel(app.RightPanel);
            app.BeatsDropDownLabel.HorizontalAlignment = 'right';
            app.BeatsDropDownLabel.Position = [333 93 36 22];
            app.BeatsDropDownLabel.Text = 'Beats';

            % Create BeatsDropDown
            app.BeatsDropDown = uidropdown(app.RightPanel);
            app.BeatsDropDown.Position = [384 93 100 22];

            % Create DropDownLabel
            app.DropDownLabel = uilabel(app.RightPanel);
            app.DropDownLabel.HorizontalAlignment = 'right';
            app.DropDownLabel.Position = [491 430 66 22];
            app.DropDownLabel.Text = 'Drop Down';

            % Create DropDown
            app.DropDown = uidropdown(app.RightPanel);
            app.DropDown.Position = [571 430 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end