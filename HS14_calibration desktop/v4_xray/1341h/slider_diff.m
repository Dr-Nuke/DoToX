function SliderDemo
% adjust lines 23+ for loading the data, and each call of imshow()

clc
clear all

NumFrames = 966; %// Check below for dummy 4D matrix/image sequence
hFig = figure('Position',[100 100 500 500],'Units','normalized');

handles.axes1 = axes('Units','normalized','Position',[.2 .2 .6 .6]);

%// Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol('Style','slider','Position',[60 20 400 50],'Min',1,'Max',NumFrames,'Value',1,'SliderStep',[1/NumFrames 2/NumFrames],'Callback',@XSliderCallback);
handles.SliderxListener = addlistener(handles.SliderFrame,'Value','PostSet',@(s,e) XListenerCallBack);

handles.Text1 = uicontrol('Style','Text','Position',[180 420 60 30],'String','Current frame');
handles.Edit1 = uicontrol('Style','Edit','Position',[250 420 100 30],'String','1');

%// Create dummy image sequence, here 4D sequence of grayscale images.
% MyImage = imread('peppers.png');
% 
% MyMatrix = cat(4,rgb2gray(MyImage),MyImage(:,:,1),MyImage(:,:,2),MyImage(:,:,3));
recon=load('diff');  % load a variable from disk
MyMatrix=recon.diff;   % rename to MyMatrix 
col=[0,0.01];           % color range for imshow()

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.

setappdata(hFig,'MyMatrix',MyMatrix); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 

%// Display 1st frame
imshow(MyMatrix(:,:,1),col)

%// IMPORTANT. Update handles structure.
guidata(hFig,handles);

%// Listener callback, executed when you drag the slider.

    function XListenerCallBack

        %// Retrieve handles structure. Used to let MATLAB recognize the
        %// edit box, slider and all UI components.
        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
MyMatrix = getappdata(hFig,'MyMatrix');

        %// Get current frame
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        %// Display appropriate frame.
        imshow(MyMatrix(:,:,CurrentFrame),col,'Parent',handles.axes1);

        guidata(hFig,handles);
    end


%// Slider callback; executed when the slider is release or you press
%// the arrows.
    function XSliderCallback(~,~)

        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
    MyMatrix = getappdata(hFig,'MyMatrix');

        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        imshow(MyMatrix(:,:,CurrentFrame),col,'Parent',handles.axes1);

        guidata(hFig,handles);
    end

end