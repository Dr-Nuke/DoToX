function SliderDemo2(stack,col,dim,F,cas)
% adjust lines 23+ for loading the data, and each call of imshow()
% stack is a 3d array
% col is a [ , ] color range for imshow()
% dim is the dimension that will be slided



%// Check below for dummy 4D matrix/image sequence
dims= ndims(stack); % number of dimensions
if dim > dims
    error('not that many dimensions')
end

for i=1:dims
    frames{i}=1:size(stack,i);
    if i==dim
        frames{i}=1;
        NumFrames = size(stack,i); 
    end
end


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
%recon=load('diff');  % load a variable from disk
MyMatrix=stack;   % rename to MyMatrix 
%col=[0,1];           % color range for imshow()

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.

setappdata(hFig,'MyMatrix',MyMatrix); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 

%// Display 1st frame
imshow(squeeze(MyMatrix(frames{1},frames{2},frames{3}))',col);set(gca,'YDir','normal')

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
        frames{dim}=CurrentFrame;
        imshow(squeeze(MyMatrix(frames{1},frames{2},frames{3}))',col,'Parent',handles.axes1);
        hold on
        scatter(squeeze(F.cfit(cas,frames{3},:,1)),...
            squeeze(F.cfit(cas,frames{3},:,2)),'+','MarkerEdgeColor',[1 1 1])
        hold off
        set(gca,'YDir','normal')
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

        imshow(squeeze(MyMatrix(frames{1},frames{2},frames{3}))',col,'Parent',handles.axes1);
        hold on
        scatter(squeeze(F.cfit(cas,frames{3},:,1)),...
            squeeze(F.cfit(cas,frames{3},:,2)),'+','MarkerEdgeColor',[1 1 1])
        hold off
        set(gca,'YDir','normal')
        guidata(hFig,handles);
    end

end