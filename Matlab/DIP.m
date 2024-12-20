function varargout = DIP(varargin)
% DIP MATLAB code for DIP.fig
%      DIP, by itself, creates a new DIP or raises the existing
%      singleton*.
%
%      H = DIP returns the handle to a new DIP or the handle to
%      the existing singleton*.
%
%      DIP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DIP.M with the given input arguments.
%
%      DIP('Property','Value',...) creates a new DIP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DIP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DIP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DIP

% Last Modified by GUIDE v2.5 19-Dec-2024 15:30:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DIP_OpeningFcn, ...
                   'gui_OutputFcn',  @DIP_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before DIP is made visible.
function DIP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DIP (see VARARGIN)

% Choose default command line output for DIP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DIP wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = DIP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end





% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
    % 创建一个文件对话框来选择图像文件
    [filename, pathname] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.gif','Image Files (*.jpg, *.jpeg, *.png, *.bmp, *.gif)'});
    % 检查用户是否选择了文件
    if isequal(filename,0) || isequal(pathname,0)
        disp('User selected Cancel');
        return;
    end
    % 构造图像的完整路径
    fullFileName = fullfile(pathname, filename);
    % 读取图像
    img = imread(fullFileName);
    % 在 GUI 的 Axes 中显示图像
    axes(handles.axes1); % 确保在正确的 Axes 中显示图像
    imshow(img);
    % 可选：设置图像标题
     title('原图像');
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
    % 获取当前在 GUI 中显示的图像
    % 假设图像显示在名为 'axes1' 的坐标轴上
    var=get(handles.popupmenu5,'value');
%获取弹出式菜单的value值
switch var
    case 1
      axesHandle = handles.axes1; % 假设你的 Axes 句柄存储在 handles 结构中
    case 2
            axesHandle = handles.axes5; % 假设你的 Axes 句柄存储在 handles 结构中
    case 3
            axesHandle = handles.axes2; % 假设你的 Axes 句柄存储在 handles 结构中
    case 4  
            axesHandle = handles.axes9; % 假设你的 Axes 句柄存储在 handles 结构中
    case 5
            axesHandle = handles.axes10; % 假设你的 Axes 句柄存储在 handles 结构中
    case 6
            axesHandle = handles.axes6; % 假设你的 Axes 句柄存储在 handles 结构中
    case 7
            axesHandle = handles.axes3; % 假设你的 Axes 句柄存储在 handles 结构中
    case 8
            axesHandle = handles.axes4; % 假设你的 Axes 句柄存储在 handles 结构中
    case 9  
            axesHandle = handles.axes7; % 假设你的 Axes 句柄存储在 handles 结构中
    case 10
            axesHandle = handles.axes8; % 假设你的 Axes 句柄存储在 handles 结构中     
end

    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        img = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    

   % 弹出保存对话框
    [filename, pathname, filterIndex] = uiputfile({'*.png;*.PNG','PNG Files (*.png)';
                                                   '*.jpg;*.JPEG','JPEG Files (*.jpg)';
                                                   '*.bmp;*.BMP','BMP Files (*.bmp)'}, 'Save As', ...
                                                fullfile(pwd, 'untitled'));
    
    % 检查用户是否取消了保存操作
    if isequal(filename,0)
        return;
    end
    
    % 根据用户选择的文件类型设置保存格式
    switch filterIndex
        case 1
            format = 'png';
        case 2
            format = 'jpg';
        case 3
            format = 'bmp';
        otherwise
            format = 'png'; % 默认格式
    end
    
    % 构造完整的文件路径
    fullFilePath = fullfile(pathname, [filename, '.', format]);
    
    % 保存图像
    imwrite(img, fullFilePath);
    
    % 可选：显示保存成功的消息框
    uiwait(msgbox('Image saved successfully!', 'Success', 'info'));
end
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
      
    axesHandle = handles.axes1; % 假设你的 Axes 句柄存储在 handles 结构中
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        img = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    

    % 将图像转换为灰度图像
    a = rgb2gray(img);
    
    
[x,y]=size(a);
h=zeros(1,256);
for i=1:256    
    h(i)=length(find(a==(i-1)))/(x*y);
    %用find依次找到a数组中0到255灰度值的向量，在用length计算出0到255灰度值向量所对应的个数，
    %在除以像素总个数，最后将得到的数值p依次存入矩阵h
end
axes(handles.axes2)
     
bar(0:255,h);
title('灰度化后图像的直方图');
%绘制原图像的直方图
s1=zeros(1,256);
s2=zeros(1,256);
t=0;
for i=1:256
    t=t+h(i);
    s1(i)=t;
    %对数值p进行累加计算出新的灰度级
end
s2=round(s1*255);
%将s1中的灰度级修正为合理的灰度级
s3=zeros(1,256);
for i=1:256
    s3(i)=sum(h(s2==(i-1)));
    %在h矩阵中找到与修正后的灰度级对应的灰度级所对应的比重值，如原始灰度级对应同样修正后的灰度级则将他们相加。
end
axes(handles.axes4)
   
bar(0:255,s3);
  title('均衡化后图像的直方图');
 %绘制均衡化后图像的直方图
c=a;
for i=1:256
    c(find(a==(i-1)))=s2(i);  
 %将均衡化后的灰度值值一一对应给图像c，从而获得直方图均衡化后的图像
end
axes(handles.axes3)
 
imshow(c) ,   title('直方图均衡化后的图像');
end
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)

    axesHandle = handles.axes1; % 假设你的 Axes 句柄存储在 handles 结构中
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        a = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    
  
[x,y,z]=size(a);
s=ones(x,y);
R=im2double(a(:,:,1));
G=im2double(a(:,:,2));
B=im2double(a(:,:,3));
%提取图像的RGB三色分量
for i=1:x
    for j=1:y
        s(i,j)=max(a(i,j,:));
    end
end
%通过遍历找到最大的像素值
        axes(handles.axes5) 
        i=round((R+B+G)/3*255);
        a(:,:,1)=i;
        a(:,:,2)=i;
        a(:,:,3)=i;
        %使分个分量值等于平均值
        imshow(a),  title('灰度化后的图像');
end

% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

    axesHandle = handles.axes5; % 假设你的 Axes 句柄存储在 handles 结构中
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        a = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    
a=a(:,:,1);
a1=double(a);
var=get(handles.popupmenu1,'value');
%获取弹出式菜单的value值
switch var
    case 1
        axes(handles.axes6)
        a2=10*log(a1+1);
        %以10为底的对数变换
        imshow(a2,[]),title("对数变换")
    case 2
        axes(handles.axes6)
         a3=100*exp(0.325*(a1-220)/30)-1;
          %以e为底的指数函数
        imshow(a3,[]),title("指数变换")
    case 3
        axes(handles.axes6)
          a3=a1.^3;
          %以像素值为底，3为指数的幂次函数
        imshow(a3,[]),title("幂数变换")
    case 4
        axes(handles.axes6)
          % 定义分段线性变换的参数
% 这里我们定义三个区间：[0, 64], [65, 128], [129, 255]
% 以及对应的线性变换参数（斜率和截距）
intervals = [0, 64; 65, 128; 129, 255];
slopes = [1.5, 0.5, 2]; % 每个区间的斜率
intercepts = [-32, 32, -64]; % 每个区间的截距
 
% 初始化输出图像
J = zeros(size(a), 'uint8');
 
% 应用分段线性变换
for y = 1:size(a, 1)
    for x = 1:size(a, 2)
        gray_value = double(a(y, x));
        % 找到gray_value所在的区间
        for i = 1:size(intervals, 1)-1
            if gray_value >= intervals(i, 1) && gray_value <= intervals(i+1, 1)
                % 计算变换后的灰度值
                new_value = slopes(i) * (gray_value - intervals(i, 1)) + intercepts(i) + intervals(i, 1);
                % 截断到有效范围 [0, 255]
                if new_value < 0
                    new_value = 0;
                elseif new_value > 255
                    new_value = 255;
                end
                J(y, x) = uint8(new_value);
                break;
            end
        end
    end
end           
        imshow(J),title("线性变换")
end
end
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end



function [image_out, T3] = histmatching(image_in, hist)
% 直方图匹配（规定化）函数
% 输入为需要进行直方图匹配的灰度图像，模板直方图
% 输出为直方图匹配后的灰度图像，进行变换的向量。

% Level为灰度级别
% T1, T2分别为输入图像，模板直方图的均衡化用到的变换向量
% T3为输入图像匹配模板直方图用到的变换向量
Level = 256;
[m, n] = size(image_in);
image_hist = imhist(image_in);
image_out = image_in;
% 求解T1
ac1 = zeros(Level, 1);
T1 = zeros(Level, 1, 'uint8');
ac1(1) = image_hist(1);
for i = 2 : Level
    ac1(i) = ac1(i - 1) + image_hist(i);
end
ac1 = ac1 * (Level - 1);
for i = 1 : 256
    T1(i) = uint8(round((ac1(i)) / (m * n)));
end

% 求解T2
ac2 = zeros(Level, 1);
T2 = zeros(Level, 1, 'uint8');
ac2(1) = hist(1);
for i = 2 : Level
    ac2(i) = ac2(i - 1) + hist(i);
end
ac2 = ac2 * (Level - 1);
hist_sum = sum(hist);
for i = 1 : 256
    T2(i) = uint8(round((ac2(i)) / hist_sum));
end

% 求解T3
% T1映射到T2^(-1)时，若有多个值，选取最小的那个值。
% 产生0 到 255 之间的256个点，即产生0,1,2,...,255的大小为256的数组
temp = zeros(Level, 1, 'uint8');
T3 = T1;
for i = 1 : 256
    for j = 1 : 256
        temp(j) = abs(T1(i) - T2(j));
    end
    [~, B] = min(temp);
    T3(i) = B - 1;
end

% 根据T3转换输入图像的值
for i = 1 : m
    for j = 1 : n
        image_out(i, j) = T3(uint32(image_in(i, j)) + 1);
    end
end
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % 创建一个文件对话框来选择图像文件
    [filename, pathname] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.gif','Image Files (*.jpg, *.jpeg, *.png, *.bmp, *.gif)'});
    
    % 检查用户是否选择了文件
    if isequal(filename,0) || isequal(pathname,0)
        disp('User selected Cancel');
        return;
    end
    
    % 构造图像的完整路径
    fullFileName = fullfile(pathname, filename);
    
    % 读取图像
    img = imread(fullFileName);
    
    % 在 GUI 的 Axes 中显示图像
    axes(handles.axes7); % 确保在正确的 Axes 中显示图像
    imshow(img);
    
    % 可选：设置图像标题
     title('模板图像');
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)

    axesHandle = handles.axes5; % 假设你的 Axes 句柄存储在 handles 结构中
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        a = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    
image2=a(:,:,1);
    axesHandle = handles.axes7; % 假设你的 Axes 句柄存储在 handles 结构中
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        image3 = get(imgObj, 'CData'); % 获取 CData
    else
        error('No image object found on the axes.');
    end    


hist3 = imhist(image3);
%match1 = histeq(image2, hist3);
[match2, T] = histmatching(image2, hist3);
axes(handles.axes7); % 确保在正确的 Axes 中显示图像
    imshow(image3);
    
    % 可选：设置图像标题
     title('模板图像');
axes(handles.axes3); imshow(match2), title('匹配后得到的图像');
axes(handles.axes8); imhist(image3), title('模板图像的直方图');
axes(handles.axes4); imhist(match2), title('匹配后的直方图');
end
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function pushbutton8_Callback(hObject, eventdata, handles)
    axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        f = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

    [h, w, c] = size(f); % 获取图像的高度、宽度和颜色通道数
    kx = str2double(get(handles.edit2, 'String')); % 获取水平缩放因子
    ky = str2double(get(handles.edit4, 'String')); % 获取垂直缩放因子
    H = ceil(ky * h); % 计算新图像的高度
    W = ceil(kx * w); % 计算新图像的宽度

    % 初始化新图像矩阵
    g1 = zeros(H, W, c, 'uint8'); % 使用uint8类型，因为图像数据通常是8位

    % 插值计算新图像的像素值
    for newx = 1:W
        for newy = 1:H
            oldx = (newx - 0.5) / kx + 0.5; % 中心点映射，避免边缘效应
            oldy = (newy - 0.5) / ky + 0.5;
            x = floor(oldx);
            y = floor(oldy);
            a = oldx - x; % 水平插值系数
            b = oldy - y; % 垂直插值系数

            % 检查边界条件
            if x >= 1 && x < w && y >= 1 && y < h
                % 双线性插值
                for ch = 1:c
                    g1(newy, newx, ch) = (1-a)*(1-b)*f(y, x, ch) + a*(1-b)*f(y, x+1, ch) + ...
                                         (1-a)*b*f(y+1, x, ch) + a*b*f(y+1, x+1, ch);
                end
            else
                % 边界处理，简单复制最近的像素值
                x_clip = max(1, min(x, w-1));
                y_clip = max(1, min(y, h-1));
                g1(newy, newx, :) = f(y_clip, x_clip, :);
            end
        end
    end

    axes(axesHandle); % 确保在正确的 Axes 中显示图像
    imshow(g1); % 显示处理后的图像
    title('处理后图像');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
    a=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(a) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','1');
    end
    guidata(hObject,handles); % 更新数据

end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
    axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        f = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

    ang = str2double(get(handles.edit3, 'String')) ;
 
      g2=imrotate(f,ang,'bilinear','crop');
    axes(axesHandle); % 确保在正确的 Axes 中显示图像
    imshow(g2), title('处理后图像'); % 显示处理后的图像
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
    a=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(a) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','0');
    end
    guidata(hObject,handles); % 更新数据
end

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
    a=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(a) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','1');
    end
    guidata(hObject,handles); % 更新数据
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        a = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

J=imnoise(a,'gauss',0,0.02);
% 添加高斯噪声
J1= imnoise(a,'salt & pepper', 0.02) ;
% 添加椒盐噪声
J3=imnoise(a, 'poisson');
var=get(handles.popupmenu2,'value');
%获取弹出式菜单的value值
switch var
    case 1
        axes(handles.axes9)
        imshow(J),title("高斯噪声")
    case 2
        axes(handles.axes9)
        imshow(J1),title("椒盐噪声")
    case 3
        axes(handles.axes9)
        imshow(J3),title("泊松噪声")
end  
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axesHandle = handles.axes9; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        i = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

I=double(i);
I=fftshift(fft2(I));
%傅里叶变换及频谱中心化
[a,b]=size(I);
g=zeros(a,b);
a0=round(a/2);
b0=round(b/2);
D0=str2double(get(handles.edit5, 'String'));


var=get(handles.popupmenu3,'value');

switch var
    case 1
        axes(handles.axes10)
        h = ones(3,3) / 9;  % 3x3的均值滤波器   
        a1 = imfilter(i, h);  % 应用滤波器  
        %a=fspecial('average',3);
        %生成3*3模板的均值滤波
        %a1=uint8(filter2(a,i));
        %进行中值滤波运算 
        imshow(a1),title('均值滤波')
    case 2 
        i = rgb2gray(i);  %转换为灰度图像（如果需要）
        b1= medfilt2(i);  % 使用2D中值滤波器
        %对图像进行3*3模板的中值滤波计算
        axes(handles.axes10)
        imshow(b1),title('中值滤波')
    case 3
        %c=fspecial('gaussian',3);
        %生成3*3模板的高斯滤波
        c1= imgaussfilt(i, 2);  % 使用标准差为2的高斯滤波器  % 应用高斯滤波器  
        %进行高斯滤波运算 
        axes(handles.axes10)
        imshow(c1),title('高斯滤波')
    case 4
 I = im2double(i);
 M = 2*size(I,1);%滤波的行数
 N = 2*size(I,2);%滤波的列数
 u = -M/2:(M/2-1);
 v = -N/2:(N/2-1);
 [u,v] = meshgrid(u,v);
 D = sqrt(u.^2+v.^2);

 H = double(D<D0);%理想低通滤波器
 J = fftshift(fft2(I,size(H,1),size(H,2)));%时域图像转换到频域
 K = J.*H;%滤波处理
 L = ifft2(ifftshift(K));%傅里叶反变换
 L = L(1:size(I,1),1:size(I,2));

    axes(handles.axes10)

    imshow(L);title("理想低通滤波器");
    case 5
I = im2double(i);
 M = 2*size(I,1);
 N = 2*size(I,2);
 u = -M/2:(M/2-1);
 v = -N/2:(N/2-1);
 [U,V] = meshgrid(u,v);
 D = sqrt(U.^2+V.^2);

 n = 6;
 H = 1-exp(-(D.^2)./(2*(D0^2)));%设计滤波器
 J = fftshift(fft2(I,size(H,1),size(H,2)));
 K = J.*H;
 L = ifft2(ifftshift(K));
 L = L(1:size(I,1),1:size(I,2));

    axes(handles.axes10)

    imshow(L);title("高斯高通滤波器");
    case 6
 I = im2double(i);
 M = 2*size(I,1);
 N = 2*size(I,2);
 u = -M/2:(M/2-1);
 v = -N/2:(N/2-1);
 [u,v] = meshgrid(u,v);
 D = sqrt(u.^2+v.^2);

 n = 6;
 H = 1./(1+(D./D0).^(2*n));%设计巴特沃斯滤波器
 J = fftshift(fft2(I,size(H,1),size(H,2)));
 K = J.*H;
 L = ifft2(ifftshift(K));
 L = L(1:size(I,1),1:size(I,2));

    axes(handles.axes10)
    imshow(L); title("巴特沃斯低通滤波器")

       
 end

end



% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
    b=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(b) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','30');
    end
    guidata(hObject,handles); % 更新数据
end

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        I = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

var=get(handles.popupmenu4,'value');
switch var
    case 1
        h=[1 0;0 -1];h1=[0 1;-1 0];
        %Roberts算子模板 
        t=imfilter(I,h);t1=imfilter(I,h1);
        %对图像进行相应滤波运算 
        J=abs(t)+abs(t1);
        %将结果取绝对值后相加 
        axes(handles.axes10)
        imshow(J),title("Robert算子 ")
        
    case 2
        h2=[-1 -2 -1;0 0 0;1 2 1];h3=[-1 0 1;-2 0 2;-1 0 1];
        %Sobel算子模板 
        t2=imfilter(I,h2);t3=imfilter(I,h3);
        %对图像进行相应滤波运算 
        J1=abs(t2)+abs(t3);
        %将结果取绝对值后相加 
        axes(handles.axes10)
        imshow(J1),title("Sober算子")
     case 3
        H1=[-1 -2 -1;0 0 0;1 2 1];
        H2=[0 -1 -2;1 0 -1; 2 1 0];
        H3=[1 0 -1;2 0 -2;1 0 -1];
        H4=[2 1 0;1 0 -1;0 -1 -2];
        H5=[1 2 1;0 0 0;-1 -2 -1];
        H6=[0 1 2;-1 0 1;-2 -1 0];
        H7=[-1 0 1;-2 0 2;-1 0 1];
        H8=[-2 -1 0;-1 0 1;0 1 2];
        %Sobel算子扩展
         %对图像进行Sobel每个扩展模板的运算 
        R1=imfilter(I,H1);
        R2=imfilter(I,H2);
        R3=imfilter(I,H3);
        R4=imfilter(I,H4);
        R5=imfilter(I,H5);
        R6=imfilter(I,H6);
        R7=imfilter(I,H7);
        R8=imfilter(I,H8);
        f1=max(max(R1,R2),max(R3,R4));
        f2=max(max(R5,R6),max(R7,R8));
        %找到每行的最大值 
        a=max(f1,f2);
        axes(handles.axes10)
        imshow(a),title("8模板Sobel算子")
     case 4
        h=[-1 -1 -1;0 0 0;1 1 1];h1=[-1 0 1;-1 0 1;-1 0 1];
        %Prewitt算子模板 
        t=imfilter(I,h);t1=imfilter(I,h1);
         %对图像进行相应滤波运算 
        J2=abs(t)+abs(t1);
         axes(handles.axes10)
        imshow(J2),title("Prewitt算子 ")
       case 5
         h=[0 -1 0;-1 4 -1;0 -1 0];
         %Laplacian算子模板
         b=imfilter(I,h);
         %对图像进行相应滤波运算 
         axes(handles.axes10)
         imshow(b),title("Laplacian算子")    
end
end




% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
end

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        i = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end
i=im2double(i);
L=superpixels(i,500);
ho=drawfreehand;
foreground=createMask(ho,i);
hb=drawfreehand;
background=createMask(hb,i);
BW=lazysnapping(i,L,foreground,background);
         axes(handles.axes5)
         imshow(BW),title("分割结果");
mask=cat(3,BW,BW,BW);
result=mask.*i;
         axes(handles.axes2)
         imshow(result),title("目标提取")    
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        i = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end
image= rgb2gray(i);
[N,M] = size(image);
lbp=zeros(N,M);
for j=2:N-1
for i=2:M-1
neighbor=[j-1 i-1;j-1 i;j-1 i+1;j i+1;j+1 i+1;j+1 i;j+1 i-1;j i-1];
count=0;
for k=1:8
if image(neighbor(k,1),neighbor(k,2))> image(j,i)
count = count + 2^(8-k);
end
end

lbp(j,i)=count;
end
end
lbp = uint8(lbp);
         axes(handles.axes3)
         imshow(lbp),title( '原图像LBP 特征图');
axesHandle = handles.axes2; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        i = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end
image= rgb2gray(i);
[N,M] = size(image);
lbp=zeros(N,M);
for j=2:N-1
for i=2:M-1
neighbor=[j-1 i-1;j-1 i;j-1 i+1;j i+1;j+1 i+1;j+1 i;j+1 i-1;j i-1];
count=0;
for k=1:8
if image(neighbor(k,1),neighbor(k,2))> image(j,i)
count = count + 2^(8-k);
end
end

lbp(j,i)=count;
end
end
lbp = uint8(lbp);
         axes(handles.axes4)
         imshow(lbp),title( '目标LBP 特征图');
end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        img = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end
% 1、%灰度化
img=rgb2gray(img);
img=double(img);

step=8;      %step*step个像素作为一个cell
[m1 n1]=size(img);
%改变图像尺寸为step的最近整数倍
img=imresize(img,[floor(m1/step)*step,floor(n1/step)*step],'nearest');
[m n]=size(img);

% 2、%伽马校正
img=sqrt(img);

% 3、求梯度和方向
fy=[-1 0 1];        %定义竖直模板
fx=fy';             %定义水平模板
Iy=imfilter(img,fy,'replicate');    %竖直梯度
Ix=imfilter(img,fx,'replicate');    %水平梯度
Ied=sqrt(Ix.^2+Iy.^2);              %梯度值
Iphase=Iy./Ix;              %边缘斜率，有些为inf,-inf,nan，其中nan需要再处理一下
the=atan(Iphase)*180/3.14159; %求梯度角度

for i=1:m
    for j=1:n
        if(Ix(i,j)>=0&Iy(i,j)>=0) %第一象限
            the(i,j)=the(i,j);
        elseif(Ix(i,j)<=0&Iy(i,j)>=0) %第二象限
            the(i,j)=the(i,j)+180;
        elseif(Ix(i,j)<=0&Iy(i,j)<=0) %第三象限
            the(i,j)=the(i,j)+180;
        elseif(Ix(i,j)>=0&Iy(i,j)<=0) %第四象限
            the(i,j)=the(i,j)+360;
        end

        if isnan(the(i,j))==1  %0/0会得到nan，如果像素是nan，重设为0
            the(i,j)=0;
        end

    end
end
the=the+0.000001; %防止角度为0

% 4、划分cell，求cell的直方图( 1 cell = 8*8 pixel )
clear i j;
%下面是求cell
step=8;                %step*step个像素作为一个cell
orient=9;               %方向直方图的方向个数
jiao=360/orient;        %每个方向包含的角度数
Cell=cell(1,1);              %所有的角度直方图,cell是可以动态增加的，所以先设了一个
ii=1;
jj=1;

for i=1:step:m
    ii=1;
    for j=1:step:n
        Hist1(1:orient)=0;
        for p=1:step
            for q=1:step
                %梯度方向直方图
                Hist1(ceil(the(i+p-1,j+q-1)/jiao))=Hist1(ceil(the(i+p-1,j+q-1)/jiao))+Ied(i+p-1,j+q-1);
            end
        end
        Cell{ii,jj}=Hist1;       %放入Cell中
        ii=ii+1;
    end
    jj=jj+1;
end

% 5、划分block，求block的特征值,使用重叠方式( 1 block = 2*2 cell )
clear m n i j;
[m n]=size(Cell);
feature=cell(1,(m-1)*(n-1));
for i=1:m-1
    for j=1:n-1
        block=[];
        block=[Cell{i,j}(:)' Cell{i,j+1}(:)' Cell{i+1,j}(:)' Cell{i+1,j+1}(:)'];
        block=block./sum(block); %归一化
        feature{(i-1)*(n-1)+j}=block;
    end
end

% 6、图像的HOG特征值
[m n]=size(feature);
l=2*2*orient;
featureVec=zeros(1,n*l);
for i=1:n
    featureVec((i-1)*l+1:i*l)=feature{i}(:);
end

%到此结束，feature即为所求
%下面是为了显示而写的
l=length(feature);
f=[];
for i=1:l
    f=[f;feature{i}(:)'];  
end 
         axes(handles.axes7)
        mesh(f) ,title( '原图像HOG特征图');
   set(handles.edit6, 'String', num2str(length(featureVec)));
axesHandle = handles.axes2; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        img = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end
% 1、%灰度化
img=rgb2gray(img);
img=double(img);

step=8;      %step*step个像素作为一个cell
[m1 n1]=size(img);
%改变图像尺寸为step的最近整数倍
img=imresize(img,[floor(m1/step)*step,floor(n1/step)*step],'nearest');
[m n]=size(img);

% 2、%伽马校正
img=sqrt(img);

% 3、求梯度和方向
fy=[-1 0 1];        %定义竖直模板
fx=fy';             %定义水平模板
Iy=imfilter(img,fy,'replicate');    %竖直梯度
Ix=imfilter(img,fx,'replicate');    %水平梯度
Ied=sqrt(Ix.^2+Iy.^2);              %梯度值
Iphase=Iy./Ix;              %边缘斜率，有些为inf,-inf,nan，其中nan需要再处理一下
the=atan(Iphase)*180/3.14159; %求梯度角度

for i=1:m
    for j=1:n
        if(Ix(i,j)>=0&Iy(i,j)>=0) %第一象限
            the(i,j)=the(i,j);
        elseif(Ix(i,j)<=0&Iy(i,j)>=0) %第二象限
            the(i,j)=the(i,j)+180;
        elseif(Ix(i,j)<=0&Iy(i,j)<=0) %第三象限
            the(i,j)=the(i,j)+180;
        elseif(Ix(i,j)>=0&Iy(i,j)<=0) %第四象限
            the(i,j)=the(i,j)+360;
        end

        if isnan(the(i,j))==1  %0/0会得到nan，如果像素是nan，重设为0
            the(i,j)=0;
        end

    end
end
the=the+0.000001; %防止角度为0

% 4、划分cell，求cell的直方图( 1 cell = 8*8 pixel )
clear i j;
%下面是求cell
step=8;                %step*step个像素作为一个cell
orient=9;               %方向直方图的方向个数
jiao=360/orient;        %每个方向包含的角度数
Cell=cell(1,1);              %所有的角度直方图,cell是可以动态增加的，所以先设了一个
ii=1;
jj=1;

for i=1:step:m
    ii=1;
    for j=1:step:n
        Hist1(1:orient)=0;
        for p=1:step
            for q=1:step
                %梯度方向直方图
                Hist1(ceil(the(i+p-1,j+q-1)/jiao))=Hist1(ceil(the(i+p-1,j+q-1)/jiao))+Ied(i+p-1,j+q-1);
            end
        end
        Cell{ii,jj}=Hist1;       %放入Cell中
        ii=ii+1;
    end
    jj=jj+1;
end

% 5、划分block，求block的特征值,使用重叠方式( 1 block = 2*2 cell )
clear m n i j;
[m n]=size(Cell);
feature=cell(1,(m-1)*(n-1));
for i=1:m-1
    for j=1:n-1
        block=[];
        block=[Cell{i,j}(:)' Cell{i,j+1}(:)' Cell{i+1,j}(:)' Cell{i+1,j+1}(:)'];
        block=block./sum(block); %归一化
        feature{(i-1)*(n-1)+j}=block;
    end
end

% 6、图像的HOG特征值
[m n]=size(feature);
l=2*2*orient;
featureVec=zeros(1,n*l);
for i=1:n
    featureVec((i-1)*l+1:i*l)=feature{i}(:);
end

%到此结束，feature即为所求
%下面是为了显示而写的
l=length(feature);
f=[];
for i=1:l
    f=[f;feature{i}(:)'];  
end 
         axes(handles.axes8)
        mesh(f) ,title( '目标HOG特征图');
           set(handles.edit7, 'String', num2str(length(featureVec)));
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
    b=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(b) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','0');
    end
    guidata(hObject,handles); % 更新数据
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
    b=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(b) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','0');
    end
    guidata(hObject,handles); % 更新数据
end
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
end

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        axesHandle = handles.axes1; % 获取 Axes 句柄
    imgObj = findobj(axesHandle, 'Type', 'image'); % 找到 Axes 上的 Image 对象
    if ~isempty(imgObj)
        img = get(imgObj, 'CData'); % 获取 CData，假设这是RGB图像
    else
        error('No image object found on the axes.');
    end

% 设置工作目录
%cd('C:\path\to\your\project'); % 替换为您的项目路径

% 导入 SavedModel 模型
modelPath = 'D://Model/final_model_savedmodel_finetuned'; % 替换为您的模型路径
net = importTensorFlowNetwork(modelPath);

% 加载类别映射文件
classesFilePath = 'C://Users/Lenovo/Desktop/archive/CUB_200_2011/classes.txt'; % 替换为您的类别映射文件路径
classMapping = readtable(classesFilePath, 'Delimiter', ' ', 'ReadVariableNames', false);
classNames = classMapping{:, 2}; % 第二列是类别名称

    imgResized = imresize(img, [224, 224]);
    
    % 归一化图像数据
    imgNormalized = double(imgResized) / 255.0;
    
    % 添加批次维度
    imgBatched = single(cat(4, imgNormalized));
    
    % 进行预测
    predictedScores = predict(net, imgBatched);
    
    % 获取预测类别索引（注意：MATLAB 索引从 1 开始）
    [~, idx] = max(predictedScores, [], 'all');
    
    % 获取预测类别名称
    predictedClass = classNames{idx};
    
   set(handles.edit8, 'String', predictedClass);
end

function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
    b=str2num(get(hObject,'String')); % 得到其中的字符串并将其转换为数字
    if isempty(b) % 判断是否为数据，若否，则将其设置为0
        set(hObject,'String','');
    end
    guidata(hObject,handles); % 更新数据
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
