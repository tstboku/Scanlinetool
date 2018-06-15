% this script enables the analysis of scanlines for rock mass and joint
% characterisations.
% you will need an input file of the same structure as "sampleinput.csv"

% Some sections of this code are based on the code presented in  Markovaara-Koivisto & Laine (2012): MATLAB script for
% analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
% This script is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation.
% This script is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHATABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.

% Thomas Strauhal (strauhal.thomas@gmx.at)
% alpS GmbH, Innsbruck, Austria


function varargout = scanlinetool(varargin)
% scanlinetool MATLAB code for scanlinetool.fig
%      scanlinetool, by itself, creates a new scanlinetool or raises the existing
%      singleton*.
%
%      H = scanlinetool returns the handle to a new scanlinetool or the handle to
%      the existing singleton*.
%
%      scanlinetool('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in scanlinetool.M with the given input arguments.
%
%      scanlinetool('Property','Value',...) creates a new scanlinetool or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scanlinetoolb_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scanlinetool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scanlinetool

% Last Modified by GUIDE v2.5 06-Mar-2018 11:23:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scanlinetool_OpeningFcn, ...
                   'gui_OutputFcn',  @scanlinetool_OutputFcn, ...
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
% --- Executes just before scanlinetool is made visible.



function scanlinetool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui03 (see VARARGIN)

% Choose default command line output for gui03
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clear global;
set(handles.panel_general,'Visible','on');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');

% --- Outputs from this function are returned to the command line.
function varargout = scanlinetool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;












%% CORE CALCULATION (TRUE BUSINESS LOGIC LIES HERE WHICH HAVE ALL THE MAJOR CALCULATIONS)
%=========================================================================================
% --- Executes on button press in pushbutton_importdata.
function pushbutton_importdata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_importdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data RESULTS sltrend slplunge sllength



%Loads information from *.csv 
%Following code section is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
% Readings of the discontinuities included in the set, their spacing and mean spacing



filepath=get(handles.edit_filelocation,'String');
    data=xlsread(filepath);
 
    dipdir=data(:,2);
    data((find(isnan(dipdir))),:)=[];
    dip=data(:,3);
    data((find(isnan(dip))),:)=[];

    
intersection=data(:,1);
dipdir=data(:,2);
dip=data(:,3);
semitracelengthabove=data(:,4);
terminationabove=data(:,5);
semitracelengthbelow=data(:,6);
terminationbelow=data(:,7);
jointroughnessfactor=data(:,8);
jointalterationfactor=data(:,9);



% Creates blank colums for angle between joint and scanline as well as the
% setname (assigned later on)
sz=size(dipdir);
slangle=zeros(sz);
setname=zeros(sz);

     
%Loads information from scanline orientation
sltrend = str2double(get (handles.edit_sltrend,'String'));
slplunge =str2double(get (handles.edit_slplunge,'String'));
slplunge=-slplunge;
sllength= str2double(get (handles.edit_sllength,'String'));



%Calculates the angle between discontinuity and scanline orientation
if dipdir>180
   dipdirn=dipdir-180;
else 
   dipdirn=dipdir+180;
end
dipn=90-dip;
slangle=90-(acosd(abs((cosd(sltrend-dipdirn)*cosd(-slplunge).*cosd(dipn))+(sind(-slplunge).*sind(dipn)))));
data(:,10)=slangle;
data(:,11)=setname;

set (handles.preview_table,'data',data);
set (handles.pushbutton_clustering,'Enable','on');
%waitbar(900 / steps)
%Calculate weighted joint density
wjd=1/sllength*sum(1./sind(slangle));
RESULTS.Rock_characteristics.wJd=wjd; %The idea of saving the results in a structure called RESULTS was adopted from Markovaara-Koivisto & Laine (2012): MATLAB script for analyzing and visualizing scanline data. Computers & Geosciences 40,185-193. doi: 10.1016/j.cageo.2011.07.010
RESULTS.Measured_area=sllength; %This parameter is important for weighting in merging several scanlines.



% --- Executes on button press in pushbutton_statistics.
function pushbutton_statistics_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS idx3

%Writes statistics into table
for n1=1:max(idx3)
setname=['Set', num2str(n1)];
corrected_meanspacings_in_sets(n1)=RESULTS.(setname).Mean_normal_spacing;
std_spacing_in_sets(n1)=RESULTS.(setname).STD_normal_spacing;
tracelengths_in_sets(n1)=RESULTS.(setname).Mean_Tracelength;
stdtracelengths_in_sets(n1)=RESULTS.(setname).STD_Tracelength;
termination_index_in_sets(n1)=RESULTS.(setname).Termination_index;
jointconditionfactor_in_sets(n1)=RESULTS.(setname).Mean_Joint_Condition_Factor;
jointroughnessfactor_in_sets(n1)=RESULTS.(setname).Mean_Joint_Roughness_Factor;
jointalterationfactor_in_sets(n1)=RESULTS.(setname).Mean_Joint_Alteration_Factor;
jointsizefactor_in_sets(n1)=RESULTS.(setname).Mean_Joint_Size_Factor;
frictionangle_in_sets(n1)=RESULTS.(setname).Mean_Peak_Friction_Angle;
stdfrictionangle_in_sets(n1)=RESULTS.(setname).STD_Peak_Friction_Angle;
tl_laslett82_in_sets(n1)=RESULTS.(setname).Mean_tracelength_Laslett1982;
tl_ph81_in_sets(n1)=RESULTS.(setname).Mean_tracelength_PriestHudson_1981;
tl_ph93_in_sets(n1)=RESULTS.(setname).Mean_tracelength_Priest_1993;
end

statistics=[corrected_meanspacings_in_sets', std_spacing_in_sets', tracelengths_in_sets',stdtracelengths_in_sets',tl_laslett82_in_sets',tl_ph81_in_sets',tl_ph93_in_sets',termination_index_in_sets', jointconditionfactor_in_sets', jointroughnessfactor_in_sets', jointalterationfactor_in_sets',jointsizefactor_in_sets',frictionangle_in_sets', stdfrictionangle_in_sets'];
set (handles.uitable_statistics,'data',statistics);

minrqd=min(RESULTS.Rock_characteristics.RQDs);
meanrqd=mean(RESULTS.Rock_characteristics.RQDs);
maxrqd=max(RESULTS.Rock_characteristics.RQDs);
allrqdvalues=length(RESULTS.Rock_characteristics.RQDs);
verypoorrqd=(sum(RESULTS.Rock_characteristics.RQDs<=25))/allrqdvalues*100;
poorrqd=(sum(RESULTS.Rock_characteristics.RQDs<=50&RESULTS.Rock_characteristics.RQDs>25))/allrqdvalues*100;
fairrqd=(sum(RESULTS.Rock_characteristics.RQDs<=75&RESULTS.Rock_characteristics.RQDs>50))/allrqdvalues*100;
goodrqd=(sum(RESULTS.Rock_characteristics.RQDs<=90&RESULTS.Rock_characteristics.RQDs>75))/allrqdvalues*100;
excellentrqd=(sum(RESULTS.Rock_characteristics.RQDs>90))/allrqdvalues*100;
set (handles.text_rqd_verypoor,'string',round(verypoorrqd));
set (handles.text_rqd_poor,'string',round(poorrqd));
set (handles.text_rqd_fair,'string',round(fairrqd));
set (handles.text_rqd_good,'string',round(goodrqd));
set (handles.text_rqd_excellent,'string',round(excellentrqd));


jv=RESULTS.Rock_characteristics.Jv;
wjd=RESULTS.Rock_characteristics.wJd;

set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','on');
set(handles.panel_volume,'Visible','off');
set(handles.panel_rmc,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');

set (handles.text_minrqd,'string',round(minrqd));
set (handles.text_meanrqd,'string',round(meanrqd));
set (handles.text_maxrqd,'string',round(maxrqd));

set (handles.text_valueJv,'string',jv);
set (handles.text_wjd,'string',wjd);





% --- Executes on button press in pushbutton_stereonet.
function pushbutton_stereonet_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_stereonet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global data

dipdir=data(:,2);
dip=data(:,3);



axesstereo=findobj('Type','axes','Tag','axes2');
axes(axesstereo); 
set(handles.axes2,'Visible','on');
set(handles.pushbutton_savefig,'Visible','on');


%Following code section (stereonet plotting and clustering) is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
    Nstereo = 50;
    cx = cos(0:pi/Nstereo:2*pi);                           % points on circle
    cy = sin(0:pi/Nstereo:2*pi);
    xh = [-1 1];                                     % horizontal axis
    yh = [0 0];
    xv = [0 0];                                      % vertical axis
    yv = [-1 1];
    axis([-1 1 -1 1]);
    axis('square');
    plot(xh,yh,'-k',xv,yv,'-k');                     %plot green axes
    axis off;
    hold on;
    plot(cx,cy,'-k');                                %plot white circle
    psi = [0:pi/Nstereo:pi];
    for i = 1:8                                      %plot great circles
        rdip = i*(pi/18);                             %at 10 deg intervals
        radip = atan(tan(rdip)*sin(psi));
        rproj = sqrt(2)*sin((pi/2 - radip)/2);
        x1 = rproj .* sin(psi);
        x2 = rproj .* (-sin(psi));
        y = rproj .* cos(psi);
        plot(x1,y,':k',x2,y,':k');
    end
    for i = 1:8                                     %plot small circles
        alpha = i*(pi/18);
        xlim = sin(alpha);
        ylim = cos(alpha);
        x = [-xlim:0.01:xlim];
        d = 1/cos(alpha);
        rd = d*sin(alpha); 				
        y0 = sqrt(rd*rd - (x .* x));
        y1 = d - y0;
        y2 = - d + y0;
        plot(x,y1,':k',x,y2,':k');
    end
    axis('square');
    hold on
    % Draws the oriented data into a Stereoplot (G. Middleton, November 1995)
    theta=(90-dipdir); % Poles are at the opposite direction to dip direction -> -180
    r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip
    [m,p]=size(data);
    % Coordinates on the strereographic projection
    for n=1:m;
    xp(n) = r(n)*cosd(theta(n));
    yp(n) = r(n)*sind(theta(n));
    plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','black','MarkerEdgeColor','black')
    end
    title('Stereoplot')

%Clustering:
set(handles.edit_clustering, 'String', 'Number of joint sets');
set(handles.popupmenu_clustersets,'Visible','on');
set(handles.pushbutton_calculate,'Visible','on');





% --- Executes on button press in pushbutton_calculate.
function pushbutton_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data dipdir dip idx3  xp yp nclu RESULTS 

dipdir=data(:,2);
dip=data(:,3);



    nclu = get(handles.popupmenu_clustersets,'Value');
    if nclu==1
    jnbarton=2;
    end
    if nclu==2
        jnbarton=3;
    end
    if nclu==3
        jnbarton=6;
    end
    if nclu==4
    jnbarton=12;
    end
    if nclu>4
        jnbarton=15;
    end
    
    set (handles.edit_rmc_jn,'string',jnbarton);
    
    
%Following code section (stereonet plotting and clustering) is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
    [m,p]=size(data);
    % Clustering is carried out on the discontinuities polevectors, which are 
    % calculated into N-matrice.
    % Calculation of the surface's normal
    % Surface is defined with two perpendicular vectors V and O
    V=zeros(1,3,m);
    for r=1:m
    % Unit vector V, which origin is at the origo, in the orientation of the 
    % dip/dip direction
        V(:,:,r)=[cosd(90-dipdir(r))*cosd(dip(r)) sind(90-dipdir(r))*cosd(dip(r)) -sind(dip(r))];
    % Unit vector O, which origin is at the origo, which is perpendicular to 
    % the dip direction and is horizontal
        O(:,:,r)=[-sind(90-dipdir(r)) cosd(90-dipdir(r)) 0];
    end
    % Normal vector N of the surface defined with the vectors V and O
    N=zeros(1,3,m);
    for k=1:m 
        N(:,:,k)=cross(V(:,:,k),O(:,:,k));
    end
    % Pole vector is counter vector of the normal vector, and it points
    % downwards 
    N=-N;
    N=N.*sign(N);
    
    % Asks to click on the approximate centroids of the clusters nclu times
    set(handles.edit_setclusters, 'String', 'Left click with mouse on the approximate centroids of the clusters!');
    % disp('Left click with mouse on the approximate centroids of the clusters.')
    %%%figure(Stereoplot);
    %%%movegui(Stereoplot,'northeast');
    [xclu yclu]=ginput(nclu);

    % Middlepoints of the clusters
    % Clusters' coordinates are transformed into Dip/Dipdirection form
    cludipdir=zeros(1,nclu);
    cludip=zeros(1,nclu);
    for s=1:nclu
        if xclu(s)>0 & yclu(s)>0
            clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
            cludipdir(s)=270-clutheta(s)*180/pi;
            cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
            cludip(s)=abs(cludipr(s)*180/pi);
        elseif xclu(s)>0 & yclu(s)<0
            clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
            cludipdir(s)=270+clutheta(s)*180/pi;
            cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
            cludip(s)=abs(cludipr(s)*180/pi);
        elseif xclu(s)<0 & yclu(s)<0
            clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
            cludipdir(s)=90-clutheta(s)*180/pi;
            cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
            cludip(s)=abs(cludipr(s)*180/pi);
        else
            clutheta(s)=asin(sqrt(1/((xclu(s)/yclu(s))^2+1)));
            cludipdir(s)=90+clutheta(s)*180/pi;
            cludipr(s)=2*atan(xclu(s)/cos(clutheta(s)));
            cludip(s)=abs(cludipr(s)*180/pi);
        end
    end

    % Dip/Dip direction forms is transformed into polevectors, which are used in the clustering
    % Calculation of the normal vectors N
    % Defining the surface with two perpendicular vectors V and O
    V=zeros(1,3,nclu);
    for r=1:nclu
        % Unit vector V, which origin is at the origo, in the orientation of the 
        % dip/dip direction
        V(:,:,r)=[cosd(90-cludipdir(r))*cosd(cludip(r)) sind(90-cludipdir(r))*cosd(cludip(r)) -sind(cludip(r))];
        % Unit vector O, which origin is at the origo, which is perpendicular to 
        % the dip diprecion and is horizontal
        O(:,:,r)=[-sind(90-cludipdir(r)) cosd(90-cludipdir(r)) 0];
    end
    % Normal vector N of the surface defined with the vectors V and O
    Nclu=zeros(1,3,s);
    for k=1:s 
        Nclu(:,:,k)=cross(V(:,:,k),O(:,:,k));
    end
    %Pole vector is counter vector of the normal vector, and it points
    %downwards 
    Nclu=-Nclu;   %Modification of the orignal code
    Nclu=Nclu.*sign(Nclu); %Modification of the orignal code

    % Iteration of the contents of the clusters: Initial centroids of the
    % clustres are given and the new ones are calculated, observations are moved
    % to the clusters and sum of distanced of centroids and the observations 
    % within the cluster are calculated
    % First iteration
    previous_clustercosum=[];
    % Observation are divided into the closest clusters
    arcdistance=zeros(s,m);
    for s=1:nclu
        for k=1:m
            arcdistance(s,k)=(acos(dot(N(:,:,k),Nclu(:,:,s))))^2; % distance between the observation and the centroids, nclu defines rows and number of observations columns
        end
    end
    for s=1:nclu
        for k=1:m
            rownumber(k)=find(arcdistance(:,k)==min(arcdistance(:,k))); % Searches the rownumber of the centroid, which has the smallest distance to the observation point.
            idx3(k)=rownumber(k); % Observation gets the rownumber as its cluster number.
        end
    end
    idx3first=idx3;
    % Calculates the sum of arcdistances within the cluster
    for s=1:nclu
        clustersum(s)=sum(arcdistance(find(idx3==s)));
    end
    % Sums up all the cluster's sums together
    cosum=sum(clustersum);

    % Iterates the cluster cosum as long as it does not get any smaller

    % Second and later iterations    
    while (isempty(previous_clustercosum) | cosum < previous_clustercosum)
        % Calculates new centroids of the clusters
        % Calculates sumvector of the clusters' polevectors
        previous_clustercosum=cosum;
        for s=1:nclu
            sum_cluster(:,s)=[sum(N(1,1,find(idx3==s))) sum(N(1,2,find(idx3==s))) sum(N(1,3,find(idx3==s)))]; %Finds the rows with the same idc3 and sums them up
            % Length of the polevectors' sum
            L(s)=sqrt(sum(N(1,1,find(idx3==s)))^2+sum(N(1,2,find(idx3==s)))^2+sum(N(1,3,find(idx3==s)))^2);
            % Sumvector is normalised into a unitvector
            new_Nclu(:,s)=sum_cluster(:,s)/L(s);
        end

        % Points are re-divided according to the closest new cluster centroids
        arcdistance=zeros(s,m);
        for s=1:nclu
            for k=1:m
                arcdistance(s,k)=(acos(dot(N(:,:,k),new_Nclu(:,s))))^2; % distance between the observation and the centroids, nclu defines rows and number of observations columns
            end
        end
        for s=1:nclu
            for k=1:m
            rownumber(k)=find(arcdistance(:,k)==min(arcdistance(:,k))); % Searches the rownumber of the centroid, which has the smallest distance to the observation point.
            end
        end

        idx3=rownumber; % Observation gets the rownumber as its cluster number.
        idx3second=idx3;
        % Calculates the sum of arcdistances within the cluster
        for s=1:nclu
            clustersum(s)=sum(arcdistance(find(idx3==s)));
        end
        % Sums up all the clustre's sums together
        cosum=sum(clustersum);
    end
    
    %Results of clustering
    close(findobj('type','figure','name','Stereoplot'))
    close(findobj('type','figure','name','Stereoplot Clusters'))
    Nstereo = 50;
    cx = cos(0:pi/Nstereo:2*pi);                           % points on circle
    cy = sin(0:pi/Nstereo:2*pi);
    xh = [-1 1];                                     % horizontal axis
    yh = [0 0];
    xv = [0 0];                                      % vertical axis
    yv = [-1 1];
    axis([-1 1 -1 1]);
    axis('square');
    plot(xh,yh,'-k',xv,yv,'-k');                     %plot green axes
    axis off;
    hold on;
    plot(cx,cy,'-k');                                %plot white circle
    psi = [0:pi/Nstereo:pi];
    for i = 1:8                                      %plot great circles
        rdip = i*(pi/18);                             %at 10 deg intervals
        radip = atan(tan(rdip)*sin(psi));
        rproj = sqrt(2)*sin((pi/2 - radip)/2);
        x1 = rproj .* sin(psi);
        x2 = rproj .* (-sin(psi));
        y = rproj .* cos(psi);
        plot(x1,y,':k',x2,y,':k');
    end
    for i = 1:8                                     %plot small circles
        alpha = i*(pi/18);
        xlim = sin(alpha);
        ylim = cos(alpha);
        x = [-xlim:0.01:xlim];
        d = 1/cos(alpha);
        rd = d*sin(alpha); 				
        y0 = sqrt(rd*rd - (x .* x));
        y1 = d - y0;
        y2 = - d + y0;
        plot(x,y1,':k',x,y2,':k');
    end
    axis('square');

    % Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
    theta=(90-dipdir); % Poles are at the opposite direction to dip direction
    r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip
    % Coordinates on the strereographic projection
    [m,p]=size(data);
    for n=1:m;
        xp(n) = r(n)*cosd(theta(n));
        yp(n) = r(n)*sind(theta(n));

        % Stereographic projection with colourcodes according to index idx3
        if idx3(n)==1 
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0 0.45 0.75','MarkerEdgeColor','black')
        elseif idx3(n)==2
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','red','MarkerEdgeColor','black')
        elseif idx3(n)==3
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','green','MarkerEdgeColor','black')
        elseif idx3(n)==4
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','magenta','MarkerEdgeColor','black')
        elseif idx3(n)==5
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.87 0.49 0','MarkerEdgeColor','black')
        elseif idx3(n)==6
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.6 0.2 0','MarkerEdgeColor','black')
        else 
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','black','MarkerEdgeColor','black')
        end
    end
    title('Stereoplot with clusters')
        
    % Asks if the user is satisfied with the result of the clustering
    % If not, the user may change the clustering in an interactive window
    Ready = [];
    button=1;
    set(handles.edit_clusterok, 'String', 'Do you want to change the results of the automated clustering:');
   
   
set(handles.popupmenu_changeclusters,'Visible','on');
set(handles.pushbutton_clusterok,'Visible','on');




% --- Executes on button press in pushbutton_clusterok.
function pushbutton_clusterok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clusterok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RESULTS idx3 data xp yp slplunge sltrend sllength

dipdir=data(:,2);
dip=data(:,3);

%Following code section (stereonet plotting and clustering) is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
    Ready = [];
    button=1;
    answer = get(handles.popupmenu_changeclusters,'string');
    V_Index = get(handles.popupmenu_changeclusters,'Value');
    while (V_Index == 1 & isempty(Ready) & button==1)
        set(handles.edit_statusmanual, 'String', 'Left click with the mouse on the individual poles as many times as it takes to change the assigned cluster! Right click into the stereoplot to finish!');
        [xcho, ycho, button]=ginput(1);
        % Calculates the minimum distance between the indicated point and the data
        % points and searches the closes data point, gives index to the row and 
        % changes the cluster group into the next one or to the first one.
        % Refreshes the figure.
        distance=sqrt((xp-xcho).^2+(yp-ycho).^2);
        cho=find(distance==min(distance));
        % Observation is indicated to the next cluster
        if (idx3(cho)<max(idx3) & button==1)
            idx3(cho)=idx3(cho)+1;
        elseif button==1
            idx3(cho)=1;
        end
        close(findobj('type','figure','name','Stereoplot Clusters'))
        Nstereo = 50;
        cx = cos(0:pi/Nstereo:2*pi);                           % points on circle
        cy = sin(0:pi/Nstereo:2*pi);
        xh = [-1 1];                                     % horizontal axis
        yh = [0 0];
        xv = [0 0];                                      % vertical axis
        yv = [-1 1];
        axis([-1 1 -1 1]);
        axis('square');
        plot(xh,yh,'-k',xv,yv,'-k');                     %plot green axes
        axis off;
        hold on;
        plot(cx,cy,'-k');                                %plot white circle
        psi = [0:pi/Nstereo:pi];
        for i = 1:8                                      %plot great circles
            rdip = i*(pi/18);                             %at 10 deg intervals
            radip = atan(tan(rdip)*sin(psi));
            rproj = sqrt(2)*sin((pi/2 - radip)/2);
            x1 = rproj .* sin(psi);
            x2 = rproj .* (-sin(psi));
            y = rproj .* cos(psi);
            plot(x1,y,':k',x2,y,':k');
        end
        for i = 1:8                                     %plot small circles
            alpha = i*(pi/18);
            xlim = sin(alpha);
            ylim = cos(alpha);
            x = [-xlim:0.01:xlim];
            d = 1/cos(alpha);
            rd = d*sin(alpha); 				
            y0 = sqrt(rd*rd - (x .* x));
            y1 = d - y0;
            y2 = - d + y0;
            plot(x,y1,':k',x,y2,':k');
        end
        axis('square');

        % Draws the oriented data into Schmidt's net (G. Middleton, November 1995)
        theta=(90-dipdir); % Poles are at the opposite direction to dip direction
        r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip

        % Coordinates on the strereographic projection
        [m,p]=size(data);
        for n=1:m;
            xp(n) = r(n)*cosd(theta(n));
            yp(n) = r(n)*sind(theta(n));
            % Stereographic projection with colourcodes according to index idx3
            if idx3(n)==1 
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0 0.45 0.75','MarkerEdgeColor','black')
            elseif idx3(n)==2
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','red','MarkerEdgeColor','black')
            elseif idx3(n)==3
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','green','MarkerEdgeColor','black')
            elseif idx3(n)==4
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','magenta','MarkerEdgeColor','black')
            elseif idx3(n)==5
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.87 0.49 0','MarkerEdgeColor','black')
            elseif idx3(n)==6
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.6 0.2 0','MarkerEdgeColor','black')
            else 
                plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','black','MarkerEdgeColor','black')
            end
        end
        title('Stereoplot with clusters')
    end
   
    data(:,11)=idx3;
set (handles.preview_table,'data',data);


%THIS SECTION IS MOST PROBABLY NOT NEEDED!! DELETE IT MAYBE???
theta=(dipdir+180); %Polevectors are at the opposite direction to the dip direction -> -180
X3=[(dip-45)/45]'; %z-coordinate on the sterepgraphic sphere
r=cos(asin(X3)); % Poles are perpendicular to the dip
[m,p]=size(data);
for n=1:m;
X1(n) = r(n)*sind(theta(n)); % x-coordinate
X2(n) = r(n)*cosd(theta(n)); % y-coordinate
end
%UNTIL HERE!!!!!

%Following code section is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010

    %Calculates the centroids of the clusters
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
   
    %Calculates products from direction cosines
    for n1=1:max(idx3)
        xx(n1)=sum(X1b(find(idx3==n1)).^2);
        xy(n1)=sum(X1b(find(idx3==n1)).*X2b(find(idx3==n1)));
        xz(n1)=sum(X1b(find(idx3==n1)).*X3b(find(idx3==n1)));
        yy(n1)=sum(X2b(find(idx3==n1)).^2);
        yz(n1)=sum(X2b(find(idx3==n1)).*X3b(find(idx3==n1)));
        zz(n1)=sum(X3b(find(idx3==n1)).^2);
        %Orientation matrix
        S(:,:,n1)=[xx(n1) xy(n1) xz(n1); xy(n1) yy(n1) yz(n1);xz(n1) yz(n1) zz(n1)];
        %Orientation matrix normalized with number of obervation in the cluster
        Sn=S./sum(idx3==n1);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,n1));
        B(:,:,n1)=Bi;
        K(:,:,n1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,n1)=max(K(:,:,n1));
        n2=find(Kmax(:,:,n1)==max(Kmax(:,:,n1)));
        mean_vector(:,n1)= B(:,n2,n1);
        %Eigenvector's x,y and z coordinates 
        xeig(n1)=mean_vector(1,n1);
        yeig(n1)=mean_vector(2,n1);
        zeig(n1)=mean_vector(3,n1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig(n1)>0 & yeig(n1)>0
            dipeigv(n1)=90-abs((asin(zeig(n1))))*(180/pi);
            dipdireigv(n1)=abs((atan(xeig(n1)/yeig(n1))))*180/pi;
        elseif xeig(n1)>0 & yeig(n1)<0
            dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
            dipdireigv(n1)=180-abs((atan(xeig(n1)/yeig(n1))))*180/pi;
        elseif xeig(n1)<0 & yeig(n1)<0
            dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
            dipdireigv(n1)=180+abs((atan(xeig(n1)/yeig(n1))))*180/pi;
        else
            dipeigv(n1)=90-abs(asin(zeig(n1)))*180/pi;
            dipdireigv(n1)=360-abs((atan(xeig(n1)/yeig(n1))))*180/pi;
        end
    end

%Following structure of saving the most important results into a
%structure RESULTS was adopted from Markovaara-Koivisto & Laine (2012): MATLAB script for analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
    for n1=1:max(idx3)
        RESULTS.(genvarname(['Set', num2str(n1)])).Number_of_joints_in_set=sum(idx3==n1);
        RESULTS.(genvarname(['Set', num2str(n1)])).Dip=dipeigv(n1)    ; 
        RESULTS.(genvarname(['Set', num2str(n1)])).Dipdir=dipdireigv(n1);
        %genvarname(['Set', num2str(n1)])=sum(idx3==n1)
    end

    
%Following code section is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
    %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(max(idx3),1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
        if idx3(n2)==1 
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
         elseif idx3(n2)==2
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(2)>180
                dipdireigvang(2)=dipdireigv(2)-180;
            else 
                dipdireigvang(2)=dipdireigv(2)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(2)=90-dipeigv(2);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(2)).*cosd(dipang(n2)).*cosd(dipeigvang(2)))+(sind(dipang(n2)).*sind(dipeigvang(2)))));
         elseif idx3(n2)==3
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(3)>180
                dipdireigvang(3)=dipdireigv(3)-180;
            else 
                dipdireigvang(3)=dipdireigv(3)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(3)=90-dipeigv(3);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(3)).*cosd(dipang(n2)).*cosd(dipeigvang(3)))+(sind(dipang(n2)).*sind(dipeigvang(3)))));
         elseif idx3(n2)==4
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(4)>180
                dipdireigvang(4)=dipdireigv(4)-180;
            else 
                dipdireigvang(4)=dipdireigv(4)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(4)=90-dipeigv(4);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(4)).*cosd(dipang(n2)).*cosd(dipeigvang(4)))+(sind(dipang(n2)).*sind(dipeigvang(4)))));
         elseif idx3(n2)==5
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(5)>180
                dipdireigvang(5)=dipdireigv(5)-180;
            else 
                dipdireigvang(5)=dipdireigv(5)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(5)=90-dipeigv(5);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(5)).*cosd(dipang(n2)).*cosd(dipeigvang(5)))+(sind(dipang(n2)).*sind(dipeigvang(5)))));
         elseif idx3(n2)==6
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(6)>180
                dipdireigvang(6)=dipdireigv(6)-180;
            else 
                dipdireigvang(6)=dipdireigv(6)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(6)=90-dipeigv(6);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(6)).*cosd(dipang(n2)).*cosd(dipeigvang(6)))+(sind(dipang(n2)).*sind(dipeigvang(6)))));
         else 
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end
       
            if dipdireigv(7)>180
                dipdireigvang(7)=dipdireigv(7)-180;
            else 
            dipdireigvang(7)=dipdireigv(7)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(7)=90-dipeigv(7);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(7)).*cosd(dipang(n2)).*cosd(dipeigvang(7)))+(sind(dipang(n2)).*sind(dipeigvang(7)))));
         end
   end
      
for n3=1:max(idx3)
	angleclu_set=angleclu(find(idx3==n3));
	stdangleclu(n3)=std(angleclu_set);
    RESULTS.(genvarname(['Set', num2str(n1)])).STD_Orientation=stdangleclu(n3);
end   


%Calculates the number of joints per cluster set
SumCluster1=sum(idx3==1);
SumCluster2=sum(idx3==2);
SumCluster3=sum(idx3==3);
SumCluster4=sum(idx3==4);
SumCluster5=sum(idx3==5);
SumCluster6=sum(idx3==6);
SumCluster7=sum(idx3==7);
SumCluster=[SumCluster1, SumCluster2, SumCluster3, SumCluster4, SumCluster5, SumCluster6, SumCluster7];
SumCluster(SumCluster==0)=[];

%Plots clustering results as a table
results=[dipdireigv',dipeigv',SumCluster', stdangleclu];
set (handles.uitable_clusterresults,'data',results);




%Following code section is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
% Draws a stereographic projection with the mean attitudes of the clusters
    Nstereo = 50;
    cx = cos(0:pi/Nstereo:2*pi);                           % points on circle
    cy = sin(0:pi/Nstereo:2*pi);
    xh = [-1 1];                                     % horizontal axis
    yh = [0 0];
    xv = [0 0];                                      % vertical axis
    yv = [-1 1];
    axis([-1 1 -1 1]);
    axis('square');
    plot(xh,yh,'-k',xv,yv,'-k');                     %plot green axes
    axis off;
    hold on;
    plot(cx,cy,'-k');                                %plot white circle
    psi = [0:pi/Nstereo:pi];
    for i = 1:8                                      %plot great circles
        rdip = i*(pi/18);                             %at 10 deg intervals
        radip = atan(tan(rdip)*sin(psi));
        rproj = sqrt(2)*sin((pi/2 - radip)/2);
        x1 = rproj .* sin(psi);
        x2 = rproj .* (-sin(psi));
        y = rproj .* cos(psi);
        plot(x1,y,':k',x2,y,':k');
    end
    for i = 1:8                                     %plot small circles
        alpha = i*(pi/18);
        xlim = sin(alpha);
        ylim = cos(alpha);
        x = [-xlim:0.01:xlim];
        d = 1/cos(alpha);
        rd = d*sin(alpha); 				
        y0 = sqrt(rd*rd - (x .* x));
        y1 = d - y0;
        y2 = - d + y0;
        plot(x,y1,':k',x,y2,':k');
    end
    axis('square');
    %Draws the oriented data into a Stereoplot (G. Middleton, November 1995)
    theta=(90-dipdir); % Poles are at the opposite direction to dip direction
    r=sqrt(2)*sind((90-dip-90)/2); % Poles are perpendicular to the dip
    %Coordinates on the strereographic projection
    [m,p]=size(data);
    for n=1:m;
        xp(n) = r(n)*cosd(theta(n));
        yp(n) = r(n)*sind(theta(n));
        %Stereographic projection with colourcodes according to index idx3
        if idx3(n)==1 
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0 0.45 0.75','MarkerEdgeColor','black')
        elseif idx3(n)==2
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','red','MarkerEdgeColor','black')
        elseif idx3(n)==3
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','green','MarkerEdgeColor','black')
        elseif idx3(n)==4
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','magenta','MarkerEdgeColor','black')
        elseif idx3(n)==5
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.87 0.49 0','MarkerEdgeColor','black')
        elseif idx3(n)==6
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','0.6 0.2 0','MarkerEdgeColor','black')
        else 
            plot(xp(n),yp(n),'o','Markersize',4,'MarkerFaceColor','black','MarkerEdgeColor','black')
        end
    end
    title('Stereoplot with clusters and centroids')

    thetaeig=(90-dipdireigv);
    reig=sqrt(2)*sind((90-dipeigv-90)/2);
    for n=1:max(idx3)
        xeigp(n) = reig(n)*cosd(thetaeig(n));
        yeigp(n) = reig(n)*sind(thetaeig(n));
    end
    plot(xeigp, yeigp, 'x','Markersize',12, 'MarkerFaceColor','black', 'linewidth',3, 'Color', 'black')
 
    
set (handles.pushbutton_statistics,'Enable','on');

set (handles.pushbutton_vol,'Enable','on');
 
% General statistical calculations
intersection=data(:,1);
dipdir=data(:,2);
dip=data(:,3);
semitracelengthabove=data(:,4);
terminationabove=data(:,5);
semitracelengthbelow=data(:,6);
terminationbelow=data(:,7);
jointroughnessfactor=data(:,8);
jointalterationfactor=data(:,9);

%Calculation of the joint size factor according to Palmstroem (2000) but
%modified for semi-tracelengths:
above(:,1)=semitracelengthabove;
above(:,2)=terminationabove;
abovenonzeros = above(all(above,2),:);

jointsemilengthfactorsabove=abovenonzeros;

for i=1:length(jointsemilengthfactorsabove)
    if jointsemilengthfactorsabove(i,1) < 0.5
        jointsemilengthfactorsabove(i,1)=2;
    elseif jointsemilengthfactorsabove(i,1)>=0.5 & jointsemilengthfactorsabove(i,1)<5
        jointsemilengthfactorsabove(i,1)=1;
    elseif jointsemilengthfactorsabove(i,1)>=5 & jointsemilengthfactorsabove(i,1)<15    
        jointsemilengthfactorsabove(i,1)=0.75;
    else % >=15
        jointsemilengthfactorsabove(i,1)=0.5;
    end 
    jointsemilengthfactorsabove(i,1)=jointsemilengthfactorsabove(i,1)*jointsemilengthfactorsabove(i,2);
end

below(:,1)=semitracelengthbelow;
below(:,2)=terminationbelow;
belownonzeros = below(all(below,2),:);

jointsemilengthfactorsbelow=belownonzeros;
for i=1:length(jointsemilengthfactorsbelow)
    if jointsemilengthfactorsbelow(i,1) < 0.5
        jointsemilengthfactorsbelow(i,1)=2;
    elseif jointsemilengthfactorsbelow(i,1)>=0.5 & jointsemilengthfactorsbelow(i,1)<5
        jointsemilengthfactorsbelow(i,1)=1;
    elseif jointsemilengthfactorsbelow(i,1)>=5 & jointsemilengthfactorsbelow(i,1)<15    
        jointsemilengthfactorsbelow(i,1)=0.75;
    else % >=15
        jointsemilengthfactorsbelow(i,1)=0.5;
    end 
    jointsemilengthfactorsbelow(i,1)=jointsemilengthfactorsbelow(i,1)*jointsemilengthfactorsbelow(i,2);
end
alljointsemisizefactors=vertcat(jointsemilengthfactorsbelow(:,1),jointsemilengthfactorsabove(:,1));
meanjointsizefactor=mean(alljointsemisizefactors);

bothterminations=vertcat(belownonzeros(:,2),abovenonzeros(:,2));
meantermination=mean(bothterminations);
bothsemitraces=vertcat(belownonzeros(:,1),abovenonzeros(:,1));
meanfulltracelength=2*mean(bothsemitraces);
stdfulltracelength=std(2.*bothsemitraces);

meanjointroughnessfactor=mean(jointroughnessfactor);
meanjointalterationfactor=mean(jointalterationfactor);
friction_angles=atand(jointroughnessfactor./jointalterationfactor);
meanfrictionangles=mean(friction_angles);
meanjointconditionfactor=meanjointroughnessfactor*meanjointsizefactor/meanjointalterationfactor;

RESULTS.Joint_characteristics.Mean_Termination=meantermination;
RESULTS.Joint_characteristics.Terminations=bothterminations';
RESULTS.Joint_characteristics.Mean_Tracelength=meanfulltracelength;
RESULTS.Joint_characteristics.STD_Tracelength=stdfulltracelength;
RESULTS.Joint_characteristics.Tracelengths=2.*bothsemitraces';
RESULTS.Joint_characteristics.Mean_Joint_Size_Factor=meanjointsizefactor;
RESULTS.Joint_characteristics.Mean_Joint_Alteration_Factor=meanjointalterationfactor;
RESULTS.Joint_characteristics.Joint_Alteration_Factor=jointalterationfactor';
RESULTS.Joint_characteristics.Mean_Joint_Roughness_Factor=meanjointroughnessfactor;
RESULTS.Joint_characteristics.Joint_Roughness_Factor=jointroughnessfactor';
RESULTS.Joint_characteristics.Mean_Joint_Condition_Factor=meanjointconditionfactor;
RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles=friction_angles';
RESULTS.Joint_characteristics.Mean_Peak_Joint_Friction_Angles=meanfrictionangles;
RESULTS.Joint_characteristics.Joint_Size_Factor=alljointsemisizefactors';
RESULTS.Joint_characteristics.Mean_Joint_Size_Factor=meanjointsizefactor;



%Loads information from scanline
for n1=1:max(idx3)
   semitracelengthabove_set=semitracelengthabove(find(idx3==n1));
   terminationabove_set=terminationabove(find(idx3==n1));
   semitracelengthbelow_set=semitracelengthbelow(find(idx3==n1));
   terminationbelow_set=terminationbelow(find(idx3==n1));
   jointroughnessfactor_set=jointroughnessfactor(find(idx3==n1));
   jointalterationfactor_set=jointalterationfactor(find(idx3==n1));
   peakfrictionangle_set=atand(jointroughnessfactor_set./jointalterationfactor_set);
   meanpeakfrictionangle_set=mean(peakfrictionangle_set);
   stdpeakfrictionangle_set=std(peakfrictionangle_set);
   Dips_set=dip(find(idx3==n1));
   Dipdirs_set=dipdir(find(idx3==n1));
   RESULTS.(genvarname(['Set', num2str(n1)])).Joint_Roughness_Factors=jointroughnessfactor_set';
   RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Joint_Roughness_Factor=mean(jointroughnessfactor_set);
   RESULTS.(genvarname(['Set', num2str(n1)])).STD_Joint_Roughness_Factor=std(jointroughnessfactor_set);
   RESULTS.(genvarname(['Set', num2str(n1)])).Joint_Alteration_Factors=jointalterationfactor_set';
   RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Joint_Alteration_Factor=mean(jointalterationfactor_set);
   RESULTS.(genvarname(['Set', num2str(n1)])).STD_Joint_Alteration_Factor=std(jointalterationfactor_set);
   RESULTS.(genvarname(['Set', num2str(n1)])).Peak_Friction_Angles=peakfrictionangle_set';
   RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Peak_Friction_Angle=meanpeakfrictionangle_set;
   RESULTS.(genvarname(['Set', num2str(n1)])).STD_Peak_Friction_Angle=stdpeakfrictionangle_set;  
   RESULTS.(genvarname(['Set', num2str(n1)])).Dips=Dips_set';
   RESULTS.(genvarname(['Set', num2str(n1)])).Dipdirs=Dipdirs_set';
   above_set(:,1)=semitracelengthabove_set;
   above_set(:,2)=terminationabove_set;
   abovenonzeros_set = above_set(all(above_set,2),:);

   below_set(:,1)=semitracelengthbelow_set;
   below_set(:,2)=terminationbelow_set;
   belownonzeros_set = below_set(all(below_set,2),:);
   
  
%Mean tracelength estimation according to Laslett (1982)
bothvisible=0;
onevisible=0;
nonevisible=0;
for i=1:length(terminationabove_set)
    if terminationabove_set(i)>0.5 & terminationbelow_set(i)>0.5
        currenttrace=semitracelengthabove_set(i)+semitracelengthbelow_set(i);
        bothvisible=vertcat(bothvisible,currenttrace);
    elseif terminationabove_set(i)>0.5 & terminationbelow_set(i)==0
        currenttrace=semitracelengthabove_set(i)+semitracelengthbelow_set(i);
        onevisible=vertcat(onevisible,currenttrace);
    elseif terminationabove_set(i)==0 & terminationbelow_set(i)>0.5
        currenttrace=semitracelengthabove_set(i)+semitracelengthbelow_set(i);
        onevisible=vertcat(onevisible,currenttrace);
    else %terminationabove_set(i)==0 & terminationbelow_set(i)==0;
        currenttrace=semitracelengthabove_set(i)+semitracelengthbelow_set(i);
        nonevisible=vertcat(nonevisible,currenttrace);    
    end  
end
bothvisible(bothvisible==0) = [];
onevisible(onevisible==0) = [];
nonevisible(nonevisible==0)=[];
meantracelaslett=(sum(bothvisible)+sum(onevisible)+sum(nonevisible))/(2*length(bothvisible)+length(onevisible));
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_tracelength_Laslett1982=meantracelaslett;
clear meantracelaslett nonevisible onevisible bothvisible currenttrace

   
%Mean tracelength estimation according to Priest (1993)
semitracelength_set=(vertcat(semitracelengthabove_set,semitracelengthbelow_set));
minbin=(min(semitracelength_set))*0.9;
maxbin=(max(semitracelength_set))*1.1;
steps=(maxbin-minbin)/14;
binstrace=(minbin:steps:maxbin);
histsemitraces = (hist(semitracelength_set,binstrace));
fitexp = fit(binstrace.',histsemitraces.','exp1');
h0 = fitexp(0)./(length(semitracelength_set).*(binstrace(2)-binstrace(1)));
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_tracelength_Priest_1993 = 1./h0;
clear semitracelength_set histsemitraces fitexp h0

%Mean tracelength estimation according to Priest & Hudson (1981)
for i=1:length(semitracelengthabove_set)
fulltracelength_set(i)=semitracelengthabove_set(i)+semitracelengthbelow_set(i);
end
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_tracelength_PriestHudson_1981 = (mean(fulltracelength_set))/2;




   %BLOCK
jointsemilengthfactorsabove_set=abovenonzeros_set;
[nrrows nrcol]=size(jointsemilengthfactorsabove_set);
for i=1:nrrows
    if jointsemilengthfactorsabove_set(i,1) < 0.5
        jointsemilengthfactorsabove_set(i,1)=2;
    elseif jointsemilengthfactorsabove_set(i,1)>=0.5 & jointsemilengthfactorsabove_set(i,1)<5
        jointsemilengthfactorsabove_set(i,1)=1;
    elseif jointsemilengthfactorsabove_set(i,1)>=5 & jointsemilengthfactorsabove_set(i,1)<15    
        jointsemilengthfactorsabove_set(i,1)=0.75;
    else % >=15
        jointsemilengthfactorsabove_set(i,1)=0.5;
    end 
    jointsemilengthfactorsabove_set(i,1)=jointsemilengthfactorsabove_set(i,1)*jointsemilengthfactorsabove_set(i,2);
end



jointsemilengthfactorsbelow_set=belownonzeros_set;
[nrrows nrcol]=size(jointsemilengthfactorsbelow_set);
for i=1:nrrows
    if jointsemilengthfactorsbelow_set(i,1) < 0.5
        jointsemilengthfactorsbelow_set(i,1)=2;
    elseif jointsemilengthfactorsbelow_set(i,1)>=0.5 & jointsemilengthfactorsbelow_set(i,1)<5
        jointsemilengthfactorsbelow_set(i,1)=1;
    elseif jointsemilengthfactorsbelow_set(i,1)>=5 & jointsemilengthfactorsbelow_set(i,1)<15    
        jointsemilengthfactorsbelow_set(i,1)=0.75;
    else % >=15
        jointsemilengthfactorsbelow_set(i,1)=0.5;
    end 
    jointsemilengthfactorsbelow_set(i,1)=jointsemilengthfactorsbelow_set(i,1)*jointsemilengthfactorsbelow_set(i,2);
end
alljointsemisizefactors_set=vertcat(jointsemilengthfactorsbelow_set(:,1),jointsemilengthfactorsabove_set(:,1));
meanjointsizefactor_set=mean(alljointsemisizefactors_set);
stdjointsizefactor_set=std(alljointsemisizefactors_set);

bothterminations_set=vertcat(belownonzeros_set(:,2),abovenonzeros_set(:,2));
meantermination_set=mean(bothterminations_set);
bothsemitraces_set=vertcat(belownonzeros_set(:,1),abovenonzeros_set(:,1));
meanfulltracelength_set=2*mean(bothsemitraces_set);
stdfulltracelength_set=std(2.*bothsemitraces_set);
jointconditionfactor_set=meanjointsizefactor_set*(mean(jointroughnessfactor_set))/(mean(jointalterationfactor_set));

RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Termination=meantermination_set;
RESULTS.(genvarname(['Set', num2str(n1)])).Terminations=bothterminations_set';
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Tracelength=meanfulltracelength_set;
RESULTS.(genvarname(['Set', num2str(n1)])).STD_Tracelength=stdfulltracelength_set;
RESULTS.(genvarname(['Set', num2str(n1)])).Tracelengths=2.*bothsemitraces_set';
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Joint_Size_Factor=meanjointsizefactor_set;
RESULTS.(genvarname(['Set', num2str(n1)])).STD_Joint_Size_Factor=stdjointsizefactor_set;
RESULTS.(genvarname(['Set', num2str(n1)])).Joint_Size_Factor=alljointsemisizefactors_set';
RESULTS.(genvarname(['Set', num2str(n1)])).Mean_Joint_Condition_Factor=jointconditionfactor_set;

 
   
   %Termination index acc. ISRM 1978 (not necessary for IBSD or rock mass
   %characterisation index)
   termination_intact_above_set=numel(terminationabove_set(terminationabove_set==2));
   termination_intact_below_set=numel(terminationbelow_set(terminationbelow_set==2));
   termination_intact_set=termination_intact_above_set+termination_intact_below_set;
   numberofelementstermination=numel(vertcat(terminationabove_set, terminationbelow_set));
   termination_index_set=termination_intact_set*100/numberofelementstermination;
   RESULTS.(genvarname(['Set', num2str(n1)])).Termination_index=termination_index_set';
   

   %Creates paramaters applicable for statistical table:
   v = genvarname(['meanfulltracelength_set', num2str(n1)]);
   tracelengths_in_sets(n1)=meanfulltracelength_set;
   v = genvarname(['stdfulltracelength_set', num2str(n1)]);
   stdtracelengths_in_sets(n1)=stdfulltracelength_set;
   v = genvarname(['termination_index_set', num2str(n1)]);
   termination_index_in_sets(n1)=termination_index_set;
   v = genvarname(['meanpeakfrictionangle_set', num2str(n1)]);
   frictionangle_in_sets(n1)=meanpeakfrictionangle_set;
   v = genvarname(['termination_index_set', num2str(n1)]);
   stdfrictionangle_in_sets(n1)=stdpeakfrictionangle_set;
   v = genvarname(['jointroughnessfactor_set', num2str(n1)]);
   jointroughnessfactor_in_sets(n1)=mean(jointroughnessfactor_set);
   v = genvarname(['jointalterationfactor_set', num2str(n1)]);
   jointalterationfactor_in_sets(n1)=mean(jointalterationfactor_set);   
   v = genvarname(['meanjointsizefactor_set', num2str(n1)]);
   jointsizefactor_in_sets(n1)=meanjointsizefactor_set;    
   v = genvarname(['jointconditionfactor_set', num2str(n1)]);
   jointconditionfactor_in_sets(n1)=jointconditionfactor_set;    
   
clear above_set below_set    
end

     
%Following code section is modifified after Markovaara-Koivisto & Laine (2012): MATLAB script for
%analyzing and visualizing scanline data. Computers & Geosciences 40,
%185-193. doi: 10.1016/j.cageo.2011.07.010
% Discontinuity spacing in sets
   intersectionb=data(:,1);
   intersection=[0; intersectionb];
   %Correction factor is calculated from the angle between the scanline and
   %the orientation of the disconintuity set.
   %Normal vector for a discontinuity set
   %Surface is defined with two perpendicular vectors V and O
   %Size of V
   V=zeros(1,3,max(idx3));
   for r=1:max(idx3)
      %Unit vector V, which origin is at the origo, in the orientation of the 
      %dip/dip direction
      V(:,:,r)=[cosd(90-dipdireigv(r))*cosd(dipeigv(r)) sind(90-dipdireigv(r))*cosd(dipeigv(r)) -sind(dipeigv(r))];
      %Unit vector O, which origin is at the origo, which is perpendicular to 
      %the dip direction and is horizontal
      O(:,:,r)=[-sind(90-dipdireigv(r)) cosd(90-dipdireigv(r)) 0];
   end

   %Normal vector of a discontinuity set is calculated with the cross product of
   %vectors V and O.
   %Length of the unit vectors is one.
   for k=1:max(idx3)
      angle_between(k)=acosd(dot(V(:,:,k),O(:,:,k))); % Vectors are perpendicular, so the angle between them is 90 degrees.
   end
   %Normal vector N of the set
   N=zeros(1,3,max(idx3));
   for k=1:max(idx3)
      N(:,:,k)=cross(V(:,:,k),O(:,:,k))/sind(angle_between(k));
   end
   %Angle between the scanline and normal of the discontinuity set is calculated with dot product -> correction factor
   %Scanline as a vector
   xA=cosd(sltrend)*cosd(slplunge);
   yA=sind(sltrend)*cosd(slplunge);
   zA=-sind(slplunge);
   A=[xA yA zA];
   %N is the normal vector of the set
   for k=1:max(idx3)
      ypsilon(k)=acosd(dot(A,N(:,:,k)));
      correctionfactor(k)=abs(cosd(ypsilon(k)));
      if correctionfactor (k)<0.2     % (Hofrichter and Winkler, 2006)
         correctionfactor(k)=0.2;
      end
   end
   %Spacing for block volumes out of drilling cores 
   for f=2:numel(data(:,1)) 
      spacingscanline(f-1)=intersection(f)-intersection(f-1);
      Spacingsforscanline=spacingscanline;
   end

   %Intersections of the discontinuities included in the set, their spacing and mean spacing
   clear intersections_set, clear spacing_set, clear f, clear n1
   for n1=1:max(idx3)
      intersections_set=intersection(find(idx3==n1));
      RESULTS.(genvarname(['Set', num2str(n1)])).Intersections=intersections_set';
      v = genvarname(['intersections_set', num2str(n1)]);
      %eval([v ' = intersections_set']);
      if (sum(idx3==n1)>1)
         for f=2:sum(idx3==n1) % Number of spacings is one less than number of discontinuities
            spacing_set(f-1)=intersections_set(f)-intersections_set(f-1);
            %RESULTS.(genvarname(['Set', num2str(n1)])).DiscontinuitySpacings=spacing_set;
            corrected_spacing_set(f-1)=spacing_set(f-1)*correctionfactor(n1);
            RESULTS.(genvarname(['Set', num2str(n1)])).Normal_spacings=corrected_spacing_set;
         end
      else
      spacing_set=[NaN]; % When a cluster contains only one discontinuity
      end
      v = genvarname(['spacing_set', num2str(n1)]);
      %eval([v ' = spacing_set']);
      mean_spacing_set=mean(spacing_set);
      median_spacing_set=median(spacing_set);
      mode_spacing_set=mode(spacing_set);
      std_spacing_set=std(spacing_set);
      std_spacing_in_sets(n1)=std_spacing_set;
      RESULTS.(genvarname(['Set', num2str(n1)])).STD_normal_spacing=std_spacing_in_sets(n1);
      v = genvarname(['mean_spacing_set', num2str(n1)]);
      %eval([v ' = mean_spacing_set']);
      corrected_mean_spacing_set=mean_spacing_set*correctionfactor(n1);
      corrected_median_spacing_set=median_spacing_set*correctionfactor(n1);
      corrected_mode_spacing_set=mode_spacing_set*correctionfactor(n1);
      v = genvarname(['corrected_mean_spacing_set', num2str(n1)]);
      %eval([v ' = corrected_mean_spacing_set']);
      corrected_meanspacings_in_sets(n1)=corrected_mean_spacing_set;
      corrected_medianspacings_in_sets(n1)=corrected_median_spacing_set;
      corrected_modespacings_in_sets(n1)=corrected_mode_spacing_set;
      clear intersections_set, clear spacing_set, clear corrected_spacing_set
   end
   for n1=1:max(idx3)
      RESULTS.(genvarname(['Set', num2str(n1)])).Mean_normal_spacing=corrected_meanspacings_in_sets(n1);     
      RESULTS.(genvarname(['Set', num2str(n1)])).Median_normal_spacing=corrected_medianspacings_in_sets(n1); 
      RESULTS.(genvarname(['Set', num2str(n1)])).Mode_normal_spacing=corrected_modespacings_in_sets(n1);
   end
   
RESULTS.Rock_characteristics.Input_parameters.Spacingsforscanline=Spacingsforscanline;
   

   
%Calculates volumetric joint count 
jv=sum(1./corrected_meanspacings_in_sets);
RESULTS.Rock_characteristics.Jv=jv;
set (handles.text_jvvalueblock,'string',RESULTS.Rock_characteristics.Jv);

%CALCULATION OF RQDs (for 1m for drilling core comparision)
intersectionb=data(:,1);
intersection=[0 intersectionb' sllength];
lengthsl=floor(sllength);
RQDsscanline=zeros(1,lengthsl);
for i=1:floor(sllength)
   h=intersection(intersection<i & intersection>(i-1));
   r=(i-1);
   fracturelocations=[r h i]-i+1;
   fracturespac=diff(fracturelocations);
   subfracspac=fracturespac(fracturespac>0.09);
   rqd=sum(subfracspac)*100;
   RQDsscanline(1,i)=rqd;
end
meanrqd=mean(RQDsscanline);

%CALCULATION OF RQDs (for full scanline length for RMC calculations)
fullscanintersections=diff(intersection);
subfullscanintersections=fullscanintersections(fullscanintersections>0.09);
fullscanrqd=sum(subfullscanintersections)*100/sllength;
RESULTS.Rock_characteristics.RQDs=RQDsscanline;


RESULTS.Rock_characteristics.Mean_RQD=mean(RQDsscanline);
RESULTS.Rock_characteristics.Min_RQD=min(RQDsscanline);
RESULTS.Rock_characteristics.Max_RQD=max(RQDsscanline);





% --- Executes on button press in pushbutton_savestatistics.
function pushbutton_savestatistics_Callback(hObject, eventdata, handles)
global RESULTS
% hObject    handle to pushbutton_savestatistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uisave('RESULTS','Results');


% --- Executes on button press in pushbutton_vol.
function pushbutton_vol_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS
set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','on');
set(handles.uipanel_merge,'Visible','off');
set(handles.panel_rmc,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');
set(handles.text_jvvalueblock,'string',RESULTS.Rock_characteristics.Jv);


% --- Executes on button press in pushbutton_drawgraphs.
function pushbutton_drawgraphs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_drawgraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RESULTS idx3


terminationsscan=RESULTS.Joint_characteristics.Terminations;
tracelengthsscan=RESULTS.Joint_characteristics.Tracelengths;
jointalterationscsan=RESULTS.Joint_characteristics.Joint_Alteration_Factor;
jointroughnessscan=RESULTS.Joint_characteristics.Joint_Roughness_Factor;
frictionanglescan=RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles;
jointsizescan=RESULTS.Joint_characteristics.Joint_Size_Factor;


if get(handles.check_box_scan_term,'Value') == 1

   figure('name', 'Terminations of scanline discontinuities', 'NumberTitle', 'off')
   title('Terminations of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(terminationsscan)
   xlabel('Scanline');
   ylabel('Joint terminations');
end


if get(handles.check_box_scan_spacing,'Value') == 1
   clear n1
   spacings_all=[];
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      spacings_all=horzcat(spacings_all,RESULTS.(setname).Normal_spacings);
   end
   figure('name', 'Spacings of scanline discontinuities', 'NumberTitle', 'off')
   title('Spacings of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(spacings_all)
   xlabel('Scanline');
   ylabel('Discontinuity spacing [m]');
end

if get(handles.check_box__scan_tracelength,'Value') == 1
   figure('name', 'Trace lengths of scanline discontinuities', 'NumberTitle', 'off')
   title('Trace lengths of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(tracelengthsscan)
   xlabel('Scanline');
   ylabel('Trace length [m]');
end

if get(handles.check_box_scan_phi,'Value') == 1

   figure('name', 'Peak friction angles of scanline discontinuities', 'NumberTitle', 'off')
   title('Peak friction angles of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(frictionanglescan)
   xlabel('Scanline');
   ylabel('Peak friction angle [°]');
end

if get(handles.check_box_scan_jointsizefactor,'Value') == 1
   figure('name', 'Joint size factor of scanline discontinuities', 'NumberTitle', 'off')
   title('Joint size factor of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(jointsizescan)
   xlabel('Scanline');
   ylabel('Joint size factor');
end

if get(handles.check_box_scan_roughness,'Value') == 1
   figure('name', 'Joint roughness factor of scanline discontinuities', 'NumberTitle', 'off')
   title('Joint roughness factor of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(jointroughnessscan)
   xlabel('Scanline');
   ylabel('Joint roughness factor');
end

if get(handles.check_box_scan_alteration,'Value') == 1
   figure('name', 'Joint alteration factor of scanline discontinuities', 'NumberTitle', 'off')
   title('Joint alteration factor of scanline discontinuities', 'fontsize', 14)
   axis off
   boxplot(jointalterationscsan)
   xlabel('Scanline');
   ylabel('Joint alteration factor');
end

if get(handles.check_box_sets_term,'Value') == 1
   figure('name', 'Terminations of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Terminations of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Terminations)
      xlabel(setname);
      ylabel('Joint termination');
   end    
end

if get(handles.check_box_sets_spacing,'Value') == 1
   figure('name', 'Spacings of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Spacings of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Normal_spacings)
      xlabel(setname);
      ylabel('Discontinuity spacing [m]');
   end    
end

if get(handles.check_box_sets_tracelengths,'Value') == 1
   figure('name', 'Trace lengths of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Trace lengths of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Tracelengths)
      xlabel(setname);
      ylabel('Discontinuity trace length [m]');
   end    
end

if get(handles.check_box_sets_phi,'Value') == 1
   figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Peak friction angle of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Peak_Friction_Angles)
      xlabel(setname);
      ylabel('Peak friction angle [°]');
   end   
end

if get(handles.check_box_sets_jointsizefactor,'Value') == 1
   figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Joint size factor of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Joint_Size_Factor)
      xlabel(setname);
      ylabel('Joint size factor');
   end    
end

if get(handles.check_box_sets_roughness,'Value') == 1
   figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Joint_Roughness_Factors)
      xlabel(setname);
      ylabel('Joint roughness factor');
   end   
end


if get(handles.check_box_sets_alteration,'Value') == 1
   figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
   sm=max(idx3);
   subplot(10,sm,[1:sm])
   title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
   axis off
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      picframes(n1,:)=[n2*sm+n1];
      subplot(10,sm,picframes(n1,:))
      boxplot(RESULTS.(setname).Joint_Alteration_Factors)
      xlabel(setname);
      ylabel('Joint alteration factor');
   end   
end

if get(handles.check_hist_scan_spacing,'Value') == 1
   spacings_all=[];
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      spacings_all=horzcat(spacings_all,RESULTS.(setname).Normal_spacings);
   end
   maxspac=max(spacings_all);
   minspac=min(spacings_all);
   bins = linspace(minspac,maxspac,12);       
   fcspac= histc(spacings_all,bins);                
   fcspacn = fcspac/sum(fcspac);   
   figure('name', 'Spacings of scanline discontinuities', 'NumberTitle', 'off') 
   bar(bins',[fcspacn'])  
   title('Spacings of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end

if get(handles.check_hist_scan_tracelength,'Value') == 1
   fulltrace=tracelengthsscan';
   mintrac=min(fulltrace);
   maxtrac=max(fulltrace);
   binstrac = linspace(mintrac,maxtrac,12);   
   fctrac= histc(fulltrace,binstrac);
   fctracn = fctrac/sum(fctrac); 
   figure('name', 'Trace lengths of scanline discontinuities', 'NumberTitle', 'off')
   bar(binstrac',[fctracn'])   
   title('Trace lengths of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end


if get(handles.check_hist_scan_phi,'Value') == 1
   histfrictscan=frictionanglescan';
   maxphi=max(histfrictscan);
   minphi=min(histfrictscan);
   binsphi = linspace(minphi,maxphi,12);   
   fcphi= histc(histfrictscan,binsphi);
   fcphin = fcphi/sum(fcphi); 
   figure('name', 'Peak friction angles of scanline discontinuities', 'NumberTitle', 'off')
   bar(binsphi',[fcphin']);
   title('Peak friction angles of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end

if get(handles.check_hist_scan_jointsizefactor,'Value') == 1
   histsizescan=jointsizescan';
   maxwav=max(histsizescan);
   minwav=min(histsizescan);
   binswav = linspace(minwav,maxwav,12);   
   fcwav= histc(histsizescan,binswav);
   fcwavn = fcwav/sum(fcwav); 
   figure('name', 'Joint size factor of scanline discontinuities', 'NumberTitle', 'off')
   bar(binswav',[fcwavn']);
   title('Joint size factor of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end


if get(handles.ckeck_hist_scan_roughness,'Value') == 1
   histroughscan=jointroughnessscan';
   maxrough=max(histroughscan);
   minrough=min(histroughscan);
   binsrough = linspace(minrough,maxrough,12);   
   fcrough= histc(histroughscan,binsrough);
   fcroughn = fcrough/sum(fcrough); 
   figure('name', 'Joint roughness factor of scanline discontinuities', 'NumberTitle', 'off')
   bar(binsrough',[fcroughn']);
   title('Joint roughness factor of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end


if get(handles.check_hist_scan_alteration,'Value') == 1
   histaltscan=jointalterationscsan';
   maxweath=max(histaltscan);
   minweath=min(histaltscan);
   binsweath = linspace(minweath,maxweath,12);   
   fcweath= histc(histaltscan,binsweath);
   fcweathn = fcweath/sum(fcweath); 
   figure('name', 'Joint alteration factor of scanline discontinuities', 'NumberTitle', 'off')
   bar(binsweath',[fcweathn']);
   title('Joint alteration of scanline discontinuities', 'fontsize', 14)
   ylabel('Frequency []');
end


if get(handles.check_hist_sets_spacing,'Value') == 1
   spacings_all=[];
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      spacings_all=horzcat(spacings_all,RESULTS.(setname).Normal_spacings);
   end
   maxspac=max(spacings_all);
   minspac=min(spacings_all);
   sm=max(idx3);
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      bins = linspace(minspac,maxspac,12);      
      fc.(setname)= histc(RESULTS.(setname).Normal_spacings,bins);               
      fcn.(setname) = fc.(setname)/sum(fc.(setname));                  
   end
   if max(idx3)==1
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1']) 
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==2
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2']) 
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==3
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2' fcn.Set3']) 
      title('Spacings of discontinuity sets', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==4
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2' fcn.Set3' fcn.Set4'])
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==5
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2' fcn.Set3' fcn.Set4' fcn.Set5']) 
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==6
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2' fcn.Set3' fcn.Set4' fcn.Set5' fcn.Set6'])
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==7
      figure('name', 'Spacings of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[fcn.Set1' fcn.Set2' fcn.Set3' fcn.Set4' fcn.Set5' fcn.Set6' fcn.Set7']) 
      title('Spacings of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
end


if get(handles.check_hist_sets_tracelength,'Value') == 1
   tracelengths_all=[];
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      tracelengths_all=horzcat(tracelengths_all,RESULTS.(setname).Tracelengths);
   end
   maxtrac=max(tracelengths_all);
   mintrac=min(tracelengths_all);
   sm=max(idx3);
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      bins = linspace(mintrac,maxtrac,12);      
      ftrac.(setname)= histc(RESULTS.(setname).Tracelengths,bins);               
      ftracn.(setname) = ftrac.(setname)/sum(ftrac.(setname));                  
   end
   if max(idx3)==1
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1']) 
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==2
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2']) 
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==3
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2' ftracn.Set3']) 
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==4
      figure('name', 'Trace lenghts of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2' ftracn.Set3' ftracn.Set4'])
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
   end 
   if max(idx3)==5
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2' ftracn.Set3' ftracn.Set4' ftracn.Set5']) 
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==6
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2' ftracn.Set3' ftracn.Set4' ftracn.Set5' ftracn.Set6'])
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==7
      figure('name', 'Trace lengths of discontinuity sets [m]', 'NumberTitle', 'off')
      bar(bins',[ftracn.Set1' ftracn.Set2' ftracn.Set3' ftracn.Set4' ftracn.Set5' ftracn.Set6' ftracn.Set7']) 
      title('Trace lengths of discontinuity sets [m]', 'fontsize', 14)
      ylabel('Frequency []');
   end 
end



if get(handles.check_hist_sets_phi,'Value') == 1
   phi_all=[];
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      phi_all=horzcat(phi_all,RESULTS.(setname).Peak_Friction_Angles);
   end
   maxphi=max(phi_all);
   minphi=min(phi_all);
   sm=max(idx3);
   clear n1
   for n1=1:max(idx3)
      setname=['Set', num2str(n1)];
      sm=max(idx3);
      n2=1:9;
      bins = linspace(minphi,maxphi,12);      
      fphi.(setname)= histc(RESULTS.(setname).Peak_Friction_Angles,bins);               
      fphin.(setname) = fphi.(setname)/sum(fphi.(setname));                  
   end
   if max(idx3)==1
      figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
      bar(bins',[fphin.Set1']) 
      title('Peak friction angle of discontinuity sets', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==2
      figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
      bar(bins',[fphin.Set1' fphin.Set2']) 
      title('Peak friction angle of discontinuity sets', 'fontsize', 14)
      ylabel('Frequency []');
   end 
   if max(idx3)==3
    figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fphin.Set1' fphin.Set2' fphin.Set3']) 
    title('Peak friction angle of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==4
    figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fphin.Set1' fphin.Set2' fphin.Set3' fphin.Set4'])
    title('Peak friction angle of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==5
    figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fphin.Set1' fphin.Set2' fphin.Set3' fphin.Set4' fphin.Set5']) 
    title('Peak friction angle of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==6
    figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fphin.Set1' fphin.Set2' fphin.Set3' fphin.Set4' fphin.Set5' fphin.Set6'])
    title('Peak friction angle of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==7
    figure('name', 'Peak friction angle of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fphin.Set1' fphin.Set2' fphin.Set3' fphin.Set4' fphin.Set5' fphin.Set6' fphin.Set7']) 
    title('Peak friction angle of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
end


if get(handles.check_hist_sets_jointsizefactor,'Value') == 1
    size_all=[];
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    size_all=horzcat(size_all,RESULTS.(setname).Joint_Size_Factor);
    end
    maxsize=max(size_all);
    minsize=min(size_all);
    sm=max(idx3);
    clear n1
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    sm=max(idx3);
    n2=1:9;
    bins = linspace(minsize,maxsize,12);      
    fsize.(setname)= histc(RESULTS.(setname).Joint_Size_Factor,bins);     
    fsizen.(setname) = fsize.(setname)/sum(fsize.(setname));                  
    end
    if max(idx3)==1
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1']) 
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==2
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2']) 
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==3
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2' fsizen.Set3']) 
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==4
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2' fsizen.Set3' fsizen.Set4'])
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==5
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2' fsizen.Set3' fsizen.Set4' fsizen.Set5']) 
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==6
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2' fsizen.Set3' fsizen.Set4' fsizen.Set5' fsizen.Set6'])
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==7
    figure('name', 'Joint size factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fsizen.Set1' fsizen.Set2' fsizen.Set3' fsizen.Set4' fsizen.Set5' fsizen.Set6' fsizen.Set7']) 
    title('Joint size factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
end


if get(handles.check_hist_sets_roughness,'Value') == 1
    roughness_all=[];
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    roughness_all=horzcat(roughness_all,RESULTS.(setname).Joint_Roughness_Factors);
    end
    maxroughness=max(roughness_all);
    minroughness=min(roughness_all);
    sm=max(idx3);
    clear n1
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    sm=max(idx3);
    n2=1:9;
    bins = linspace(minroughness,maxroughness,12);      
    froughness.(setname)= histc(RESULTS.(setname).Joint_Roughness_Factors,bins);     
    froughnessn.(setname) = froughness.(setname)/sum(froughness.(setname));                  
    end
    if max(idx3)==1
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1']) 
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==2
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2']) 
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==3
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2' froughnessn.Set3']) 
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==4
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2' froughnessn.Set3' froughnessn.Set4'])
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==5
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2' froughnessn.Set3' froughnessn.Set4' froughnessn.Set5']) 
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==6
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2' froughnessn.Set3' froughnessn.Set4' froughnessn.Set5' froughnessn.Set6'])
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==7
    figure('name', 'Joint roughness factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[froughnessn.Set1' froughnessn.Set2' froughnessn.Set3' froughnessn.Set4' froughnessn.Set5' froughnessn.Set6' froughnessn.Set7']) 
    title('Joint roughness factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end
end


if get(handles.check_hist_sets_alteration,'Value') == 1
    weathering_all=[];
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    weathering_all=horzcat(weathering_all,RESULTS.(setname).Joint_Alteration_Factors);
    end
    maxweathering=max(weathering_all);
    minweathering=min(weathering_all);
    sm=max(idx3);
    clear n1
    for n1=1:max(idx3)
    setname=['Set', num2str(n1)];
    sm=max(idx3);
    n2=1:9;
    bins = linspace(minweathering,maxweathering,12);      
    fweathering.(setname)= histc(RESULTS.(setname).Joint_Alteration_Factors,bins);     
    fweatheringn.(setname) = fweathering.(setname)/sum(fweathering.(setname));                  
    end
    if max(idx3)==1
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1']) 
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==2
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2']) 
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==3
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2' fweatheringn.Set3']) 
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==4
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2' fweatheringn.Set3' fweatheringn.Set4'])
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==5
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2' fweatheringn.Set3' fweatheringn.Set4' fweatheringn.Set5']) 
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==6
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2' fweatheringn.Set3' fweatheringn.Set4' fweatheringn.Set5' fweatheringn.Set6'])
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end 
    if max(idx3)==7
    figure('name', 'Joint alteration factor of discontinuity sets', 'NumberTitle', 'off')
    bar(bins',[fweatheringn.Set1' fweatheringn.Set2' fweatheringn.Set3' fweatheringn.Set4' fweatheringn.Set5' fweatheringn.Set6' fweatheringn.Set7']) 
    title('Joint alteration factor of discontinuity sets', 'fontsize', 14)
    ylabel('Frequency []');
    end
end



if get(handles.check_rqd_hist,'Value') == 1
    rqds_all=RESULTS.Rock_characteristics.RQDs;
    maxrqd=max(rqds_all);
    minrqd=min(rqds_all);
    bins = linspace(minrqd,maxrqd,6);       
fcrqd= histc(rqds_all,bins);                
fcrqdn = fcrqd/sum(fcrqd);   
 figure('name', 'RQDs of scanline', 'NumberTitle', 'off') 
bar(bins',[fcrqdn'])  
 title('RQDs of scanline', 'fontsize', 14)
 ylabel('Frequency []');
end


if get(handles.check_rqd_distance,'Value') == 1
      rqds_all=RESULTS.Rock_characteristics.RQDs;
   numberrqds=length(rqds_all)-0.5;
   depth=0.5:1:numberrqds;
    rqddistfig=figure('name', 'RQDs vs. distance', 'NumberTitle', 'off');
   set(rqddistfig, 'Position', [100 100 250 600])
     title('RQDs vs. distance', 'fontsize', 14)
    axis off
    barh(depth,rqds_all);
set(gca,'YDir','reverse');
    xlabel('RQD [%]');
    ylabel('Distance [m]');
end


% --- Executes on selection change in popup_sets_c.
function popup_sets_c_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sets_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sets_c contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sets_c
global  RESULTS 

switch get(handles.popup_sets_c,'Value')  
    case 2
         spacesetc=RESULTS.Set1.Mean_normal_spacing;
         dipdirc=RESULTS.Set1.Dipdir;
         dipc=RESULTS.Set1.Dip;
         dipdirsc=RESULTS.Set1.Dipdirs;
         dipsc=RESULTS.Set1.Dips;
         spacingsofsetc=RESULTS.Set1.Normal_spacings;
    case 3
         spacesetc=RESULTS.Set2.Mean_normal_spacing;
         dipdirc=RESULTS.Set2.Dipdir;
         dipc=RESULTS.Set2.Dip;
         dipdirsc=RESULTS.Set2.Dipdirs;
         dipsc=RESULTS.Set2.Dips;
         spacingsofsetc=RESULTS.Set2.Normal_spacings;
    case 4
         spacesetc=RESULTS.Set3.Mean_normal_spacing;
         dipdirc=RESULTS.Set3.Dipdir;
         dipc=RESULTS.Set3.Dip;
         dipdirsc=RESULTS.Set3.Dipdirs;
         dipsc=RESULTS.Set3.Dips;
         spacingsofsetc=RESULTS.Set3.Normal_spacings;
    case 5
         spacesetc=RESULTS.Set4.Mean_normal_spacing;
         dipdirc=RESULTS.Set4.Dipdir;
         dipc=RESULTS.Set4.Dip;
         dipdirsc=RESULTS.Set4.Dipdirs;
         dipsc=RESULTS.Set4.Dips;
         spacingsofsetc=RESULTS.Set4.Normal_spacings;
    case 6
         spacesetc=RESULTS.Set5.Mean_normal_spacing;
         dipdirc=RESULTS.Set5.Dipdir;
         dipc=RESULTS.Set5.Dip;
         dipdirsc=RESULTS.Set5.Dipdirs;
         dipsc=RESULTS.Set5.Dips;
         spacingsofsetc=RESULTS.Set5.Normal_spacings;
    case 7
         spacesetc=RESULTS.Set6.Mean_normal_spacing;
         dipdirc=RESULTS.Set6.Dipdir;
         dipc=RESULTS.Set6.Dip;
         dipdirsc=RESULTS.Set6.Dipdirs;
         dipsc=RESULTS.Set6.Dips;
         spacingsofsetc=RESULTS.Set6.Normal_spacings;
    case 8
         spacesetc=RESULTS.Set7.Mean_normal_spacing;
         dipdirc=RESULTS.Set7.Dipdir;
         dipc=RESULTS.Set7.Dip;
         dipdirsc=RESULTS.Set7.Dipdirs;
         dipsc=RESULTS.Set7.Dips;
         spacingsofsetc=RESULTS.Set7.Normal_spacings;
end
set(handles.edit_spacing_setc,'Visible','on');
set(handles.text_spacings_blockvol,'Visible','on');


set(handles.edit_spacing_setc, 'String', spacesetc);
RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetC=spacingsofsetc;
RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetC=dipdirc;
RESULTS.Rock_characteristics.Input_parameters.Dip_SetC=dipc;
RESULTS.Rock_characteristics.Input_parameters.Mean_normal_spacing_SetC=spacesetc;
RESULTS.Rock_characteristics.Input_parameters.Dips_SetC=dipsc;
RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetC=dipdirsc;




% --- Executes on selection change in popup_sets_b.
function popup_sets_b_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sets_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sets_b contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sets_b
global  RESULTS 
switch get(handles.popup_sets_b,'Value')  
    case 2
         spacesetb=RESULTS.Set1.Mean_normal_spacing;
         dipdirb=RESULTS.Set1.Dipdir;
         dipb=RESULTS.Set1.Dip;
         dipdirsb=RESULTS.Set1.Dipdirs;
         dipsb=RESULTS.Set1.Dips;
         spacingsofsetb=RESULTS.Set1.Normal_spacings;
    case 3
         spacesetb=RESULTS.Set2.Mean_normal_spacing;
         dipdirb=RESULTS.Set2.Dipdir;
         dipb=RESULTS.Set2.Dip;
         dipdirsb=RESULTS.Set2.Dipdirs;
         dipsb=RESULTS.Set2.Dips;
         spacingsofsetb=RESULTS.Set2.Normal_spacings;
    case 4
         spacesetb=RESULTS.Set3.Mean_normal_spacing;
         dipdirb=RESULTS.Set3.Dipdir;
         dipb=RESULTS.Set3.Dip;
         dipdirsb=RESULTS.Set3.Dipdirs;
         dipsb=RESULTS.Set3.Dips;
         spacingsofsetb=RESULTS.Set3.Normal_spacings;
    case 5
         spacesetb=RESULTS.Set4.Mean_normal_spacing;
         dipdirb=RESULTS.Set4.Dipdir;
         dipb=RESULTS.Set4.Dip;
         dipdirsb=RESULTS.Set4.Dipdirs;
         dipsb=RESULTS.Set4.Dips;
         spacingsofsetb=RESULTS.Set4.Normal_spacings;
    case 6
         spacesetb=RESULTS.Set5.Mean_normal_spacing;
         dipdirb=RESULTS.Set5.Dipdir;
         dipb=RESULTS.Set5.Dip;
         dipdirsb=RESULTS.Set5.Dipdirs;
         dipsb=RESULTS.Set5.Dips;
         spacingsofsetb=RESULTS.Set5.Normal_spacings;
    case 7
         spacesetb=RESULTS.Set6.Mean_normal_spacing;
         dipdirb=RESULTS.Set6.Dipdir;
         dipb=RESULTS.Set6.Dip;
         dipdirsb=RESULTS.Set6.Dipdirs;
         dipsb=RESULTS.Set6.Dips;
         spacingsofsetb=RESULTS.Set6.Normal_spacings;
    case 8
         spacesetb=RESULTS.Set7.Mean_normal_spacing;
         dipdirb=RESULTS.Set7.Dipdir;
         dipb=RESULTS.Set7.Dip;
         dipdirsb=RESULTS.Set7.Dipdirs;
         dipsb=RESULTS.Set7.Dips;
         spacingsofsetb=RESULTS.Set7.Normal_spacings;
end

set(handles.edit_spacing_setb,'Visible','on');
set(handles.text_spacings_blockvol,'Visible','on');


set(handles.edit_spacing_setb, 'String', spacesetb);
RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetB=spacingsofsetb;
RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetB=dipdirb;
RESULTS.Rock_characteristics.Input_parameters.Dip_SetB=dipb;
RESULTS.Rock_characteristics.Input_parameters.Mean_normal_spacing_SetB=spacesetb;
RESULTS.Rock_characteristics.Input_parameters.Dips_SetB=dipsb;
RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetB=dipdirsb;






% --- Executes on selection change in popup_sets_a.
function popup_sets_a_Callback(hObject, eventdata, handles)
% hObject    handle to popup_sets_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_sets_a contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_sets_a
global RESULTS 
switch get(handles.popup_sets_a,'Value')  
    case 2
         spaceseta=RESULTS.Set1.Mean_normal_spacing;
         dipdira=RESULTS.Set1.Dipdir;
         dipa=RESULTS.Set1.Dip;
         dipdirsa=RESULTS.Set1.Dipdirs;
         dipsa=RESULTS.Set1.Dips;
         spacingsofseta=RESULTS.Set1.Normal_spacings;
    case 3
         spaceseta=RESULTS.Set2.Mean_normal_spacing;
         dipdira=RESULTS.Set2.Dipdir;
         dipa=RESULTS.Set2.Dip;
         dipdirsa=RESULTS.Set2.Dipdirs;
         dipsa=RESULTS.Set2.Dips;
         spacingsofseta=RESULTS.Set2.Normal_spacings;
    case 4
         spaceseta=RESULTS.Set3.Mean_normal_spacing;
         dipdira=RESULTS.Set3.Dipdir;
         dipa=RESULTS.Set3.Dip;
         dipdirsa=RESULTS.Set3.Dipdirs;
         dipsa=RESULTS.Set3.Dips;
         spacingsofseta=RESULTS.Set3.Normal_spacings;
    case 5
         spaceseta=RESULTS.Set4.Mean_normal_spacing;
         dipdira=RESULTS.Set4.Dipdir;
         dipa=RESULTS.Set4.Dip;
         dipdirsa=RESULTS.Set4.Dipdirs;
         dipsa=RESULTS.Set4.Dips;
         spacingsofseta=RESULTS.Set4.Normal_spacings;
    case 6
         spaceseta=RESULTS.Set5.Mean_normal_spacing;
         dipdira=RESULTS.Set5.Dipdir;
         dipa=RESULTS.Set5.Dip;
         dipdirsa=RESULTS.Set5.Dipdirs;
         dipsa=RESULTS.Set5.Dips;
         spacingsofseta=RESULTS.Set5.Normal_spacings;
    case 7
         spaceseta=RESULTS.Set6.Mean_normal_spacing;
         dipdira=RESULTS.Set6.Dipdir;
         dipa=RESULTS.Set6.Dip;
         dipdirsa=RESULTS.Set6.Dipdirs;
         dipsa=RESULTS.Set6.Dips;
         spacingsofseta=RESULTS.Set6.Normal_spacings;
    case 8
         spaceseta=RESULTS.Set7.Mean_normal_spacing;
         dipdira=RESULTS.Set7.Dipdir;
         dipa=RESULTS.Set7.Dip;
         dipdirsa=RESULTS.Set7.Dipdirs;
         dipsa=RESULTS.Set7.Dips;
         spacingsofseta=RESULTS.Set7.Normal_spacings;
end
set(handles.edit_spacing_seta,'Visible','on');
set(handles.text_spacings_blockvol,'Visible','on');


set(handles.edit_spacing_seta, 'String', spaceseta);

RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetA=spacingsofseta;
RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetA=dipdira;
RESULTS.Rock_characteristics.Input_parameters.Dip_SetA=dipa;
RESULTS.Rock_characteristics.Input_parameters.Mean_normal_spacing_SetA=spaceseta;
RESULTS.Rock_characteristics.Input_parameters.Dips_SetA=dipsa;
RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetA=dipdirsa;



% --- Executes on button press in pushbutton_saveblockvolume.
function pushbutton_saveblockvolume_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveblockvolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global vnegexp vfractal vuniform vexpuni vlognormal RESULTS
Blockvol=figure('name', 'Block Volume', 'NumberTitle', 'off');
passing=[10;20;30;40;50;60;70;80;90;100];
vfractal=RESULTS.Block_volume_calculations.ISBDs.Fractal;
vlognormal=RESULTS.Block_volume_calculations.ISBDs.Lognormal;
vexpuni=RESULTS.Block_volume_calculations.ISBDs.TwoNegExp_OneUniform;
vuniform=RESULTS.Block_volume_calculations.ISBDs.Uniform;
vnegexp=RESULTS.Block_volume_calculations.ISBDs.Negativeexponential;


if get(handles.radio_3negexp,'Value') == 1
    semilogx(vnegexp,passing,'LineWidth',3,'Color','red')
    hold on
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    %hold on
    legend('3 x Neg. exponential (persistent)','Location','northwest');
end
if get(handles.radio_3uniform,'Value') == 1
    semilogx(vuniform,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    legend('3 x Uniform (persistent)','Location','northwest');
end
if get(handles.radio_2negexp1uni,'Value') == 1
    semilogx(vexpuni,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    legend('2 x Neg. Exponential, 1 x Uniform (persistent)','Location','northwest');
end
if get(handles.radio_3lognormal,'Value') == 1
    semilogx(vlognormal,passing,'LineWidth',3,'Color','red') 
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    legend('3 x Lognormal (persistent)','Location','northwest');
end
if get(handles.radio_3fractal,'Value') == 1
    semilogx(vfractal,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    legend('3 x Fractal (persistent)','Location','northwest');
end
hold off


% --- Executes on button press in pushbutton_create_bsd.
function pushbutton_create_bsd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_create_bsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global    Spacingsforscanline RESULTS sort_ibsd_palm  yaxis_palm isbdvdrill perc_isbd_vdrill vnegexp vfractal vuniform vexpuni vlognormal angleab angleac anglebc spacingsofseta spacingsofsetb spacingsofsetc jv beta s

spaceseta =str2double(get (handles.edit_spacing_seta,'String'));
spacesetb =str2double(get (handles.edit_spacing_setb,'String'));
spacesetc =str2double(get (handles.edit_spacing_setc,'String'));

spacingsofsetc=RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetC;
dipdirc=RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetC;
dipc=RESULTS.Rock_characteristics.Input_parameters.Dip_SetC;
spacingsofsetb=RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetB;
dipdirb=RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetB;
dipb=RESULTS.Rock_characteristics.Input_parameters.Dip_SetB;
spacingsofseta=RESULTS.Rock_characteristics.Input_parameters.Normal_spacings_SetA;
dipdira=RESULTS.Rock_characteristics.Input_parameters.Dipdir_SetA;
dipa=RESULTS.Rock_characteristics.Input_parameters.Dip_SetA;
dipdirsa=RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetA;
dipsa=RESULTS.Rock_characteristics.Input_parameters.Dips_SetA;
dipdirsb=RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetB;
dipsb=RESULTS.Rock_characteristics.Input_parameters.Dips_SetB;
dipdirsc=RESULTS.Rock_characteristics.Input_parameters.Dipdirs_SetC;
dipsc=RESULTS.Rock_characteristics.Input_parameters.Dips_SetC;

%For Mean Block size calculation
if dipdira>180
    dipdiran=dipdira-180;
else 
    dipdiran=dipdira+180;
end

if dipdirb>180
    dipdirbn=dipdirb-180;
else 
    dipdirbn=dipdirb+180;
end

if dipdirc>180
    dipdircn=dipdirc-180;
else 
    dipdircn=dipdirc+180;
end

dipan=90-dipa;
dipbn=90-dipb;
dipcn=90-dipc;

angleab=90-(acosd(abs((cosd(dipdiran-dipdirbn)*cosd(dipan)*cosd(dipbn))+(sind(dipan)*sind(dipbn)))));
angleac=90-(acosd(abs((cosd(dipdirbn-dipdircn)*cosd(dipbn)*cosd(dipcn))+(sind(dipbn)*sind(dipcn)))));
anglebc=90-(acosd(abs((cosd(dipdiran-dipdircn)*cosd(dipan)*cosd(dipcn))+(sind(dipan)*sind(dipcn)))));


%For  ISBD calculation
randpicka=length(dipdirsa);
for i=1:5000
    randorient=randi([1 randpicka],1,1);
    rorientations_a(1,i)=dipdirsa(1,randorient);
    rorientations_a(2,i)=dipsa(1,randorient);
end

randpickb=length(dipdirsb);
for i=1:5000
    randorient=randi([1 randpickb],1,1);
    rorientations_b(1,i)=dipdirsb(1,randorient);
    rorientations_b(2,i)=dipsb(1,randorient);
end

randpickc=length(dipdirsc);
for i=1:5000
    randorient=randi([1 randpickc],1,1);
    rorientations_c(1,i)=dipdirsc(1,randorient);
    rorientations_c(2,i)=dipsc(1,randorient);
end



for i=1:length(rorientations_a)
if rorientations_a(1,i)>180
    randdipdirsan(1,i)=rorientations_a(1,i)-180;
else 
    randdipdirsan(1,i)=rorientations_a(1,i)+180;
end
if rorientations_b(1,i)>180
    randdipdirsbn(1,i)=rorientations_b(1,i)-180;
else 
    randdipdirsbn(1,i)=rorientations_b(1,i)+180;
end
if rorientations_c(1,i)>180
    randdipdirscn(1,i)=rorientations_c(1,i)-180;
else 
    randdipdirscn(1,i)=rorientations_c(1,i)+180;
end

randdipsan(1,i)=90-rorientations_a(2,i);
randdipsbn(1,i)=90-rorientations_b(2,i);
randdipscn(1,i)=90-rorientations_c(2,i);
end


randangleab=90-(acosd(abs((cosd(randdipdirsan-randdipdirsbn).*cosd(randdipsan).*cosd(randdipsbn))+(sind(randdipsan).*sind(randdipsbn)))));
randangleac=90-(acosd(abs((cosd(randdipdirsbn-randdipdirscn).*cosd(randdipsbn).*cosd(randdipscn))+(sind(randdipsbn).*sind(randdipscn)))));
randanglebc=90-(acosd(abs((cosd(randdipdirsan-randdipdirscn).*cosd(randdipsan).*cosd(randdipscn))+(sind(randdipsan).*sind(randdipscn)))));


mainspacesets=[spaceseta,spacesetb,spacesetc];
sortmainspacings=sort(mainspacesets);
alphatwo=sortmainspacings(:,2)/sortmainspacings(:,1);
alphathree=sortmainspacings(:,3)/sortmainspacings(:,1);
betashapefactor=round(((alphatwo+alphatwo*alphathree+alphathree)^3)/((alphatwo*alphathree)^2));
set(handles.text_betafromspacings,'Visible','on');
set(handles.text_betafromspacings, 'String', betashapefactor);

%Palmström 1982 Mean Block Volume calculation acc. JV
if get(handles.radio_beta_predefined,'Value') == 1
    beta =str2double(get (handles.edit_blockshapefactor,'String'));
end
if get(handles.radio_beta_spacings,'Value') == 1
    beta =str2double(get (handles.text_betafromspacings,'String'));
end
RESULTS.Rock_characteristics.Rock_Shape_Beta=beta;
%ISBD from drilling cores (unoriented spacing)
Spacingsforscanline=RESULTS.Rock_characteristics.Input_parameters.Spacingsforscanline;
allscanspacings=Spacingsforscanline;
vdrill=beta.*(allscanspacings./2).^3;
isbdvdrill=sort(vdrill);
numvdrill=numel(vdrill);
perc_isbd_vdrill=100/numvdrill:100/numvdrill:100;
drill10=isbdvdrill(1,round(numvdrill*0.1));
drill20=isbdvdrill(1,round(numvdrill*0.2));
drill30=isbdvdrill(1,round(numvdrill*0.3));
drill40=isbdvdrill(1,round(numvdrill*0.4));
drill50=isbdvdrill(1,round(numvdrill*0.5));
drill60=isbdvdrill(1,round(numvdrill*0.6));
drill70=isbdvdrill(1,round(numvdrill*0.7));
drill80=isbdvdrill(1,round(numvdrill*0.8));
drill90=isbdvdrill(1,round(numvdrill*0.9));
drill100=isbdvdrill(1,round(numvdrill));
drilltable=[drill10 drill20 drill30 drill40 drill50 drill60 drill70 drill80 drill90 drill100];




cnegexp=[0.332;0.71;1.207;1.852;2.708;3.98;5.867;8.948;15.332;38.922];
cuniform=[0.375;0.700;1.052;1.460;1.939;2.548;3.343;4.495;6.623;17.772];
cexpuni=[0.420;0.825;1.282;1.824;2.487;3.325;4.439;6.151;9.144;24.905];
clognormal=[0.469;0.949;1.511;2.225;3.094;4.283;5.949;8.498;13.376;38.207];
vnegexp50=2.708*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vuniform50=1.939*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vexpuni50=2.487*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vlognormal50=3.094*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));

vnegexp=cnegexp.*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vuniform=cuniform.*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vexpuni=cexpuni.*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
vlognormal=clognormal.*((spaceseta*spacesetb*spacesetc)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));


cfraccoeff=[0.4649;1.1685;2.1606;3.5458;5.3165;8.0903;13.392;22.607;39.666;108.97];
cfracbi=[0.7882;0.72;0.6719;0.6433;0.644;0.6053;0.5874;0.59;0.5335;0.4675];
vfractal=cfraccoeff.*(spaceseta*spacesetb*spacesetc).^cfracbi;
vfractal50=5.3165*(spaceseta*spacesetb*spacesetc)^0.644;


%Calculates block shape factor ß
jv=RESULTS.Rock_characteristics.Jv;
vol_jv=beta*(jv^-3)/(cosd(angleab)*cosd(angleac)*cosd(anglebc));
set (handles.pushbutton_rmc,'Enable','on');
set(handles.text_volfromjv, 'String', vol_jv);


blockvol=zeros(10,2);
blockvol(:,1)=vnegexp;

blockvol(:,2)=vuniform;

blockvol(:,3)=vexpuni;

blockvol(:,4)=vlognormal;

blockvol(:,5)=vfractal;

blockvol(:,7)=drilltable;

set(handles.text74,'Visible','on');
set(handles.table_blockvolume,'Visible','on');
%set (handles.table_blockvolume,'data',blockvol);


axesblockvol=findobj('Type','axes','Tag','axes_blockvol');
axes(axesblockvol); 
set(handles.textdeviation,'Visible','on');
set(handles.axes_blockvol,'Visible','on');
set(handles.pushbutton_saveblockvolume,'Visible','on');
set(handles.text_spacings_blockvol,'Visible','on');
set(handles.text_angles,'Visible','on');



%end

set(handles.edit_angleab,'Visible','on');
set(handles.edit_angleab, 'String', angleab);
set(handles.edit_angleac,'Visible','on');
set(handles.edit_angleac, 'String', angleac);
set(handles.edit_anglebc,'Visible','on');
set(handles.edit_anglebc, 'String', anglebc);

passing=[10;20;30;40;50;60;70;80;90;100];
angleparameter=sqrt(1+2*sind(angleab)*sind(angleac)*sind(anglebc)-sind(anglebc)^2-sind(angleac)^2-sind(angleab)^2);
vol_palmmean=((spaceseta*spacesetb*spacesetc*angleparameter)/(cosd(angleab)*cosd(angleac)*cosd(anglebc)));
set(handles.text_volpalmmean, 'String', vol_palmmean);

randangleparameter=sqrt(1+2.*sind(randangleab).*sind(randangleac).*sind(randanglebc)-sind(randanglebc).^2-sind(randangleac).^2-sind(randangleab).^2);
palmranda=datasample(spacingsofseta,5000);
palmrandb=datasample(spacingsofsetb,5000);
palmrandc=datasample(spacingsofsetc,5000);
randanglespalm=cosd(randangleab).*cosd(randangleac).*cosd(randanglebc);
%anglespalm=cosd(angleab)*cosd(angleac)*cosd(anglebc);
ibsd_palm=((palmranda.*palmrandb.*palmrandc).*randangleparameter)./(randanglespalm);

sort_ibsd_palm=sort(ibsd_palm);
cumsumisbd=cumsum(sort_ibsd_palm)';
maxisbd=max(cumsumisbd);
[n,p]=size(cumsumisbd);
disp(cumsumisbd)
for i=1:n
    yaxis_palm(1,i)=(cumsumisbd(i,1)./maxisbd)*100;
end
[~,palmy10] = min(abs(yaxis_palm-10));
palm10=sort_ibsd_palm(1,palmy10);
[~,palmy20] = min(abs(yaxis_palm-20));
palm20=sort_ibsd_palm(1,palmy20);
[~,palmy30] = min(abs(yaxis_palm-30));
palm30=sort_ibsd_palm(1,palmy30);
[~,palmy40] = min(abs(yaxis_palm-40));
palm40=sort_ibsd_palm(1,palmy40);
[~,palmy50] = min(abs(yaxis_palm-50));
palm50=sort_ibsd_palm(1,palmy50);
[~,palmy60] = min(abs(yaxis_palm-60));
palm60=sort_ibsd_palm(1,palmy60);
[~,palmy70] = min(abs(yaxis_palm-70));
palm70=sort_ibsd_palm(1,palmy70);
[~,palmy80] = min(abs(yaxis_palm-80));
palm80=sort_ibsd_palm(1,palmy80);
[~,palmy90] = min(abs(yaxis_palm-90));
palm90=sort_ibsd_palm(1,palmy90);
palm100=max(sort_ibsd_palm);

palmtable=[palm10 palm20 palm30 palm40 palm50 palm60 palm70 palm80 palm90 palm100];
blockvol(:,6)=palmtable;
set (handles.table_blockvolume,'data',blockvol);


if get(handles.radio_3negexp,'Value') == 1
    
    semilogx(vnegexp,passing,'LineWidth',3,'Color','blue')
    hold on
    
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    %hold on
    semilogx(sort_ibsd_palm,yaxis_palm,'LineWidth',3,'Color','red')
    legend('3 x Neg. exponential (persistent)','Palmstroem 2005','Location','northwest');
end
if get(handles.radio_3uniform,'Value') == 1
    semilogx(vuniform,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
        semilogx(sort_ibsd_palm,yaxis_palm,'LineWidth',3,'Color','blue')
    legend('3 x Uniform (persistent)','Palmstroem 2005','Location','northwest');

end
if get(handles.radio_2negexp1uni,'Value') == 1
    semilogx(vexpuni,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
        semilogx(sort_ibsd_palm,yaxis_palm,'LineWidth',3,'Color','blue')
    legend('2 x Neg. Exponential, 1 x Uniform (persistent)','Palmstroem 2005','Location','northwest');
end
if get(handles.radio_3lognormal,'Value') == 1
    semilogx(vlognormal,passing,'LineWidth',3,'Color','red') 
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    semilogx(sort_ibsd_palm,yaxis_palm,'LineWidth',3,'Color','blue')
    legend('3 x Lognormal (persistent)','Palmstroem 2005','Location','northwest');
end
if get(handles.radio_3fractal,'Value') == 1
    semilogx(vfractal,passing,'LineWidth',3,'Color','red')
    set(gca,'FontSize',12);
    title('Block size distribution', 'fontsize', 14)
    xlabel('Block size [m³]');
    ylabel('Percentage');
    hold on
    semilogx(sort_ibsd_palm,yaxis_palm,'LineWidth',3,'Color','blue')
    legend('3 x Fractal (persistent)','Palmstroem 2005','Location','northwest');

end
hold off



RESULTS.Block_volume_calculations.ISBDs.Fractal=vfractal;
RESULTS.Block_volume_calculations.ISBDs.Lognormal=vlognormal;
RESULTS.Block_volume_calculations.ISBDs.TwoNegExp_OneUniform=vexpuni;
RESULTS.Block_volume_calculations.ISBDs.Uniform=vuniform;
RESULTS.Block_volume_calculations.ISBDs.Negativeexponential=vnegexp;
RESULTS.Block_volume_calculations.ISBDs.Palmstroem_volume=palmtable;
RESULTS.Block_volume_calculations.ISBDs.Latham_et_al_drillingcores=drilltable;


RESULTS.Block_volume_calculations.Means.Palmstroem_volume=vol_palmmean;
RESULTS.Block_volume_calculations.Means.Jv=vol_jv;




% --- Executes on button press in pushbutton_mergedata.
function pushbutton_mergedata_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_mergedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS idx3
%clear global
filepathsl1=get(handles.edit_data1,'String');
sl1=load(filepathsl1);
filepathsl2=get(handles.edit_data2,'String');
sl2=load(filepathsl2);
sl1length=cat(1,sl1.RESULTS.Measured_area);
sl2length=cat(1,sl2.RESULTS.Measured_area);
RESULTS.Measured_area=sl1length+sl2length;

sl1spacingsscan=cat(1,sl1.RESULTS.Rock_characteristics.Input_parameters.Spacingsforscanline);
sl2spacingsscan=cat(1,sl2.RESULTS.Rock_characteristics.Input_parameters.Spacingsforscanline);
RESULTS.Rock_characteristics.Input_parameters.Spacingsforscanline=[sl1spacingsscan sl2spacingsscan];

%Joint characteristics
sl1terminations=cat(1,sl1.RESULTS.Joint_characteristics.Terminations);
sl2terminations=cat(1,sl2.RESULTS.Joint_characteristics.Terminations);
RESULTS.Joint_characteristics.Terminations=[sl1terminations sl2terminations];
RESULTS.Joint_characteristics.Mean_Termination=mean(RESULTS.Joint_characteristics.Terminations);

sl1TraceLengths=cat(1,sl1.RESULTS.Joint_characteristics.Tracelengths);
sl2TraceLengths=cat(1,sl2.RESULTS.Joint_characteristics.Tracelengths);
RESULTS.Joint_characteristics.Tracelengths=[sl1TraceLengths sl2TraceLengths];
RESULTS.Joint_characteristics.Mean_Tracelength=mean(RESULTS.Joint_characteristics.Tracelengths);
RESULTS.Joint_characteristics.STD_Tracelength=std(RESULTS.Joint_characteristics.Tracelengths);

sl1jointsize=cat(1,sl1.RESULTS.Joint_characteristics.Joint_Size_Factor);
sl2jointsize=cat(1,sl2.RESULTS.Joint_characteristics.Joint_Size_Factor);
RESULTS.Joint_characteristics.Joint_Size_Factor=[sl1jointsize sl2jointsize];
RESULTS.Joint_characteristics.Mean_Joint_Size_Factor=mean(RESULTS.Joint_characteristics.Joint_Size_Factor);

sl1jointalteration=cat(1,sl1.RESULTS.Joint_characteristics.Joint_Alteration_Factor);
sl2jointalteration=cat(1,sl2.RESULTS.Joint_characteristics.Joint_Alteration_Factor);
RESULTS.Joint_characteristics.Joint_Alteration_Factor=[sl1jointalteration sl2jointalteration];
RESULTS.Joint_characteristics.Mean_Joint_Alteration_Factor=mean(RESULTS.Joint_characteristics.Joint_Alteration_Factor);

sl1jointroughness=cat(1,sl1.RESULTS.Joint_characteristics.Joint_Roughness_Factor);
sl2jointroughness=cat(1,sl2.RESULTS.Joint_characteristics.Joint_Roughness_Factor);
RESULTS.Joint_characteristics.Joint_Roughness_Factor=[sl1jointroughness sl2jointroughness];
RESULTS.Joint_characteristics.Mean_Joint_Roughness_Factor=mean(RESULTS.Joint_characteristics.Joint_Roughness_Factor);

sl1Phi_peak_jointsurface=cat(1,sl1.RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles);
sl2Phi_peak_jointsurface=cat(1,sl2.RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles);
RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles=[sl1Phi_peak_jointsurface sl2Phi_peak_jointsurface];
RESULTS.Joint_characteristics.Mean_Peak_Joint_Friction_Angles=mean(RESULTS.Joint_characteristics.Peak_Joint_Friction_Angles);


sl1Joint_Condition_Factor=cat(1,sl1.RESULTS.Joint_characteristics.Mean_Joint_Condition_Factor);
sl2Joint_Condition_Factor=cat(1,sl2.RESULTS.Joint_characteristics.Mean_Joint_Condition_Factor);
sl1observations=numel(sl1jointsize);
sl2observations=numel(sl2jointsize);
weigthedsl1condition=(sl1observations)/(sl1observations+sl2observations)*sl1Joint_Condition_Factor;
weigthedsl2condition=(sl2observations)/(sl1observations+sl2observations)*sl2Joint_Condition_Factor;
RESULTS.Joint_characteristics.Mean_Joint_Condition_Factor=weigthedsl1condition+weigthedsl2condition;



%Set characteristics
%Set1
setnumber=1;
sl1Set1observations=cat(1,sl1.RESULTS.Set1.Number_of_joints_in_set);
sl2Set1observations=cat(1,sl2.RESULTS.Set1.Number_of_joints_in_set);
RESULTS.Set1.Number_of_joints_in_set=[sl1Set1observations sl2Set1observations];
RESULTS.Set1.Number_of_joints_in_set=sum(RESULTS.Set1.Number_of_joints_in_set);

sl1Set1dips=cat(1,sl1.RESULTS.Set1.Dips);
sl2Set1dips=cat(1,sl2.RESULTS.Set1.Dips);
RESULTS.Set1.Dips=[sl1Set1dips sl2Set1dips];

sl1Set1dipdirs=cat(1,sl1.RESULTS.Set1.Dipdirs);
sl2Set1dipdirs=cat(1,sl2.RESULTS.Set1.Dipdirs);
RESULTS.Set1.Dipdirs=[sl1Set1dipdirs sl2Set1dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set1.Dips;
    dipdir=RESULTS.Set1.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set1.Dip=dipeigv;
RESULTS.Set1.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set1.STD_Orientation=stdangleclu;


sl1Set1TraceLengths=cat(1,sl1.RESULTS.Set1.Tracelengths);
sl2Set1TraceLengths=cat(1,sl2.RESULTS.Set1.Tracelengths);
RESULTS.Set1.Tracelengths=[sl1Set1TraceLengths sl2Set1TraceLengths];
RESULTS.Set1.Mean_Tracelength=mean(RESULTS.Set1.Tracelengths);
RESULTS.Set1.STD_Tracelength=std(RESULTS.Set1.Tracelengths);

sl1Set1Joint_Condition_Factor=cat(1,sl1.RESULTS.Set1.Mean_Joint_Condition_Factor);
sl2Set1Joint_Condition_Factor=cat(1,sl2.RESULTS.Set1.Mean_Joint_Condition_Factor);
weigthedsl1Set1condition=(sl1Set1observations)/(sl1Set1observations+sl2Set1observations)*sl1Set1Joint_Condition_Factor;
weigthedsl2Set1condition=(sl2Set1observations)/(sl1Set1observations+sl2Set1observations)*sl2Set1Joint_Condition_Factor;
RESULTS.Set1.Mean_Joint_Condition_Factor=weigthedsl1Set1condition+weigthedsl2Set1condition;


sl1Set1Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set1.Peak_Friction_Angles);
sl2Set1Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set1.Peak_Friction_Angles);
RESULTS.Set1.Peak_Friction_Angles=[sl1Set1Phi_peak_jointsurface sl2Set1Phi_peak_jointsurface];
RESULTS.Set1.Mean_Peak_Friction_Angle=mean(RESULTS.Set1.Peak_Friction_Angles);
RESULTS.Set1.STD_Peak_Friction_Angle=std(RESULTS.Set1.Peak_Friction_Angles);

sl1Set1jointalteration=cat(1,sl1.RESULTS.Set1.Joint_Alteration_Factors);
sl2Set1jointalteration=cat(1,sl2.RESULTS.Set1.Joint_Alteration_Factors);
RESULTS.Set1.Joint_Alteration_Factors=[sl1Set1jointalteration sl2Set1jointalteration];
RESULTS.Set1.Mean_Joint_Alteration_Factor=mean(RESULTS.Set1.Joint_Alteration_Factors);
RESULTS.Set1.STD_Joint_Alteration_Factor=std(RESULTS.Set1.Joint_Alteration_Factors);

sl1Set1jointroughness=cat(1,sl1.RESULTS.Set1.Joint_Roughness_Factors);
sl2Set1jointroughness=cat(1,sl2.RESULTS.Set1.Joint_Roughness_Factors);
RESULTS.Set1.Joint_Roughness_Factors=[sl1Set1jointroughness sl2Set1jointroughness];
RESULTS.Set1.Mean_Joint_Roughness_Factor=mean(RESULTS.Set1.Joint_Roughness_Factors);
RESULTS.Set1.STD_Joint_Roughness_Factor=std(RESULTS.Set1.Joint_Roughness_Factors);

sl1Set1jointsize=cat(1,sl1.RESULTS.Set1.Joint_Size_Factor);
sl2Set1jointsize=cat(1,sl2.RESULTS.Set1.Joint_Size_Factor);
RESULTS.Set1.Joint_Size_Factor=[sl1Set1jointsize sl2Set1jointsize];
RESULTS.Set1.Mean_Joint_Size_Factor=mean(RESULTS.Set1.Joint_Size_Factor);
RESULTS.Set1.STD_Joint_Size_Factor=std(RESULTS.Set1.Joint_Size_Factor);

sl1Set1Intersections=cat(1,sl1.RESULTS.Set1.Intersections);
sl2Set1Intersections=cat(1,sl2.RESULTS.Set1.Intersections);
RESULTS.Set1.Intersections=[sl1Set1Intersections sl2Set1Intersections];

sl1Set1terminations=cat(1,sl1.RESULTS.Set1.Terminations);
sl2Set1terminations=cat(1,sl2.RESULTS.Set1.Terminations);
RESULTS.Set1.Terminations=[sl1Set1terminations sl2Set1terminations];
RESULTS.Set1.Mean_Termination=mean(RESULTS.Set1.Terminations);

sl1Set1termination=cat(1,sl1.RESULTS.Set1.Termination_index);
sl2Set1termination=cat(1,sl2.RESULTS.Set1.Termination_index);
weigthedsl1Set1termination=(sl1Set1observations)/(sl1Set1observations+sl2Set1observations)*sl1Set1termination;
weigthedsl2Set1termination=(sl2Set1observations)/(sl1Set1observations+sl2Set1observations)*sl2Set1termination;
RESULTS.Set1.Termination_index=weigthedsl1Set1termination+weigthedsl2Set1termination;


sl1Set1DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set1.Normal_spacings);
sl2Set1DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set1.Normal_spacings);
RESULTS.Set1.Normal_spacings=[sl1Set1DiscontinuityNormalSpacings sl2Set1DiscontinuityNormalSpacings];
RESULTS.Set1.Mean_normal_spacing=mean(RESULTS.Set1.Normal_spacings);
RESULTS.Set1.Median_normal_spacing=median(RESULTS.Set1.Normal_spacings);
RESULTS.Set1.Mode_normal_spacing=mode(RESULTS.Set1.Normal_spacings);
RESULTS.Set1.STD_normal_spacing=std(RESULTS.Set1.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set1.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set1observations)/(sl1Set1observations+sl2Set1observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set1.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set1observations)/(sl1Set1observations+sl2Set1observations)*sl2MeanLaslettTL;
RESULTS.Set1.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set1.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set1observations)/(sl1Set1observations+sl2Set1observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set1.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set1observations)/(sl1Set1observations+sl2Set1observations)*sl2MeanPriest93TL;
RESULTS.Set1.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set1.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set1observations)/(sl1Set1observations+sl2Set1observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set1.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set1observations)/(sl1Set1observations+sl2Set1observations)*sl2MeanPriest81TL;
RESULTS.Set1.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;

%Set2
if isfield(sl1.RESULTS,'Set2')==1
setnumber=2;

sl1Set2observations=cat(1,sl1.RESULTS.Set2.Number_of_joints_in_set);
sl2Set2observations=cat(1,sl2.RESULTS.Set2.Number_of_joints_in_set);
RESULTS.Set2.Number_of_joints_in_set=[sl1Set2observations sl2Set2observations];
RESULTS.Set2.Number_of_joints_in_set=sum(RESULTS.Set2.Number_of_joints_in_set);

sl1Set2dips=cat(1,sl1.RESULTS.Set2.Dips);
sl2Set2dips=cat(1,sl2.RESULTS.Set2.Dips);
RESULTS.Set2.Dips=[sl1Set2dips sl2Set2dips];

sl1Set2dipdirs=cat(1,sl1.RESULTS.Set2.Dipdirs);
sl2Set2dipdirs=cat(1,sl2.RESULTS.Set2.Dipdirs);
RESULTS.Set2.Dipdirs=[sl1Set2dipdirs sl2Set2dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set2.Dips;
    dipdir=RESULTS.Set2.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set2.Dip=dipeigv;
RESULTS.Set2.Dipdir=dipdireigv;

    %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set2.STD_Orientation=stdangleclu;
 



sl1Set2TraceLengths=cat(1,sl1.RESULTS.Set2.Tracelengths);
sl2Set2TraceLengths=cat(1,sl2.RESULTS.Set2.Tracelengths);
RESULTS.Set2.Tracelengths=[sl1Set2TraceLengths sl2Set2TraceLengths];
RESULTS.Set2.Mean_Tracelength=mean(RESULTS.Set2.Tracelengths);
RESULTS.Set2.STD_Tracelength=std(RESULTS.Set2.Tracelengths);

sl1Set2Joint_Condition_Factor=cat(1,sl1.RESULTS.Set2.Mean_Joint_Condition_Factor);
sl2Set2Joint_Condition_Factor=cat(1,sl2.RESULTS.Set2.Mean_Joint_Condition_Factor);
weigthedsl1Set2condition=(sl1Set2observations)/(sl1Set2observations+sl2Set2observations)*sl1Set2Joint_Condition_Factor;
weigthedsl2Set2condition=(sl2Set2observations)/(sl1Set2observations+sl2Set2observations)*sl2Set2Joint_Condition_Factor;
RESULTS.Set2.Mean_Joint_Condition_Factor=weigthedsl1Set2condition+weigthedsl2Set2condition;


sl1Set2Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set2.Peak_Friction_Angles);
sl2Set2Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set2.Peak_Friction_Angles);
RESULTS.Set2.Peak_Friction_Angles=[sl1Set2Phi_peak_jointsurface sl2Set2Phi_peak_jointsurface];
RESULTS.Set2.Mean_Peak_Friction_Angle=mean(RESULTS.Set2.Peak_Friction_Angles);
RESULTS.Set2.STD_Peak_Friction_Angle=std(RESULTS.Set2.Peak_Friction_Angles);

sl1Set2jointalteration=cat(1,sl1.RESULTS.Set2.Joint_Alteration_Factors);
sl2Set2jointalteration=cat(1,sl2.RESULTS.Set2.Joint_Alteration_Factors);
RESULTS.Set2.Joint_Alteration_Factors=[sl1Set2jointalteration sl2Set2jointalteration];
RESULTS.Set2.Mean_Joint_Alteration_Factor=mean(RESULTS.Set2.Joint_Alteration_Factors);
RESULTS.Set2.STD_Joint_Alteration_Factor=std(RESULTS.Set2.Joint_Alteration_Factors);

sl1Set2jointroughness=cat(1,sl1.RESULTS.Set2.Joint_Roughness_Factors);
sl2Set2jointroughness=cat(1,sl2.RESULTS.Set2.Joint_Roughness_Factors);
RESULTS.Set2.Joint_Roughness_Factors=[sl1Set2jointroughness sl2Set2jointroughness];
RESULTS.Set2.Mean_Joint_Roughness_Factor=mean(RESULTS.Set2.Joint_Roughness_Factors);
RESULTS.Set2.STD_Joint_Roughness_Factor=std(RESULTS.Set2.Joint_Roughness_Factors);

sl1Set2jointsize=cat(1,sl1.RESULTS.Set2.Joint_Size_Factor);
sl2Set2jointsize=cat(1,sl2.RESULTS.Set2.Joint_Size_Factor);
RESULTS.Set2.Joint_Size_Factor=[sl1Set2jointsize sl2Set2jointsize];
RESULTS.Set2.Mean_Joint_Size_Factor=mean(RESULTS.Set2.Joint_Size_Factor);
RESULTS.Set2.STD_Joint_Size_Factor=std(RESULTS.Set2.Joint_Size_Factor);

sl1Set2Intersections=cat(1,sl1.RESULTS.Set2.Intersections);
sl2Set2Intersections=cat(1,sl2.RESULTS.Set2.Intersections);
RESULTS.Set2.Intersections=[sl1Set2Intersections sl2Set2Intersections];

sl1Set2terminations=cat(1,sl1.RESULTS.Set2.Terminations);
sl2Set2terminations=cat(1,sl2.RESULTS.Set2.Terminations);
RESULTS.Set2.Terminations=[sl1Set2terminations sl2Set2terminations];
RESULTS.Set2.Mean_Termination=mean(RESULTS.Set2.Terminations);

sl1Set2termination=cat(1,sl1.RESULTS.Set2.Termination_index);
sl2Set2termination=cat(1,sl2.RESULTS.Set2.Termination_index);
weigthedsl1Set2termination=(sl1Set2observations)/(sl1Set2observations+sl2Set2observations)*sl1Set2termination;
weigthedsl2Set2termination=(sl2Set2observations)/(sl1Set2observations+sl2Set2observations)*sl2Set2termination;
RESULTS.Set2.Termination_index=weigthedsl1Set2termination+weigthedsl2Set2termination;


sl1Set2DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set2.Normal_spacings);
sl2Set2DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set2.Normal_spacings);
RESULTS.Set2.Normal_spacings=[sl1Set2DiscontinuityNormalSpacings sl2Set2DiscontinuityNormalSpacings];
RESULTS.Set2.Mean_normal_spacing=mean(RESULTS.Set2.Normal_spacings);
RESULTS.Set2.Median_normal_spacing=median(RESULTS.Set2.Normal_spacings);
RESULTS.Set2.Mode_normal_spacing=mode(RESULTS.Set2.Normal_spacings);
RESULTS.Set2.STD_normal_spacing=std(RESULTS.Set2.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set2.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set2observations)/(sl1Set2observations+sl2Set2observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set2.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set2observations)/(sl1Set2observations+sl2Set2observations)*sl2MeanLaslettTL;
RESULTS.Set2.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set2.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set2observations)/(sl1Set2observations+sl2Set2observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set2.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set2observations)/(sl1Set2observations+sl2Set2observations)*sl2MeanPriest93TL;
RESULTS.Set2.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set2.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set2observations)/(sl1Set2observations+sl2Set2observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set2.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set2observations)/(sl1Set2observations+sl2Set2observations)*sl2MeanPriest81TL;
RESULTS.Set2.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end

%Set3
if isfield(sl1.RESULTS,'Set3')==1
setnumber=3;
sl1Set3observations=cat(1,sl1.RESULTS.Set3.Number_of_joints_in_set);
sl2Set3observations=cat(1,sl2.RESULTS.Set3.Number_of_joints_in_set);
RESULTS.Set3.Number_of_joints_in_set=[sl1Set3observations sl2Set3observations];
RESULTS.Set3.Number_of_joints_in_set=sum(RESULTS.Set3.Number_of_joints_in_set);

sl1Set3dips=cat(1,sl1.RESULTS.Set3.Dips);
sl2Set3dips=cat(1,sl2.RESULTS.Set3.Dips);
RESULTS.Set3.Dips=[sl1Set3dips sl2Set3dips];

sl1Set3dipdirs=cat(1,sl1.RESULTS.Set3.Dipdirs);
sl2Set3dipdirs=cat(1,sl2.RESULTS.Set3.Dipdirs);
RESULTS.Set3.Dipdirs=[sl1Set3dipdirs sl2Set3dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set3.Dips;
    dipdir=RESULTS.Set3.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set3.Dip=dipeigv;
RESULTS.Set3.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set3.STD_Orientation=stdangleclu;

sl1Set3TraceLengths=cat(1,sl1.RESULTS.Set3.Tracelengths);
sl2Set3TraceLengths=cat(1,sl2.RESULTS.Set3.Tracelengths);
RESULTS.Set3.Tracelengths=[sl1Set3TraceLengths sl2Set3TraceLengths];
RESULTS.Set3.Mean_Tracelength=mean(RESULTS.Set3.Tracelengths);
RESULTS.Set3.STD_Tracelength=std(RESULTS.Set3.Tracelengths);

sl1Set3Joint_Condition_Factor=cat(1,sl1.RESULTS.Set3.Mean_Joint_Condition_Factor);
sl2Set3Joint_Condition_Factor=cat(1,sl2.RESULTS.Set3.Mean_Joint_Condition_Factor);
weigthedsl1Set3condition=(sl1Set3observations)/(sl1Set3observations+sl2Set3observations)*sl1Set3Joint_Condition_Factor;
weigthedsl2Set3condition=(sl2Set3observations)/(sl1Set3observations+sl2Set3observations)*sl2Set3Joint_Condition_Factor;
RESULTS.Set3.Mean_Joint_Condition_Factor=weigthedsl1Set3condition+weigthedsl2Set3condition;


sl1Set3Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set3.Peak_Friction_Angles);
sl2Set3Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set3.Peak_Friction_Angles);
RESULTS.Set3.Peak_Friction_Angles=[sl1Set3Phi_peak_jointsurface sl2Set3Phi_peak_jointsurface];
RESULTS.Set3.Mean_Peak_Friction_Angle=mean(RESULTS.Set3.Peak_Friction_Angles);
RESULTS.Set3.STD_Peak_Friction_Angle=std(RESULTS.Set3.Peak_Friction_Angles);

sl1Set3jointalteration=cat(1,sl1.RESULTS.Set3.Joint_Alteration_Factors);
sl2Set3jointalteration=cat(1,sl2.RESULTS.Set3.Joint_Alteration_Factors);
RESULTS.Set3.Joint_Alteration_Factors=[sl1Set3jointalteration sl2Set3jointalteration];
RESULTS.Set3.Mean_Joint_Alteration_Factor=mean(RESULTS.Set3.Joint_Alteration_Factors);
RESULTS.Set3.STD_Joint_Alteration_Factor=std(RESULTS.Set3.Joint_Alteration_Factors);

sl1Set3jointroughness=cat(1,sl1.RESULTS.Set3.Joint_Roughness_Factors);
sl2Set3jointroughness=cat(1,sl2.RESULTS.Set3.Joint_Roughness_Factors);
RESULTS.Set3.Joint_Roughness_Factors=[sl1Set3jointroughness sl2Set3jointroughness];
RESULTS.Set3.Mean_Joint_Roughness_Factor=mean(RESULTS.Set3.Joint_Roughness_Factors);
RESULTS.Set3.STD_Joint_Roughness_Factor=std(RESULTS.Set3.Joint_Roughness_Factors);

sl1Set3jointsize=cat(1,sl1.RESULTS.Set3.Joint_Size_Factor);
sl2Set3jointsize=cat(1,sl2.RESULTS.Set3.Joint_Size_Factor);
RESULTS.Set3.Joint_Size_Factor=[sl1Set3jointsize sl2Set3jointsize];
RESULTS.Set3.Mean_Joint_Size_Factor=mean(RESULTS.Set3.Joint_Size_Factor);
RESULTS.Set3.STD_Joint_Size_Factor=std(RESULTS.Set3.Joint_Size_Factor);

sl1Set3Intersections=cat(1,sl1.RESULTS.Set3.Intersections);
sl2Set3Intersections=cat(1,sl2.RESULTS.Set3.Intersections);
RESULTS.Set3.Intersections=[sl1Set3Intersections sl2Set3Intersections];

sl1Set3terminations=cat(1,sl1.RESULTS.Set3.Terminations);
sl2Set3terminations=cat(1,sl2.RESULTS.Set3.Terminations);
RESULTS.Set3.Terminations=[sl1Set3terminations sl2Set3terminations];
RESULTS.Set3.Mean_Termination=mean(RESULTS.Set3.Terminations);

sl1Set3termination=cat(1,sl1.RESULTS.Set3.Termination_index);
sl2Set3termination=cat(1,sl2.RESULTS.Set3.Termination_index);
weigthedsl1Set3termination=(sl1Set3observations)/(sl1Set3observations+sl2Set3observations)*sl1Set3termination;
weigthedsl2Set3termination=(sl2Set3observations)/(sl1Set3observations+sl2Set3observations)*sl2Set3termination;
RESULTS.Set3.Termination_index=weigthedsl1Set3termination+weigthedsl2Set3termination;


sl1Set3DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set3.Normal_spacings);
sl2Set3DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set3.Normal_spacings);
RESULTS.Set3.Normal_spacings=[sl1Set3DiscontinuityNormalSpacings sl2Set3DiscontinuityNormalSpacings];
RESULTS.Set3.Mean_normal_spacing=mean(RESULTS.Set3.Normal_spacings);
RESULTS.Set3.Median_normal_spacing=median(RESULTS.Set3.Normal_spacings);
RESULTS.Set3.Mode_normal_spacing=mode(RESULTS.Set3.Normal_spacings);
RESULTS.Set3.STD_normal_spacing=std(RESULTS.Set3.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set3.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set3observations)/(sl1Set3observations+sl2Set3observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set3.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set3observations)/(sl1Set3observations+sl2Set3observations)*sl2MeanLaslettTL;
RESULTS.Set3.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set3.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set3observations)/(sl1Set3observations+sl2Set3observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set3.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set3observations)/(sl1Set3observations+sl2Set3observations)*sl2MeanPriest93TL;
RESULTS.Set3.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set3.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set3observations)/(sl1Set3observations+sl2Set3observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set3.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set3observations)/(sl1Set3observations+sl2Set3observations)*sl2MeanPriest81TL;
RESULTS.Set3.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end

%Set4
if isfield(sl1.RESULTS,'Set4')==1
setnumber=4;
sl1Set4observations=cat(1,sl1.RESULTS.Set4.Number_of_joints_in_set);
sl2Set4observations=cat(1,sl2.RESULTS.Set4.Number_of_joints_in_set);
RESULTS.Set4.Number_of_joints_in_set=[sl1Set4observations sl2Set4observations];
RESULTS.Set4.Number_of_joints_in_set=sum(RESULTS.Set4.Number_of_joints_in_set);

sl1Set4dips=cat(1,sl1.RESULTS.Set4.Dips);
sl2Set4dips=cat(1,sl2.RESULTS.Set4.Dips);
RESULTS.Set4.Dips=[sl1Set4dips sl2Set4dips];

sl1Set4dipdirs=cat(1,sl1.RESULTS.Set4.Dipdirs);
sl2Set4dipdirs=cat(1,sl2.RESULTS.Set4.Dipdirs);
RESULTS.Set4.Dipdirs=[sl1Set4dipdirs sl2Set4dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set4.Dips;
    dipdir=RESULTS.Set4.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set4.Dip=dipeigv;
RESULTS.Set4.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set4.STD_Orientation=stdangleclu;

sl1Set4TraceLengths=cat(1,sl1.RESULTS.Set4.Tracelengths);
sl2Set4TraceLengths=cat(1,sl2.RESULTS.Set4.Tracelengths);
RESULTS.Set4.Tracelengths=[sl1Set4TraceLengths sl2Set4TraceLengths];
RESULTS.Set4.Mean_Tracelength=mean(RESULTS.Set4.Tracelengths);
RESULTS.Set4.STD_Tracelength=std(RESULTS.Set4.Tracelengths);

sl1Set4Joint_Condition_Factor=cat(1,sl1.RESULTS.Set4.Mean_Joint_Condition_Factor);
sl2Set4Joint_Condition_Factor=cat(1,sl2.RESULTS.Set4.Mean_Joint_Condition_Factor);
weigthedsl1Set4condition=(sl1Set4observations)/(sl1Set4observations+sl2Set4observations)*sl1Set4Joint_Condition_Factor;
weigthedsl2Set4condition=(sl2Set4observations)/(sl1Set4observations+sl2Set4observations)*sl2Set4Joint_Condition_Factor;
RESULTS.Set4.Mean_Joint_Condition_Factor=weigthedsl1Set4condition+weigthedsl2Set4condition;


sl1Set4Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set4.Peak_Friction_Angles);
sl2Set4Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set4.Peak_Friction_Angles);
RESULTS.Set4.Peak_Friction_Angles=[sl1Set4Phi_peak_jointsurface sl2Set4Phi_peak_jointsurface];
RESULTS.Set4.Mean_Peak_Friction_Angle=mean(RESULTS.Set4.Peak_Friction_Angles);
RESULTS.Set4.STD_Peak_Friction_Angle=std(RESULTS.Set4.Peak_Friction_Angles);

sl1Set4jointalteration=cat(1,sl1.RESULTS.Set4.Joint_Alteration_Factors);
sl2Set4jointalteration=cat(1,sl2.RESULTS.Set4.Joint_Alteration_Factors);
RESULTS.Set4.Joint_Alteration_Factors=[sl1Set4jointalteration sl2Set4jointalteration];
RESULTS.Set4.Mean_Joint_Alteration_Factor=mean(RESULTS.Set4.Joint_Alteration_Factors);
RESULTS.Set4.STD_Joint_Alteration_Factor=std(RESULTS.Set4.Joint_Alteration_Factors);

sl1Set4jointroughness=cat(1,sl1.RESULTS.Set4.Joint_Roughness_Factors);
sl2Set4jointroughness=cat(1,sl2.RESULTS.Set4.Joint_Roughness_Factors);
RESULTS.Set4.Joint_Roughness_Factors=[sl1Set4jointroughness sl2Set4jointroughness];
RESULTS.Set4.Mean_Joint_Roughness_Factor=mean(RESULTS.Set4.Joint_Roughness_Factors);
RESULTS.Set4.STD_Joint_Roughness_Factor=std(RESULTS.Set4.Joint_Roughness_Factors);

sl1Set4jointsize=cat(1,sl1.RESULTS.Set4.Joint_Size_Factor);
sl2Set4jointsize=cat(1,sl2.RESULTS.Set4.Joint_Size_Factor);
RESULTS.Set4.Joint_Size_Factor=[sl1Set4jointsize sl2Set4jointsize];
RESULTS.Set4.Mean_Joint_Size_Factor=mean(RESULTS.Set4.Joint_Size_Factor);
RESULTS.Set4.STD_Joint_Size_Factor=std(RESULTS.Set4.Joint_Size_Factor);

sl1Set4Intersections=cat(1,sl1.RESULTS.Set4.Intersections);
sl2Set4Intersections=cat(1,sl2.RESULTS.Set4.Intersections);
RESULTS.Set4.Intersections=[sl1Set4Intersections sl2Set4Intersections];

sl1Set4terminations=cat(1,sl1.RESULTS.Set4.Terminations);
sl2Set4terminations=cat(1,sl2.RESULTS.Set4.Terminations);
RESULTS.Set4.Terminations=[sl1Set4terminations sl2Set4terminations];
RESULTS.Set4.Mean_Termination=mean(RESULTS.Set4.Terminations);

sl1Set4termination=cat(1,sl1.RESULTS.Set4.Termination_index);
sl2Set4termination=cat(1,sl2.RESULTS.Set4.Termination_index);
weigthedsl1Set4termination=(sl1Set4observations)/(sl1Set4observations+sl2Set4observations)*sl1Set4termination;
weigthedsl2Set4termination=(sl2Set4observations)/(sl1Set4observations+sl2Set4observations)*sl2Set4termination;
RESULTS.Set4.Termination_index=weigthedsl1Set4termination+weigthedsl2Set4termination;


sl1Set4DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set4.Normal_spacings);
sl2Set4DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set4.Normal_spacings);
RESULTS.Set4.Normal_spacings=[sl1Set4DiscontinuityNormalSpacings sl2Set4DiscontinuityNormalSpacings];
RESULTS.Set4.Mean_normal_spacing=mean(RESULTS.Set4.Normal_spacings);
RESULTS.Set4.Median_normal_spacing=median(RESULTS.Set4.Normal_spacings);
RESULTS.Set4.Mode_normal_spacing=mode(RESULTS.Set4.Normal_spacings);
RESULTS.Set4.STD_normal_spacing=std(RESULTS.Set4.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set4.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set4observations)/(sl1Set4observations+sl2Set4observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set4.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set4observations)/(sl1Set4observations+sl2Set4observations)*sl2MeanLaslettTL;
RESULTS.Set4.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set4.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set4observations)/(sl1Set4observations+sl2Set4observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set4.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set4observations)/(sl1Set4observations+sl2Set4observations)*sl2MeanPriest93TL;
RESULTS.Set4.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set4.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set4observations)/(sl1Set4observations+sl2Set4observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set4.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set4observations)/(sl1Set4observations+sl2Set4observations)*sl2MeanPriest81TL;
RESULTS.Set4.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end

%Set5
  
if isfield(sl1.RESULTS,'Set5')==1
setnumber=5;
sl1Set5observations=cat(1,sl1.RESULTS.Set5.Number_of_joints_in_set);
sl2Set5observations=cat(1,sl2.RESULTS.Set5.Number_of_joints_in_set);
RESULTS.Set5.Number_of_joints_in_set=[sl1Set5observations sl2Set5observations];
RESULTS.Set5.Number_of_joints_in_set=sum(RESULTS.Set5.Number_of_joints_in_set);

sl1Set5dips=cat(1,sl1.RESULTS.Set5.Dips);
sl2Set5dips=cat(1,sl2.RESULTS.Set5.Dips);
RESULTS.Set5.Dips=[sl1Set5dips sl2Set5dips];

sl1Set5dipdirs=cat(1,sl1.RESULTS.Set5.Dipdirs);
sl2Set5dipdirs=cat(1,sl2.RESULTS.Set5.Dipdirs);
RESULTS.Set5.Dipdirs=[sl1Set5dipdirs sl2Set5dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set5.Dips;
    dipdir=RESULTS.Set5.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set5.Dip=dipeigv;
RESULTS.Set5.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set5.STD_Orientation=stdangleclu;

sl1Set5TraceLengths=cat(1,sl1.RESULTS.Set5.Tracelengths);
sl2Set5TraceLengths=cat(1,sl2.RESULTS.Set5.Tracelengths);
RESULTS.Set5.Tracelengths=[sl1Set5TraceLengths sl2Set5TraceLengths];
RESULTS.Set5.Mean_Tracelength=mean(RESULTS.Set5.Tracelengths);
RESULTS.Set5.STD_Tracelength=std(RESULTS.Set5.Tracelengths);

sl1Set5Joint_Condition_Factor=cat(1,sl1.RESULTS.Set5.Mean_Joint_Condition_Factor);
sl2Set5Joint_Condition_Factor=cat(1,sl2.RESULTS.Set5.Mean_Joint_Condition_Factor);
weigthedsl1Set5condition=(sl1Set5observations)/(sl1Set5observations+sl2Set5observations)*sl1Set5Joint_Condition_Factor;
weigthedsl2Set5condition=(sl2Set5observations)/(sl1Set5observations+sl2Set5observations)*sl2Set5Joint_Condition_Factor;
RESULTS.Set5.Mean_Joint_Condition_Factor=weigthedsl1Set5condition+weigthedsl2Set5condition;


sl1Set5Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set5.Peak_Friction_Angles);
sl2Set5Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set5.Peak_Friction_Angles);
RESULTS.Set5.Peak_Friction_Angles=[sl1Set5Phi_peak_jointsurface sl2Set5Phi_peak_jointsurface];
RESULTS.Set5.Mean_Peak_Friction_Angle=mean(RESULTS.Set5.Peak_Friction_Angles);
RESULTS.Set5.STD_Peak_Friction_Angle=std(RESULTS.Set5.Peak_Friction_Angles);

sl1Set5jointalteration=cat(1,sl1.RESULTS.Set5.Joint_Alteration_Factors);
sl2Set5jointalteration=cat(1,sl2.RESULTS.Set5.Joint_Alteration_Factors);
RESULTS.Set5.Joint_Alteration_Factors=[sl1Set5jointalteration sl2Set5jointalteration];
RESULTS.Set5.Mean_Joint_Alteration_Factor=mean(RESULTS.Set5.Joint_Alteration_Factors);
RESULTS.Set5.STD_Joint_Alteration_Factor=std(RESULTS.Set5.Joint_Alteration_Factors);

sl1Set5jointroughness=cat(1,sl1.RESULTS.Set5.Joint_Roughness_Factors);
sl2Set5jointroughness=cat(1,sl2.RESULTS.Set5.Joint_Roughness_Factors);
RESULTS.Set5.Joint_Roughness_Factors=[sl1Set5jointroughness sl2Set5jointroughness];
RESULTS.Set5.Mean_Joint_Roughness_Factor=mean(RESULTS.Set5.Joint_Roughness_Factors);
RESULTS.Set5.STD_Joint_Roughness_Factor=std(RESULTS.Set5.Joint_Roughness_Factors);

sl1Set5jointsize=cat(1,sl1.RESULTS.Set5.Joint_Size_Factor);
sl2Set5jointsize=cat(1,sl2.RESULTS.Set5.Joint_Size_Factor);
RESULTS.Set5.Joint_Size_Factor=[sl1Set5jointsize sl2Set5jointsize];
RESULTS.Set5.Mean_Joint_Size_Factor=mean(RESULTS.Set5.Joint_Size_Factor);
RESULTS.Set5.STD_Joint_Size_Factor=std(RESULTS.Set5.Joint_Size_Factor);

sl1Set5Intersections=cat(1,sl1.RESULTS.Set5.Intersections);
sl2Set5Intersections=cat(1,sl2.RESULTS.Set5.Intersections);
RESULTS.Set5.Intersections=[sl1Set5Intersections sl2Set5Intersections];

sl1Set5terminations=cat(1,sl1.RESULTS.Set5.Terminations);
sl2Set5terminations=cat(1,sl2.RESULTS.Set5.Terminations);
RESULTS.Set5.Terminations=[sl1Set5terminations sl2Set5terminations];
RESULTS.Set5.Mean_Termination=mean(RESULTS.Set5.Terminations);

sl1Set5termination=cat(1,sl1.RESULTS.Set5.Termination_index);
sl2Set5termination=cat(1,sl2.RESULTS.Set5.Termination_index);
weigthedsl1Set5termination=(sl1Set5observations)/(sl1Set5observations+sl2Set5observations)*sl1Set5termination;
weigthedsl2Set5termination=(sl2Set5observations)/(sl1Set5observations+sl2Set5observations)*sl2Set5termination;
RESULTS.Set5.Termination_index=weigthedsl1Set5termination+weigthedsl2Set5termination;


sl1Set5DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set5.Normal_spacings);
sl2Set5DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set5.Normal_spacings);
RESULTS.Set5.Normal_spacings=[sl1Set5DiscontinuityNormalSpacings sl2Set5DiscontinuityNormalSpacings];
RESULTS.Set5.Mean_normal_spacing=mean(RESULTS.Set5.Normal_spacings);
RESULTS.Set5.Median_normal_spacing=median(RESULTS.Set5.Normal_spacings);
RESULTS.Set5.Mode_normal_spacing=mode(RESULTS.Set5.Normal_spacings);
RESULTS.Set5.STD_normal_spacing=std(RESULTS.Set5.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set5.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set5observations)/(sl1Set5observations+sl2Set5observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set5.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set5observations)/(sl1Set5observations+sl2Set5observations)*sl2MeanLaslettTL;
RESULTS.Set5.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set5.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set5observations)/(sl1Set5observations+sl2Set5observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set5.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set5observations)/(sl1Set5observations+sl2Set5observations)*sl2MeanPriest93TL;
RESULTS.Set5.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set5.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set5observations)/(sl1Set5observations+sl2Set5observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set5.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set5observations)/(sl1Set5observations+sl2Set5observations)*sl2MeanPriest81TL;
RESULTS.Set5.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end

%Set6
if isfield(sl1.RESULTS,'Set6')==1
setnumber=6;  
sl1Set6observations=cat(1,sl1.RESULTS.Set6.Number_of_joints_in_set);
sl2Set6observations=cat(1,sl2.RESULTS.Set6.Number_of_joints_in_set);
RESULTS.Set6.Number_of_joints_in_set=[sl1Set6observations sl2Set6observations];
RESULTS.Set6.Number_of_joints_in_set=sum(RESULTS.Set6.Number_of_joints_in_set);

sl1Set6dips=cat(1,sl1.RESULTS.Set6.Dips);
sl2Set6dips=cat(1,sl2.RESULTS.Set6.Dips);
RESULTS.Set6.Dips=[sl1Set6dips sl2Set6dips];

sl1Set6dipdirs=cat(1,sl1.RESULTS.Set6.Dipdirs);
sl2Set6dipdirs=cat(1,sl2.RESULTS.Set6.Dipdirs);
RESULTS.Set6.Dipdirs=[sl1Set6dipdirs sl2Set6dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set6.Dips;
    dipdir=RESULTS.Set6.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set6.Dip=dipeigv;
RESULTS.Set6.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set6.STD_Orientation=stdangleclu;

sl1Set6TraceLengths=cat(1,sl1.RESULTS.Set6.Tracelengths);
sl2Set6TraceLengths=cat(1,sl2.RESULTS.Set6.Tracelengths);
RESULTS.Set6.Tracelengths=[sl1Set6TraceLengths sl2Set6TraceLengths];
RESULTS.Set6.Mean_Tracelength=mean(RESULTS.Set6.Tracelengths);
RESULTS.Set6.STD_Tracelength=std(RESULTS.Set6.Tracelengths);

sl1Set6Joint_Condition_Factor=cat(1,sl1.RESULTS.Set6.Mean_Joint_Condition_Factor);
sl2Set6Joint_Condition_Factor=cat(1,sl2.RESULTS.Set6.Mean_Joint_Condition_Factor);
weigthedsl1Set6condition=(sl1Set6observations)/(sl1Set6observations+sl2Set6observations)*sl1Set6Joint_Condition_Factor;
weigthedsl2Set6condition=(sl2Set6observations)/(sl1Set6observations+sl2Set6observations)*sl2Set6Joint_Condition_Factor;
RESULTS.Set6.Mean_Joint_Condition_Factor=weigthedsl1Set6condition+weigthedsl2Set6condition;


sl1Set6Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set6.Peak_Friction_Angles);
sl2Set6Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set6.Peak_Friction_Angles);
RESULTS.Set6.Peak_Friction_Angles=[sl1Set6Phi_peak_jointsurface sl2Set6Phi_peak_jointsurface];
RESULTS.Set6.Mean_Peak_Friction_Angle=mean(RESULTS.Set6.Peak_Friction_Angles);
RESULTS.Set6.STD_Peak_Friction_Angle=std(RESULTS.Set6.Peak_Friction_Angles);

sl1Set6jointalteration=cat(1,sl1.RESULTS.Set6.Joint_Alteration_Factors);
sl2Set6jointalteration=cat(1,sl2.RESULTS.Set6.Joint_Alteration_Factors);
RESULTS.Set6.Joint_Alteration_Factors=[sl1Set6jointalteration sl2Set6jointalteration];
RESULTS.Set6.Mean_Joint_Alteration_Factor=mean(RESULTS.Set6.Joint_Alteration_Factors);
RESULTS.Set6.STD_Joint_Alteration_Factor=std(RESULTS.Set6.Joint_Alteration_Factors);

sl1Set6jointroughness=cat(1,sl1.RESULTS.Set6.Joint_Roughness_Factors);
sl2Set6jointroughness=cat(1,sl2.RESULTS.Set6.Joint_Roughness_Factors);
RESULTS.Set6.Joint_Roughness_Factors=[sl1Set6jointroughness sl2Set6jointroughness];
RESULTS.Set6.Mean_Joint_Roughness_Factor=mean(RESULTS.Set6.Joint_Roughness_Factors);
RESULTS.Set6.STD_Joint_Roughness_Factor=std(RESULTS.Set6.Joint_Roughness_Factors);

sl1Set6jointsize=cat(1,sl1.RESULTS.Set6.Joint_Size_Factor);
sl2Set6jointsize=cat(1,sl2.RESULTS.Set6.Joint_Size_Factor);
RESULTS.Set6.Joint_Size_Factor=[sl1Set6jointsize sl2Set6jointsize];
RESULTS.Set6.Mean_Joint_Size_Factor=mean(RESULTS.Set6.Joint_Size_Factor);
RESULTS.Set6.STD_Joint_Size_Factor=std(RESULTS.Set6.Joint_Size_Factor);

sl1Set6Intersections=cat(1,sl1.RESULTS.Set6.Intersections);
sl2Set6Intersections=cat(1,sl2.RESULTS.Set6.Intersections);
RESULTS.Set6.Intersections=[sl1Set6Intersections sl2Set6Intersections];

sl1Set6terminations=cat(1,sl1.RESULTS.Set6.Terminations);
sl2Set6terminations=cat(1,sl2.RESULTS.Set6.Terminations);
RESULTS.Set6.Terminations=[sl1Set6terminations sl2Set6terminations];
RESULTS.Set6.Mean_Termination=mean(RESULTS.Set6.Terminations);

sl1Set6termination=cat(1,sl1.RESULTS.Set6.Termination_index);
sl2Set6termination=cat(1,sl2.RESULTS.Set6.Termination_index);
weigthedsl1Set6termination=(sl1Set6observations)/(sl1Set6observations+sl2Set6observations)*sl1Set6termination;
weigthedsl2Set6termination=(sl2Set6observations)/(sl1Set6observations+sl2Set6observations)*sl2Set6termination;
RESULTS.Set6.Termination_index=weigthedsl1Set6termination+weigthedsl2Set6termination;


sl1Set6DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set6.Normal_spacings);
sl2Set6DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set6.Normal_spacings);
RESULTS.Set6.Normal_spacings=[sl1Set6DiscontinuityNormalSpacings sl2Set6DiscontinuityNormalSpacings];
RESULTS.Set6.Mean_normal_spacing=mean(RESULTS.Set6.Normal_spacings);
RESULTS.Set6.Median_normal_spacing=median(RESULTS.Set6.Normal_spacings);
RESULTS.Set6.Mode_normal_spacing=mode(RESULTS.Set6.Normal_spacings);
RESULTS.Set6.STD_normal_spacing=std(RESULTS.Set6.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set6.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set6observations)/(sl1Set6observations+sl2Set6observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set6.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set6observations)/(sl1Set6observations+sl2Set6observations)*sl2MeanLaslettTL;
RESULTS.Set6.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set6.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set6observations)/(sl1Set6observations+sl2Set6observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set6.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set6observations)/(sl1Set6observations+sl2Set6observations)*sl2MeanPriest93TL;
RESULTS.Set6.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set6.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set6observations)/(sl1Set6observations+sl2Set6observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set6.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set6observations)/(sl1Set6observations+sl2Set6observations)*sl2MeanPriest81TL;
RESULTS.Set6.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end

%Set7
if isfield(sl1.RESULTS,'Set7')==1
setnumber=7;    
sl1Set7observations=cat(1,sl1.RESULTS.Set7.Number_of_joints_in_set);
sl2Set7observations=cat(1,sl2.RESULTS.Set7.Number_of_joints_in_set);
RESULTS.Set7.Number_of_joints_in_set=[sl1Set7observations sl2Set7observations];
RESULTS.Set7.Number_of_joints_in_set=sum(RESULTS.Set7.Number_of_joints_in_set);

sl1Set7dips=cat(1,sl1.RESULTS.Set7.Dips);
sl2Set7dips=cat(1,sl2.RESULTS.Set7.Dips);
RESULTS.Set7.Dips=[sl1Set7dips sl2Set7dips];

sl1Set7dipdirs=cat(1,sl1.RESULTS.Set7.Dipdirs);
sl2Set7dipdirs=cat(1,sl2.RESULTS.Set7.Dipdirs);
RESULTS.Set7.Dipdirs=[sl1Set7dipdirs sl2Set7dipdirs];

    %Calculates the centroids of the clusters
    dip=RESULTS.Set7.Dips;
    dipdir=RESULTS.Set7.Dipdirs;
    for n=1:length(dip);
        if dipdir(n)>180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)-180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)-180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        elseif dipdir(n)<=180
            X1b(n) = cosd(90-dip(n))*sind(dipdir(n)+180); %x-coordinate
            X2b(n) = cosd(90-dip(n))*cosd(dipdir(n)+180); %y-coordinate
            X3b(n) = [-sind(90-dip(n))]'; % z-coordinate negative because is on the lower hemisphere projection
        end
    end
    %Calculates products from direction cosines
        xx=sum(X1b.^2);
        xy=sum(X1b.*X2b);
        xz=sum(X1b.*X3b);
        yy=sum(X2b.^2);
        yz=sum(X2b.*X3b);
        zz=sum(X3b.^2);
        %Orientation matrix
        S(:,:,1)=[xx xy xz; xy yy yz;xz yz zz];
        %Orientation matrix normalized with number of observation in the cluster
        Sn=S./numel(dip);
        %B eigenvectors, K eigenvalues
        [Bi,Ki]=eig(Sn(:,:,1));
        B(:,:,1)=Bi;
        K(:,:,1)=Ki;
        %Eigenvector associated with the highest eigenvalue-> mean vector of the
        %group of N vectors
        Kmax(:,:,1)=max(K(:,:,1));
        n2=find(Kmax(:,:,1)==max(Kmax(:,:,1)));
        mean_vector(:,1)= B(:,n2,1);
        %Eigenvector's x,y and z coordinates 
        xeig=mean_vector(1,1);
        yeig=mean_vector(2,1);
        zeig=mean_vector(3,1);
        %Transforms the coordinates into dip/dipdir form
        %Pole vector's orientations in quarters
        if xeig>0 & yeig>0
            dipeigv=90-abs((asin(zeig)))*(180/pi);
            dipdireigv=abs((atan(xeig/yeig)))*180/pi;
        elseif xeig>0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180-abs((atan(xeig/yeig)))*180/pi;
        elseif xeig<0 & yeig<0
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=180+abs((atan(xeig/yeig)))*180/pi;
        else
            dipeigv=90-abs(asin(zeig))*180/pi;
            dipdireigv=360-abs((atan(xeig/yeig)))*180/pi;
        end
 
RESULTS.Set7.Dip=dipeigv;
RESULTS.Set7.Dipdir=dipdireigv;

  %Calculation of Standard deviations of the Joint orientations within the Clustersets
    stdangleclu=zeros(1,1);
    dipang=dip;
    dipdirang=dipdir;
    dipdireigvang=dipdireigv;
    dipeigvang=dipeigv;
    for n2=1:length(dip);
            if dipdir(n2)>180
                dipdirang(n2)=dipdir(n2)-180;
            else 
                dipdirang(n2)=dipdir(n2)+180;
            end

            if dipdireigv(1)>180
                dipdireigvang(1)=dipdireigv(1)-180;
            else 
                dipdireigvang(1)=dipdireigv(1)+180;
            end
            dipang(n2)=90-dip(n2);
            dipeigvang(1)=90-dipeigv(1);
            angleclu(n2)=acosd(abs((cosd(dipdirang(n2)-dipdireigvang(1)).*cosd(dipang(n2)).*cosd(dipeigvang(1)))+(sind(dipang(n2)).*sind(dipeigvang(1)))));
        
         
      end
   	angleclu_set=angleclu;
	stdangleclu=std(angleclu_set);
    RESULTS.Set7.STD_Orientation=stdangleclu;

sl1Set7TraceLengths=cat(1,sl1.RESULTS.Set7.Tracelengths);
sl2Set7TraceLengths=cat(1,sl2.RESULTS.Set7.Tracelengths);
RESULTS.Set7.Tracelengths=[sl1Set7TraceLengths sl2Set7TraceLengths];
RESULTS.Set7.Mean_Tracelength=mean(RESULTS.Set7.Tracelengths);
RESULTS.Set7.STD_Tracelength=std(RESULTS.Set7.Tracelengths);

sl1Set7Joint_Condition_Factor=cat(1,sl1.RESULTS.Set7.Mean_Joint_Condition_Factor);
sl2Set7Joint_Condition_Factor=cat(1,sl2.RESULTS.Set7.Mean_Joint_Condition_Factor);
weigthedsl1Set7condition=(sl1Set7observations)/(sl1Set7observations+sl2Set7observations)*sl1Set7Joint_Condition_Factor;
weigthedsl2Set7condition=(sl2Set7observations)/(sl1Set7observations+sl2Set7observations)*sl2Set7Joint_Condition_Factor;
RESULTS.Set7.Mean_Joint_Condition_Factor=weigthedsl1Set7condition+weigthedsl2Set7condition;


sl1Set7Phi_peak_jointsurface=cat(1,sl1.RESULTS.Set7.Peak_Friction_Angles);
sl2Set7Phi_peak_jointsurface=cat(1,sl2.RESULTS.Set7.Peak_Friction_Angles);
RESULTS.Set7.Peak_Friction_Angles=[sl1Set7Phi_peak_jointsurface sl2Set7Phi_peak_jointsurface];
RESULTS.Set7.Mean_Peak_Friction_Angle=mean(RESULTS.Set7.Peak_Friction_Angles);
RESULTS.Set7.STD_Peak_Friction_Angle=std(RESULTS.Set7.Peak_Friction_Angles);

sl1Set7jointalteration=cat(1,sl1.RESULTS.Set7.Joint_Alteration_Factors);
sl2Set7jointalteration=cat(1,sl2.RESULTS.Set7.Joint_Alteration_Factors);
RESULTS.Set7.Joint_Alteration_Factors=[sl1Set7jointalteration sl2Set7jointalteration];
RESULTS.Set7.Mean_Joint_Alteration_Factor=mean(RESULTS.Set7.Joint_Alteration_Factors);
RESULTS.Set7.STD_Joint_Alteration_Factor=std(RESULTS.Set7.Joint_Alteration_Factors);

sl1Set7jointroughness=cat(1,sl1.RESULTS.Set7.Joint_Roughness_Factors);
sl2Set7jointroughness=cat(1,sl2.RESULTS.Set7.Joint_Roughness_Factors);
RESULTS.Set7.Joint_Roughness_Factors=[sl1Set7jointroughness sl2Set7jointroughness];
RESULTS.Set7.Mean_Joint_Roughness_Factor=mean(RESULTS.Set7.Joint_Roughness_Factors);
RESULTS.Set7.STD_Joint_Roughness_Factor=std(RESULTS.Set7.Joint_Roughness_Factors);

sl1Set7jointsize=cat(1,sl1.RESULTS.Set7.Joint_Size_Factor);
sl2Set7jointsize=cat(1,sl2.RESULTS.Set7.Joint_Size_Factor);
RESULTS.Set7.Joint_Size_Factor=[sl1Set7jointsize sl2Set7jointsize];
RESULTS.Set7.Mean_Joint_Size_Factor=mean(RESULTS.Set7.Joint_Size_Factor);
RESULTS.Set7.STD_Joint_Size_Factor=std(RESULTS.Set7.Joint_Size_Factor);

sl1Set7Intersections=cat(1,sl1.RESULTS.Set7.Intersections);
sl2Set7Intersections=cat(1,sl2.RESULTS.Set7.Intersections);
RESULTS.Set7.Intersections=[sl1Set7Intersections sl2Set7Intersections];

sl1Set7terminations=cat(1,sl1.RESULTS.Set7.Terminations);
sl2Set7terminations=cat(1,sl2.RESULTS.Set7.Terminations);
RESULTS.Set7.Terminations=[sl1Set7terminations sl2Set7terminations];
RESULTS.Set7.Mean_Termination=mean(RESULTS.Set7.Terminations);

sl1Set7termination=cat(1,sl1.RESULTS.Set7.Termination_index);
sl2Set7termination=cat(1,sl2.RESULTS.Set7.Termination_index);
weigthedsl1Set7termination=(sl1Set7observations)/(sl1Set7observations+sl2Set7observations)*sl1Set7termination;
weigthedsl2Set7termination=(sl2Set7observations)/(sl1Set7observations+sl2Set7observations)*sl2Set7termination;
RESULTS.Set7.Termination_index=weigthedsl1Set7termination+weigthedsl2Set7termination;


sl1Set7DiscontinuityNormalSpacings=cat(1,sl1.RESULTS.Set7.Normal_spacings);
sl2Set7DiscontinuityNormalSpacings=cat(1,sl2.RESULTS.Set7.Normal_spacings);
RESULTS.Set7.Normal_spacings=[sl1Set7DiscontinuityNormalSpacings sl2Set7DiscontinuityNormalSpacings];
RESULTS.Set7.Mean_normal_spacing=mean(RESULTS.Set7.Normal_spacings);
RESULTS.Set7.Median_normal_spacing=median(RESULTS.Set7.Normal_spacings);
RESULTS.Set7.Mode_normal_spacing=mode(RESULTS.Set7.Normal_spacings);
RESULTS.Set7.STD_normal_spacing=std(RESULTS.Set7.Normal_spacings);

sl1MeanLaslettTL=cat(1,sl1.RESULTS.Set7.Mean_tracelength_Laslett1982);
weightedsl1MeanLaslettTL=(sl1Set7observations)/(sl1Set7observations+sl2Set7observations)*sl1MeanLaslettTL;
sl2MeanLaslettTL=cat(1,sl2.RESULTS.Set7.Mean_tracelength_Laslett1982);
weightedsl2MeanLaslettTL=(sl2Set7observations)/(sl1Set7observations+sl2Set7observations)*sl2MeanLaslettTL;
RESULTS.Set7.Mean_tracelength_Laslett1982=weightedsl1MeanLaslettTL+weightedsl2MeanLaslettTL;

sl1MeanPriest93TL=cat(1,sl2.RESULTS.Set7.Mean_tracelength_Priest_1993);
weightedsl1MeanPriest93TL=(sl1Set7observations)/(sl1Set7observations+sl2Set7observations)*sl1MeanPriest93TL;
sl2MeanPriest93TL=cat(1,sl2.RESULTS.Set7.Mean_tracelength_Priest_1993);
weightedsl2MeanPriest93TL=(sl2Set7observations)/(sl1Set7observations+sl2Set7observations)*sl2MeanPriest93TL;
RESULTS.Set7.Mean_tracelength_Priest_1993=weightedsl1MeanPriest93TL+weightedsl2MeanPriest93TL;

sl1MeanPriest81TL=cat(1,sl2.RESULTS.Set7.Mean_tracelength_PriestHudson_1981);
weightedsl1MeanPriest81TL=(sl1Set7observations)/(sl1Set7observations+sl2Set7observations)*sl1MeanPriest81TL;
sl2MeanPriest81TL=cat(1,sl2.RESULTS.Set7.Mean_tracelength_PriestHudson_1981);
weightedsl2MeanPriest81TL=(sl2Set7observations)/(sl1Set7observations+sl2Set7observations)*sl2MeanPriest81TL;
RESULTS.Set7.Mean_tracelength_PriestHudson_1981=weightedsl1MeanPriest81TL+weightedsl2MeanPriest81TL;
end


%Block volume
sl1MeanPalmstroemvol=cat(1,sl1.RESULTS.Block_volume_calculations.Means.Palmstroem_volume);
weightedsl1MeanPalmstroemvol=(sl1length/(sl1length+sl2length)*sl1MeanPalmstroemvol);
sl2MeanPalmstroemvol=cat(1,sl2.RESULTS.Block_volume_calculations.Means.Palmstroem_volume);
weightedsl1MeanPalmstroemvol=(sl2length/(sl1length+sl2length))*sl2MeanPalmstroemvol;
RESULTS.Block_volume_calculations.Means.Palmstroem_volume=weightedsl1MeanPalmstroemvol+weightedsl1MeanPalmstroemvol;

sl1MeanJvmvol=cat(1,sl1.RESULTS.Block_volume_calculations.Means.Jv);
weightedsl1MeanJvvol=(sl1length/(sl1length+sl2length)*sl1MeanJvmvol);
sl2MeanJvvol=cat(1,sl2.RESULTS.Block_volume_calculations.Means.Jv);
weightedsl1MeanJvvol=(sl2length/(sl1length+sl2length))*sl2MeanJvvol;
RESULTS.Block_volume_calculations.Means.Jv=weightedsl1MeanJvvol+weightedsl1MeanJvvol;

sl1isbdfractal=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Fractal);
weightedsl1isbdfractal=(sl1length/(sl1length+sl2length).*sl1isbdfractal);
sl2isbdfractal=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Fractal);
weightedsl2isbdfractal=(sl2length/(sl1length+sl2length).*sl2isbdfractal);
RESULTS.Block_volume_calculations.ISBDs.Fractal=weightedsl1isbdfractal+weightedsl2isbdfractal;

sl1isbdlognormal=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Lognormal);
weightedsl1isbdlognormal=(sl1length/(sl1length+sl2length).*sl1isbdlognormal);
sl2isbdlognormal=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Lognormal);
weightedsl2isbdlognormal=(sl2length/(sl1length+sl2length).*sl2isbdlognormal);
RESULTS.Block_volume_calculations.ISBDs.Lognormal=weightedsl1isbdlognormal+weightedsl2isbdlognormal;

sl1isbdtwoneg=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.TwoNegExp_OneUniform);
weightedsl1isbdtwoneg=(sl1length/(sl1length+sl2length).*sl1isbdtwoneg);
sl2isbdtwoneg=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.TwoNegExp_OneUniform);
weightedsl2isbdtwoneg=(sl2length/(sl1length+sl2length).*sl2isbdtwoneg);
RESULTS.Block_volume_calculations.ISBDs.TwoNegExp_OneUniform=weightedsl1isbdtwoneg+weightedsl2isbdtwoneg;

sl1isbduniform=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Uniform);
weightedsl1isbduniform=(sl1length/(sl1length+sl2length).*sl1isbduniform);
sl2isbduniform=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Uniform);
weightedsl2isbduniform=(sl2length/(sl1length+sl2length).*sl2isbduniform);
RESULTS.Block_volume_calculations.ISBDs.Uniform=weightedsl1isbduniform+weightedsl2isbduniform;

sl1isbdnegexp=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Negativeexponential);
weightedsl1isbdnegexp=(sl1length/(sl1length+sl2length).*sl1isbdnegexp);
sl2isbdnegexp=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Negativeexponential);
weightedsl2isbdnegexp=(sl2length/(sl1length+sl2length).*sl2isbdnegexp);
RESULTS.Block_volume_calculations.ISBDs.Negativeexponential=weightedsl1isbdnegexp+weightedsl2isbdnegexp;

sl1isbdpalm=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Palmstroem_volume);
weightedsl1isbdpalm=(sl1length/(sl1length+sl2length).*sl1isbdpalm);
sl2isbdpalm=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Palmstroem_volume);
weightedsl2isbdpalm=(sl2length/(sl1length+sl2length).*sl2isbdpalm);
RESULTS.Block_volume_calculations.ISBDs.Palmstroem_volume=weightedsl1isbdpalm+weightedsl2isbdpalm;

sl1isbdlatham=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Latham_et_al_drillingcores);
weightedsl1isbdlatham=(sl1length/(sl1length+sl2length).*sl1isbdlatham);
sl2isbdlatham=cat(1,sl1.RESULTS.Block_volume_calculations.ISBDs.Latham_et_al_drillingcores);
weightedsl2isbdlatham=(sl2length/(sl1length+sl2length).*sl2isbdlatham);
RESULTS.Block_volume_calculations.ISBDs.Latham_et_al_drillingcores=weightedsl1isbdlatham+weightedsl2isbdlatham;


%Rock mass characteristics
sl1RQDs=cat(1,sl1.RESULTS.Rock_characteristics.RQDs);
sl2RQDs=cat(1,sl2.RESULTS.Rock_characteristics.RQDs);
RESULTS.Rock_characteristics.RQDs=[sl1RQDs sl2RQDs];
RESULTS.Rock_characteristics.Mean_RQD=mean(RESULTS.Rock_characteristics.RQDs);
RESULTS.Rock_characteristics.Min_RQD=min(RESULTS.Rock_characteristics.RQDs);
RESULTS.Rock_characteristics.Max_RQD=max(RESULTS.Rock_characteristics.RQDs);

sl1wJd=cat(1,sl1.RESULTS.Rock_characteristics.wJd);
weightedsl1wJd=(sl1length/(sl1length+sl2length)*sl1wJd);
sl2wJd=cat(1,sl2.RESULTS.Rock_characteristics.wJd);
weightedsl2wJd=(sl2length/(sl1length+sl2length))*sl2wJd;
RESULTS.Rock_characteristics.wJd=weightedsl1wJd+weightedsl2wJd;

sl1Jv=cat(1,sl1.RESULTS.Rock_characteristics.Jv);
weightedsl1Jv=(sl1length/(sl1length+sl2length)*sl1Jv);
sl2Jv=cat(1,sl2.RESULTS.Rock_characteristics.Jv);
weightedsl2Jv=(sl2length/(sl1length+sl2length))*sl2Jv;
RESULTS.Rock_characteristics.Jv=weightedsl1Jv+weightedsl2Jv;

%rmicheck=isfield(RESULTS.Rock_characteristics,'RMi');
%if rmicheck==1
%	sl1rmi=cat(1,sl1.RESULTS.Rock_characteristics.RMi);
%	weightedsl1rmi=(sl1length/(sl1length+sl2length))*sl1rmi;
%	sl2rmi=cat(1,sl2.RESULTS.Rock_characteristics.RMi);
%	weightedsl2rmi=(sl2length/(sl1length+sl2length))*sl2rmi;
%	RESULTS.Rock_characteristics.RMi=weightedsl1rmi+weightedsl2rmi;
%end

gsicheck=isfield(RESULTS.Rock_characteristics,'Mean_GSI');
if gsicheck==1
	sl1gsimean=cat(1,sl1.RESULTS.Rock_characteristics.Mean_GSI);
	weightedsl1gsimean=(sl1length/(sl1length+sl2length))*sl1gsimean;
	sl2gsimean=cat(1,sl2.RESULTS.Rock_characteristics.Mean_GSI);
	weightedsl2gsimean=(sl2length/(sl1length+sl2length))*sl2gsimean;
	RESULTS.Rock_characteristics.Mean_GSI=weightedsl1gsimean+weightedsl2gsimean;
    
    sl1gsicai=cat(1,sl1.RESULTS.Rock_characteristics.GSI_CaiandKaiser);
    weightedsl1gsicai=(sl1length/(sl1length+sl2length))*sl1gsicai;
    sl2gsicai=cat(1,sl2.RESULTS.Rock_characteristics.GSI_CaiandKaiser);
    weightedsl2gsicai=(sl2length/(sl1length+sl2length))*sl2gsicai;
    RESULTS.Rock_characteristics.GSI_CaiandKaiser=weightedsl1gsicai+weightedsl2gsicai;
    
    
    sl1gsirusso=cat(1,sl1.RESULTS.Rock_characteristics.GSI_Russo);
    weightedsl1gsirusso=(sl1length/(sl1length+sl2length))*sl1gsirusso;
    sl2gsirusso=cat(1,sl2.RESULTS.Rock_characteristics.GSI_Russo);
    weightedsl2gsirusso=(sl2length/(sl1length+sl2length))*sl2gsirusso;
    RESULTS.Rock_characteristics.GSI_Russo=weightedsl1gsirusso+weightedsl2gsirusso;    
    
    sl1gsibarton=cat(1,sl1.RESULTS.Rock_characteristics.GSI_Barton);
    weightedsl1gsibarton=(sl1length/(sl1length+sl2length))*sl1gsibarton;
    sl2gsibarton=cat(1,sl2.RESULTS.Rock_characteristics.GSI_Barton);
    weightedsl2gsibarton=(sl2length/(sl1length+sl2length))*sl2gsibarton;
    RESULTS.Rock_characteristics.GSI_Barton=weightedsl1gsibarton+weightedsl2gsibarton;      
end

%qapostrophcheck=isfield(RESULTS.Rock_characteristics,'Q_Apostroph');
%if qapostrophcheck==1
%	sl1qapo=cat(1,sl1.RESULTS.Rock_characteristics.Q_Apostroph);
%	weightedsl1qapo=(sl1length/(sl1length+sl2length))*sl1qapo;
%	sl2qapo=cat(1,sl2.RESULTS.Rock_characteristics.Q_Apostroph);
%	weightedsl2qapo=(sl2length/(sl1length+sl2length))*sl2qapo;
%	RESULTS.Rock_characteristics.Q_Apostroph=weightedsl1qapo+weightedsl2qapo;
%end

%qcheck=isfield(RESULTS.Rock_characteristics,'Q');
%if qcheck==1
%	sl1q=cat(1,sl1.RESULTS.Rock_characteristics.Q);
%	weightedsl1q=(sl1length/(sl1length+sl2length))*sl1q;
%	sl2q=cat(1,sl2.RESULTS.Rock_characteristics.Q);
%	weightedsl2q=(sl2length/(sl1length+sl2length))*sl2q;
%	RESULTS.Rock_characteristics.Q=weightedsl1q+weightedsl2q;
%end


%Creates new *.mat file
uisave('RESULTS','Results');

idx3=setnumber;


 set (handles.pushbutton_statistics,'Enable','on');
 set (handles.pushbutton_vol,'Enable','on');
 set (handles.pushbutton_general,'Enable','off');
 set (handles.pushbutton_clustering,'Enable','off');
 set (handles.pushbutton_rmc,'Enable','on');




% --- Executes on button press in pushbutton_rmc.
function pushbutton_rmc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_rmc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS 

set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');
set(handles.panel_rmc,'Visible','on');
set(handles.edit_rmc_rqd, 'String', RESULTS.Rock_characteristics.Mean_RQD); 

set (handles.edit_rmc_ja,'string',RESULTS.Joint_characteristics.Mean_Joint_Alteration_Factor);
set (handles.edit_rmc_jr,'string',RESULTS.Joint_characteristics.Mean_Joint_Roughness_Factor);
set (handles.edit_rmc_jl,'string',RESULTS.Joint_characteristics.Mean_Joint_Size_Factor);
set (handles.edit_rmc_jc,'string',RESULTS.Joint_characteristics.Mean_Joint_Condition_Factor);

set (handles.edit_rmc_vb,'String',RESULTS.Block_volume_calculations.Means.Palmstroem_volume);
 



% --- Executes on button press in push_rmc_rmi.
function push_rmc_rmi_Callback(hObject, eventdata, handles)
% hObject    handle to push_rmc_rmi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS 
UCSrmi=str2double(get(handles.edit_rmc_rmr_ucs,'String'));
Volrmi=str2double(get(handles.edit_rmc_vb,'String'));

jcrmi=str2double(get(handles.edit_rmc_jc,'String'));

%Following RMi calculation is based on Palmström (1996) and other
%publications by Palmström
RMi=0.2*UCSrmi*sqrt(jcrmi)*Volrmi^(0.37*jcrmi^-0.2);
JP=0.2*sqrt(jcrmi)*Volrmi^(0.37*jcrmi^-0.2);
if RMi<0.1
    rmiclass='Extremely low';;
elseif RMi>=0.1 & RMi<0.4
    rmiclass='Very low';
elseif RMi>=0.4 & RMi<1
    rmiclass='Low';
elseif RMi>=1 & RMi<10
    rmiclass='Moderate';
elseif RMi>=10 & RMi<40
    rmiclass='High';
else %RMi>=40
    rmiclass='Very high';
end
set(handles.text_rmc_rmi, 'String', RMi);
set(handles.text_class_rmi, 'String', rmiclass);
RESULTS.Rock_characteristics.RMi=RMi;



% --- Executes on button press in push_rmc_gsi.
function push_rmc_gsi_Callback(hObject, eventdata, handles)
% hObject    handle to push_rmc_gsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS

%Reads the parameters for GSI estimations:
VolGSI=str2double(get(handles.edit_rmc_vb,'String'));
VolGSICai=VolGSI*1000000;
jnGSI=str2double(get(handles.edit_rmc_jn,'String'));
RQDGSI=str2double(get(handles.edit_rmc_rqd,'String'));

%Following GSI estimation is made according to Cai & Kaiser 2006
jccaikaiserGSI=str2double(get(handles.edit_rmc_jc,'String'));
GSIcaikaiser=round((26.5+8.79*log(jccaikaiserGSI)+0.9*log(VolGSICai))/(1+0.0151*log(jccaikaiserGSI)-0.0253*log(VolGSICai)));
set(handles.text_rmc_gsicaikaiser, 'String', GSIcaikaiser);
RESULTS.Rock_characteristics.GSI_CaiandKaiser=GSIcaikaiser;

%Following GSI estimation is made according to Hoek et al. 2013
%jchoek_q=str2double(get(handles.edit_rmc_jc,'String'));
%GSIhoek_q=round(((52*jchoek_q)/(1+jchoek_q))+0.5*RQDGSI);
%set(handles.text_rmc_gsi_hoekq, 'String', GSIhoek_q);
%RESULTS.Rock_characteristics.GSI_HoeketalQ=GSIhoek_q;

%Following GSI estimation is made according to Russo 2009
jcpalmgsi=str2double(get(handles.edit_rmc_jc,'String'));
JP=0.2*sqrt(jcpalmgsi)*VolGSI^(0.37*jcpalmgsi^-0.2);
GSIrusso=round(153-(165/(1+((JP/0.19)^0.44))));
set(handles.text_rmc_gsi_russo, 'String', GSIrusso);
RESULTS.Rock_characteristics.GSI_Russo=GSIrusso;

%Following GSI estimation is made according to Barton 1995
jrbarton=str2double(get(handles.edit_rmc_jr,'String'));
jabarton=str2double(get(handles.edit_rmc_ja,'String'));
jcbarton=jrbarton/jabarton;
GSIbarton=round(15*log10((RQDGSI/jnGSI)*(jcbarton))+50);
set(handles.text_rmc_gsi_bartonq, 'String', GSIbarton);
RESULTS.Rock_characteristics.GSI_Barton=GSIbarton;

meanGSI=round((GSIcaikaiser+GSIrusso+GSIbarton)/3);
set(handles.text_meangsi, 'String', meanGSI);
RESULTS.Rock_characteristics.Mean_GSI=meanGSI;



% --- Executes on button press in push_rmc_qapo.
function push_rmc_qapo_Callback(hObject, eventdata, handles)
% hObject    handle to push_rmc_qapo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  RESULTS

jcGSI=str2double(get(handles.edit_rmc_jc,'String'));
jnGSI=str2double(get(handles.edit_rmc_jn,'String'));
RQDGSI=str2double(get(handles.edit_rmc_rqd,'String'));
Qapostroph=round((RQDGSI/jnGSI)*(jcGSI));
set(handles.text_rmc_qapo, 'String', Qapostroph);
RESULTS.Rock_characteristics.Q_Apostroph=Qapostroph;

% --- Executes on button press in push_rmc_q.
function push_rmc_q_Callback(hObject, eventdata, handles)
% hObject    handle to push_rmc_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  RESULTS
Qapostroph=RESULTS.Rock_characteristics.Q_Apostroph;

%Calculates Q
popupqwaterfactor=get(handles.popup_rmc_water,'Value');
if popupqwaterfactor==2 %dry
    q_jw=1;
end
if popupqwaterfactor==3 %dry
    q_jw=1;
end
if popupqwaterfactor==4 %wet
    q_jw=0.66;
end
if popupqwaterfactor==5 %dripping
    q_jw=0.5;
end
if popupqwaterfactor==6 %gushing
    q_jw=0.3;
end
if popupqwaterfactor==7 %flowing
    q_jw=0.15;
end
if popupqwaterfactor==8 %heavily flowing
    q_jw=0.08;
end

popupqstressfactor=get(handles.popup_rmc_q_srf,'Value');
if popupqstressfactor==2 %low stress
    q_srf=2.5;
end
if popupqstressfactor==3 %medium stress
    q_srf=1;
end
if popupqstressfactor==4 %high stress
    q_srf=0.67;
end
if popupqstressfactor==5 %moderate slabbing
    q_srf=25;
end
if popupqstressfactor==6 %slabbing and rock burst
    q_srf=100;
end
if popupqstressfactor==7 %heavy rock burst
    q_srf=300;
end
if popupqstressfactor==8 %mild squeezing
    q_srf=10;
end
if popupqstressfactor==9 %heavy squeezing
    q_srf=20;
end



%Following ratings and classifications are based on those of the computer
%spreadsheet http://rockmass.net/files/Q-RMR-RMi_v3.xls by RockMass as. (10.04.2016)
Q=Qapostroph*q_jw/q_srf;
set(handles.text_rmc_q, 'String', Q);
if Q<0.01
    qclass='Exceptionally poor';
elseif Q>=0.01 & Q<0.1
    qclass='Extremely poor';
elseif Q>=0.1 & Q<1
    qclass='Very poor';
elseif Q>=1 & Q<4
    qclass='Poor';
elseif Q>=4 & Q<10
    qclass='Fair';
elseif Q>=10 & Q<40
    qclass='Good';
elseif Q>=40 & Q<100
    qclass='Very good';
elseif Q>=100 & Q<400
    qclass='Extremely good';
else %Q>=400
    qclass='Except, good';
end
set(handles.text_class_q, 'String', qclass);
RESULTS.Rock_characteristics.Q=Q;



% --- Executes on button press in pushbutton_saveresults.
function pushbutton_saveresults_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global RESULTS
uisave('RESULTS','Results');



% --- Executes on selection change in popupDFN_number.
function popupDFN_number_Callback(hObject, eventdata, handles)
% hObject    handle to popupDFN_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupDFN_number contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupDFN_number
global numbersets


switch get(handles.popupDFN_number,'Value')  
    case 2
        numbersets=1;
        set(handles.editDFN_dipdirA,'Visible','on'); 
        set(handles.editDFN_dipdirB,'Visible','off'); 
        set(handles.editDFN_dipdirC,'Visible','off'); 
        set(handles.editDFN_dipdirD,'Visible','off'); 
        set(handles.editDFN_dipdirE,'Visible','off'); 
        set(handles.editDFN_dipA,'Visible','on'); 
        set(handles.editDFN_dipB,'Visible','off'); 
        set(handles.editDFN_dipC,'Visible','off'); 
        set(handles.editDFN_dipD,'Visible','off'); 
        set(handles.editDFN_dipE,'Visible','off'); 
        set(handles.editDFN_spaceA,'Visible','on'); 
        set(handles.editDFN_spaceB,'Visible','off'); 
        set(handles.editDFN_spaceC,'Visible','off'); 
        set(handles.editDFN_spaceD,'Visible','off'); 
        set(handles.editDFN_spaceE,'Visible','off'); 
        set(handles.editDFN_diaA,'Visible','on'); 
        set(handles.editDFN_diaB,'Visible','off'); 
        set(handles.editDFN_diaC,'Visible','off');
        set(handles.editDFN_diaD,'Visible','off'); 
        set(handles.editDFN_diaE,'Visible','off'); 
        set(handles.editDFN_stdA,'Visible','on');
        set(handles.editDFN_stdB,'Visible','off'); 
        set(handles.editDFN_stdC,'Visible','off'); 
        set(handles.editDFN_stdD,'Visible','off'); 
        set(handles.editDFN_stdE,'Visible','off'); 
        
    case 3
        set(handles.editDFN_dipdirA,'Visible','on'); 
        set(handles.editDFN_dipdirB,'Visible','on'); 
        set(handles.editDFN_dipdirC,'Visible','off'); 
        set(handles.editDFN_dipdirD,'Visible','off'); 
        set(handles.editDFN_dipdirE,'Visible','off'); 
        set(handles.editDFN_dipA,'Visible','on'); 
        set(handles.editDFN_dipB,'Visible','on'); 
        set(handles.editDFN_dipC,'Visible','off'); 
        set(handles.editDFN_dipD,'Visible','off'); 
        set(handles.editDFN_dipE,'Visible','off'); 
        set(handles.editDFN_spaceA,'Visible','on'); 
        set(handles.editDFN_spaceB,'Visible','on'); 
        set(handles.editDFN_spaceC,'Visible','off'); 
        set(handles.editDFN_spaceD,'Visible','off'); 
        set(handles.editDFN_spaceE,'Visible','off'); 
        set(handles.editDFN_diaA,'Visible','on'); 
        set(handles.editDFN_diaB,'Visible','on'); 
        set(handles.editDFN_diaC,'Visible','off');
        set(handles.editDFN_diaD,'Visible','off'); 
        set(handles.editDFN_diaE,'Visible','off'); 
        set(handles.editDFN_stdA,'Visible','on');
        set(handles.editDFN_stdB,'Visible','on'); 
        set(handles.editDFN_stdC,'Visible','off'); 
        set(handles.editDFN_stdD,'Visible','off'); 
        set(handles.editDFN_stdE,'Visible','off'); 
    case 4
        set(handles.editDFN_dipdirA,'Visible','on'); 
        set(handles.editDFN_dipdirB,'Visible','on'); 
        set(handles.editDFN_dipdirC,'Visible','on'); 
        set(handles.editDFN_dipdirD,'Visible','off'); 
        set(handles.editDFN_dipdirE,'Visible','off'); 
        set(handles.editDFN_dipA,'Visible','on'); 
        set(handles.editDFN_dipB,'Visible','on'); 
        set(handles.editDFN_dipC,'Visible','on'); 
        set(handles.editDFN_dipD,'Visible','off'); 
        set(handles.editDFN_dipE,'Visible','off'); 
        set(handles.editDFN_spaceA,'Visible','on'); 
        set(handles.editDFN_spaceB,'Visible','on'); 
        set(handles.editDFN_spaceC,'Visible','on'); 
        set(handles.editDFN_spaceD,'Visible','off'); 
        set(handles.editDFN_spaceE,'Visible','off'); 
        set(handles.editDFN_diaA,'Visible','on'); 
        set(handles.editDFN_diaB,'Visible','on'); 
        set(handles.editDFN_diaC,'Visible','on');
        set(handles.editDFN_diaD,'Visible','off'); 
        set(handles.editDFN_diaE,'Visible','off'); 
        set(handles.editDFN_stdA,'Visible','on');
        set(handles.editDFN_stdB,'Visible','on'); 
        set(handles.editDFN_stdC,'Visible','on'); 
        set(handles.editDFN_stdD,'Visible','off'); 
        set(handles.editDFN_stdE,'Visible','off'); 
    case 5
        set(handles.editDFN_dipdirA,'Visible','on'); 
        set(handles.editDFN_dipdirB,'Visible','on'); 
        set(handles.editDFN_dipdirC,'Visible','on'); 
        set(handles.editDFN_dipdirD,'Visible','on'); 
        set(handles.editDFN_dipdirE,'Visible','off'); 
        set(handles.editDFN_dipA,'Visible','on'); 
        set(handles.editDFN_dipB,'Visible','on'); 
        set(handles.editDFN_dipC,'Visible','on'); 
        set(handles.editDFN_dipD,'Visible','on'); 
        set(handles.editDFN_dipE,'Visible','off'); 
        set(handles.editDFN_spaceA,'Visible','on'); 
        set(handles.editDFN_spaceB,'Visible','on'); 
        set(handles.editDFN_spaceC,'Visible','on'); 
        set(handles.editDFN_spaceD,'Visible','on'); 
        set(handles.editDFN_spaceE,'Visible','off'); 
        set(handles.editDFN_diaA,'Visible','on'); 
        set(handles.editDFN_diaB,'Visible','on'); 
        set(handles.editDFN_diaC,'Visible','on');
        set(handles.editDFN_diaD,'Visible','on'); 
        set(handles.editDFN_diaE,'Visible','off'); 
        set(handles.editDFN_stdA,'Visible','on');
        set(handles.editDFN_stdB,'Visible','on'); 
        set(handles.editDFN_stdC,'Visible','on'); 
        set(handles.editDFN_stdD,'Visible','on'); 
        set(handles.editDFN_stdE,'Visible','off'); 
    case 6
        set(handles.editDFN_dipdirA,'Visible','on'); 
        set(handles.editDFN_dipdirB,'Visible','on'); 
        set(handles.editDFN_dipdirC,'Visible','on'); 
        set(handles.editDFN_dipdirD,'Visible','on'); 
        set(handles.editDFN_dipdirE,'Visible','on'); 
        set(handles.editDFN_dipA,'Visible','on'); 
        set(handles.editDFN_dipB,'Visible','on'); 
        set(handles.editDFN_dipC,'Visible','on'); 
        set(handles.editDFN_dipD,'Visible','on'); 
        set(handles.editDFN_dipE,'Visible','on'); 
        set(handles.editDFN_spaceA,'Visible','on'); 
        set(handles.editDFN_spaceB,'Visible','on'); 
        set(handles.editDFN_spaceC,'Visible','on'); 
        set(handles.editDFN_spaceD,'Visible','on'); 
        set(handles.editDFN_spaceE,'Visible','on'); 
        set(handles.editDFN_diaA,'Visible','on'); 
        set(handles.editDFN_diaB,'Visible','on'); 
        set(handles.editDFN_diaC,'Visible','on');
        set(handles.editDFN_diaD,'Visible','on'); 
        set(handles.editDFN_diaE,'Visible','on'); 
        set(handles.editDFN_stdA,'Visible','on');
        set(handles.editDFN_stdB,'Visible','on'); 
        set(handles.editDFN_stdC,'Visible','on'); 
        set(handles.editDFN_stdD,'Visible','on'); 
        set(handles.editDFN_stdE,'Visible','on');
end



% --- Executes on button press in pushbuttonDFN_calculate.
function pushbuttonDFN_calculate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDFN_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numbersets 
        set(handles.textDFN_calc,'Visible','on');
set(handles.text296,'Visible','off');
set(handles.text297,'Visible','off');
dimrast = str2double(get (handles.editDFN_raster,'String'));
dimensions = str2double(get (handles.editDFN_size,'String'));
xlo=0;
xup=dimensions;
ylo=0;
yup=dimensions;
zlo=0;
zup=dimensions;

switch get(handles.popupDFN_number,'Value')  
    case 2
        numbersets=1;
    case 3
        numbersets=2;
         case 4
        numbersets=3;
    case 5
        numbersets=4;
     case 6
        numbersets=5;  
end
        
        
if numbersets==1 || numbersets==2 || numbersets==3 || numbersets==4 || numbersets==5        
DFN.SetA.dipdirmean=str2double(get (handles.editDFN_dipdirA,'String'));
DFN.SetA.dipmean=str2double(get (handles.editDFN_dipA,'String'));
DFN.SetA.orientstandard=str2double(get (handles.editDFN_stdA,'String'));
DFN.SetA.rmean=(str2double(get (handles.editDFN_diaA,'String')))/2;
DFN.SetA.onedimfreq=1/(str2double(get (handles.editDFN_spaceA,'String')));

%Calculation of number of joints
mean_radius=DFN.SetA.rmean;
onedimfreq=DFN.SetA.onedimfreq;       
% tolerance value within which the simulated 1D frequency should fall.
tolval = 1e-3;
% no. of imaginary scanlines for which 1D frequency is calculated. the larger the number
% the longer the simulation will take.
nsims = 100;
% two ends of search space for number of circles. 
maxguess=(onedimfreq*2.5)/((mean_radius^2)*pi)*xup*yup*zup;  %usually approx. 5 times larger than the final number of circles
ncircles = [1 maxguess];
freq = 0;

while abs(onedimfreq - freq) > tolval && ncircles(1)~=ncircles(2)
    nc = round(mean(ncircles));
    pproc = xup*rand(nc, 1) + 1i*zup*rand(nc, 1);
    circle_r = exprnd(mean_radius, nc, 1);
    
    lines_loc = xup*rand(1,nsims) + 1i*zup*rand(1, nsims);
    dist_mat = abs(repmat(pproc, 1, nsims) - repmat(lines_loc, nc, 1));
    
    intersect_circles = sum((dist_mat-repmat(circle_r, 1, nsims))<0);
    freq = mean(intersect_circles/zup);
    
    if freq>onedimfreq
        ncircles(end) = nc;
    else
        ncircles(1) = nc;
    end
end
DFN.SetA.njoint=ncircles(1);


% 3D joint point generation (Poisson process)
DFN.SetA.xcenter=xlo+(xup-xlo)*rand(DFN.SetA.njoint,1);
DFN.SetA.ycenter=ylo+(yup-ylo)*rand(DFN.SetA.njoint,1);
DFN.SetA.zcenter=zlo+(zup-ylo)*rand(DFN.SetA.njoint,1);

%Normal vectors
DFN.SetA.meannormalx=sind(DFN.SetA.dipdirmean).*cosd(90-DFN.SetA.dipmean);
DFN.SetA.meannormaly=cosd(DFN.SetA.dipdirmean).*cosd(90-DFN.SetA.dipmean);
DFN.SetA.meannormalz=sind(90-DFN.SetA.dipmean);
DFN.SetA.normalx=normrnd(DFN.SetA.meannormalx,(cosd(90-DFN.SetA.orientstandard)),DFN.SetA.njoint,1);
DFN.SetA.normaly=normrnd(DFN.SetA.meannormaly,(cosd(90-DFN.SetA.orientstandard)),DFN.SetA.njoint,1);
DFN.SetA.dip=normrnd(DFN.SetA.dipmean,DFN.SetA.orientstandard,DFN.SetA.njoint,1);
DFN.SetA.normalz=cosd(DFN.SetA.dip);
radii=exprnd(DFN.SetA.rmean,DFN.SetA.njoint,1);
radii_t=(randsample(radii(radii>0.1),DFN.SetA.njoint,1));  %Truncation of joints smaller than 20 cm in diameter
meanRT=mean(radii_t);
meanR=DFN.SetA.rmean;

while meanRT > (meanR*1.025)
radii=radii.*0.99;
    radii_t=(randsample(radii(radii>0.1),DFN.SetA.njoint,1));
meanRT=mean(radii_t);    
end
DFN.SetA.radius=radii_t;





clear radii radii_t meanR meanRT
u=numbersets;
for z=1:u
for p=1:length(DFN.SetA.normalx)
    if DFN.SetA.normalx(p,:)>1
        DFN.SetA.normalx(p,:)=1-(DFN.SetA.normalx(p,:)-1);
    end
    if DFN.SetA.normalx(p,:)<-1
        DFN.SetA.normalx(p,:)=-1+(-1-DFN.SetA.normalx(p,:));
    end
    if DFN.SetA.normaly(p,:)>1
        DFN.SetA.normaly(p,:)=1-(DFN.SetA.normaly(p,:)-1);
    end
    if DFN.SetA.normaly(p,:)<-1
        DFN.SetA.normaly(p,:)=-1+(-1-DFN.SetA.normaly(p,:));
    end
        if DFN.SetA.normalz(p,:)>1
        DFN.SetA.normalz(p,:)=-1+(DFN.SetA.normalz(p,:)-1);
       
        end
    if DFN.SetA.normalz(p,:)<-1
        DFN.SetA.normalz(p,:)=1-(DFN.SetA.normalz(p,:)+1);

    end
	
	end
end


z=1;
for z=1:u

DFN.SetA.TEMP.lowerx=DFN.SetA.normalx.*sign(DFN.SetA.normalz);
DFN.SetA.TEMP.lowery=DFN.SetA.normaly.*sign(DFN.SetA.normalz);
DFN.SetA.TEMP.lowerz=DFN.SetA.normalz.*sign(DFN.SetA.normalz);



if DFN.SetA.TEMP.lowerx==0 & DFN.SetA.TEMP.lowery==0
    DFN.SetA.TEMP.poleazimuth=0
elseif DFN.SetA.TEMP.lowerx<0 & DFN.SetA.TEMP.lowery>=0
    DFN.SetA.TEMP.poleazimuth=450-(atan2d(DFN.SetA.TEMP.lowery,DFN.SetA.TEMP.lowerx));
else
    DFN.SetA.TEMP.poleazimuth=90-(atan2d(DFN.SetA.TEMP.lowery,DFN.SetA.TEMP.lowerx));
end

DFN.SetA.TEMP.poleplunge=90-acosd((DFN.SetA.TEMP.lowerz));
DFN.SetA.dipdir=DFN.SetA.TEMP.poleazimuth;


% Simulation of the Radius distribution
% uniform distribution
DFN.SetA.dip=90-DFN.SetA.TEMP.poleplunge;
end
% convert dip dir/dip angle of joint in normal vectors
u=numbersets;
for z=1:u
for p=1:length(DFN.SetA.dip)
    if DFN.SetA.dipdir(p,:)>=360
        DFN.SetA.dipdir(p,:)=DFN.SetA.dipdir(p,:)-360;
    end
    if DFN.SetA.dipdir(p,:)<0
        DFN.SetA.dipdir(p,:)=360-abs(DFN.SetA.dipdir(p,:));
    end
end
end
DFN.SetA.normalx=sind(DFN.SetA.dipdir).*cosd(90-DFN.SetA.dip);
DFN.SetA.normaly=cosd(DFN.SetA.dipdir).*cosd(90-DFN.SetA.dip);
DFN.SetA.normalz=sind(90-DFN.SetA.dip);

    circles1(:,1)=DFN.SetA.xcenter;
    circles1(:,2)=DFN.SetA.ycenter;
	circles1(:,3)=DFN.SetA.zcenter;
    circles1(:,4)=DFN.SetA.radius;
    circles1(:,5)=DFN.SetA.normalx;
    circles1(:,6)=DFN.SetA.normaly;
    circles1(:,7)=DFN.SetA.normalz;
    


else
end
if numbersets==2 || numbersets==3 || numbersets==4 || numbersets==5    
DFN.SetB.dipdirmean=str2double(get (handles.editDFN_dipdirB,'String'));
DFN.SetB.dipmean=str2double(get (handles.editDFN_dipB,'String'));
DFN.SetB.orientstandard=str2double(get (handles.editDFN_stdB,'String'));
DFN.SetB.rmean=(str2double(get (handles.editDFN_diaB,'String')))/2;
DFN.SetB.onedimfreq=1/(str2double(get (handles.editDFN_spaceB,'String')));

%Calculation of number of joints
mean_radius=DFN.SetB.rmean;
onedimfreq=DFN.SetB.onedimfreq;       
% tolerance value within which the simulated 1D frequency should fall.
tolval = 1e-3;
% no. of imaginary scanlines for which 1D frequency is calculated. the larger the number
% the longer the simulation will take.
nsims = 100;
% two ends of search space for number of circles. 
maxguess=(onedimfreq*2.5)/((mean_radius^2)*pi)*xup*yup*zup;  %usually approx. 5 times larger than the final number of circles
ncircles = [1 maxguess];
freq = 0;

while abs(onedimfreq - freq) > tolval && ncircles(1)~=ncircles(2)
    nc = round(mean(ncircles));
    pproc = xup*rand(nc, 1) + 1i*zup*rand(nc, 1);
    circle_r = exprnd(mean_radius, nc, 1);
    
    lines_loc = xup*rand(1,nsims) + 1i*zup*rand(1, nsims);
    dist_mat = abs(repmat(pproc, 1, nsims) - repmat(lines_loc, nc, 1));
    
    intersect_circles = sum((dist_mat-repmat(circle_r, 1, nsims))<0);
    freq = mean(intersect_circles/zup);
    
    if freq>onedimfreq
        ncircles(end) = nc;
    else
        ncircles(1) = nc;
    end
end
DFN.SetB.njoint=ncircles(1);


% 3D joint point generation (Poisson process)
DFN.SetB.xcenter=xlo+(xup-xlo)*rand(DFN.SetB.njoint,1);
DFN.SetB.ycenter=ylo+(yup-ylo)*rand(DFN.SetB.njoint,1);
DFN.SetB.zcenter=zlo+(zup-ylo)*rand(DFN.SetB.njoint,1);

%Normal vectors
DFN.SetB.meannormalx=sind(DFN.SetB.dipdirmean).*cosd(90-DFN.SetB.dipmean);
DFN.SetB.meannormaly=cosd(DFN.SetB.dipdirmean).*cosd(90-DFN.SetB.dipmean);
DFN.SetB.meannormalz=sind(90-DFN.SetB.dipmean);
DFN.SetB.normalx=normrnd(DFN.SetB.meannormalx,(cosd(90-DFN.SetB.orientstandard)),DFN.SetB.njoint,1);
DFN.SetB.normaly=normrnd(DFN.SetB.meannormaly,(cosd(90-DFN.SetB.orientstandard)),DFN.SetB.njoint,1);
DFN.SetB.dip=normrnd(DFN.SetB.dipmean,DFN.SetB.orientstandard,DFN.SetB.njoint,1);
DFN.SetB.normalz=cosd(DFN.SetB.dip);
radii=exprnd(DFN.SetB.rmean,DFN.SetB.njoint,1);
radii_t=(randsample(radii(radii>0.1),DFN.SetB.njoint,1));  %Truncation of joints smaller than 20 cm in diameter
meanRT=mean(radii_t);
meanR=DFN.SetB.rmean;

while meanRT > (meanR*1.025)
radii=radii.*0.99;
    radii_t=(randsample(radii(radii>0.1),DFN.SetB.njoint,1));
meanRT=mean(radii_t);    
end
DFN.SetB.radius=radii_t;





clear radii radii_t meanR meanRT
u=numbersets;
for z=1:u
for p=1:length(DFN.SetB.normalx)
    if DFN.SetB.normalx(p,:)>1
        DFN.SetB.normalx(p,:)=1-(DFN.SetB.normalx(p,:)-1);
    end
    if DFN.SetB.normalx(p,:)<-1
        DFN.SetB.normalx(p,:)=-1+(-1-DFN.SetB.normalx(p,:));
    end
    if DFN.SetB.normaly(p,:)>1
        DFN.SetB.normaly(p,:)=1-(DFN.SetB.normaly(p,:)-1);
    end
    if DFN.SetB.normaly(p,:)<-1
        DFN.SetB.normaly(p,:)=-1+(-1-DFN.SetB.normaly(p,:));
    end
        if DFN.SetB.normalz(p,:)>1
        DFN.SetB.normalz(p,:)=-1+(DFN.SetB.normalz(p,:)-1);
       
        end
    if DFN.SetB.normalz(p,:)<-1
        DFN.SetB.normalz(p,:)=1-(DFN.SetB.normalz(p,:)+1);

    end
	
	end
end


z=1;
for z=1:u

DFN.SetB.TEMP.lowerx=DFN.SetB.normalx.*sign(DFN.SetB.normalz);
DFN.SetB.TEMP.lowery=DFN.SetB.normaly.*sign(DFN.SetB.normalz);
DFN.SetB.TEMP.lowerz=DFN.SetB.normalz.*sign(DFN.SetB.normalz);



if DFN.SetB.TEMP.lowerx==0 & DFN.SetB.TEMP.lowery==0
    DFN.SetB.TEMP.poleazimuth=0
elseif DFN.SetB.TEMP.lowerx<0 & DFN.SetB.TEMP.lowery>=0
    DFN.SetB.TEMP.poleazimuth=450-(atan2d(DFN.SetB.TEMP.lowery,DFN.SetB.TEMP.lowerx));
else
    DFN.SetB.TEMP.poleazimuth=90-(atan2d(DFN.SetB.TEMP.lowery,DFN.SetB.TEMP.lowerx));
end

DFN.SetB.TEMP.poleplunge=90-acosd((DFN.SetB.TEMP.lowerz));
DFN.SetB.dipdir=DFN.SetB.TEMP.poleazimuth;


% Simulation of the Radius distribution
% uniform distribution
DFN.SetB.dip=90-DFN.SetB.TEMP.poleplunge;
end
% convert dip dir/dip angle of joint in normal vectors
u=numbersets;
for z=1:u
for p=1:length(DFN.SetB.dip)
    if DFN.SetB.dipdir(p,:)>=360
        DFN.SetB.dipdir(p,:)=DFN.SetB.dipdir(p,:)-360;
    end
    if DFN.SetB.dipdir(p,:)<0
        DFN.SetB.dipdir(p,:)=360-abs(DFN.SetB.dipdir(p,:));
    end
end
end
DFN.SetB.normalx=sind(DFN.SetB.dipdir).*cosd(90-DFN.SetB.dip);
DFN.SetB.normaly=cosd(DFN.SetB.dipdir).*cosd(90-DFN.SetB.dip);
DFN.SetB.normalz=sind(90-DFN.SetB.dip);

    circles2(:,1)=DFN.SetB.xcenter;
    circles2(:,2)=DFN.SetB.ycenter;
	circles2(:,3)=DFN.SetB.zcenter;
    circles2(:,4)=DFN.SetB.radius;
    circles2(:,5)=DFN.SetB.normalx;
    circles2(:,6)=DFN.SetB.normaly;
    circles2(:,7)=DFN.SetB.normalz;    
else
end
if numbersets==3 || numbersets==4 || numbersets==5      
    DFN.SetC.dipdirmean=str2double(get (handles.editDFN_dipdirC,'String'));
DFN.SetC.dipmean=str2double(get (handles.editDFN_dipC,'String'));
DFN.SetC.orientstandard=str2double(get (handles.editDFN_stdC,'String'));
DFN.SetC.rmean=(str2double(get (handles.editDFN_diaC,'String')))/2;
DFN.SetC.onedimfreq=1/(str2double(get (handles.editDFN_spaceC,'String')));

%Calculation of number of joints
mean_radius=DFN.SetC.rmean;
onedimfreq=DFN.SetC.onedimfreq;       
% tolerance value within which the simulated 1D frequency should fall.
tolval = 1e-3;
% no. of imaginary scanlines for which 1D frequency is calculated. the larger the number
% the longer the simulation will take.
nsims = 100;
% two ends of search space for number of circles. 
maxguess=(onedimfreq*2.5)/((mean_radius^2)*pi)*xup*yup*zup;  %usually approx. 5 times larger than the final number of circles
ncircles = [1 maxguess];
freq = 0;

while abs(onedimfreq - freq) > tolval && ncircles(1)~=ncircles(2)
    nc = round(mean(ncircles));
    pproc = xup*rand(nc, 1) + 1i*zup*rand(nc, 1);
    circle_r = exprnd(mean_radius, nc, 1);
    
    lines_loc = xup*rand(1,nsims) + 1i*zup*rand(1, nsims);
    dist_mat = abs(repmat(pproc, 1, nsims) - repmat(lines_loc, nc, 1));
    
    intersect_circles = sum((dist_mat-repmat(circle_r, 1, nsims))<0);
    freq = mean(intersect_circles/zup);
    
    if freq>onedimfreq
        ncircles(end) = nc;
    else
        ncircles(1) = nc;
    end
end
DFN.SetC.njoint=ncircles(1);


% 3D joint point generation (Poisson process)
DFN.SetC.xcenter=xlo+(xup-xlo)*rand(DFN.SetC.njoint,1);
DFN.SetC.ycenter=ylo+(yup-ylo)*rand(DFN.SetC.njoint,1);
DFN.SetC.zcenter=zlo+(zup-ylo)*rand(DFN.SetC.njoint,1);

%Normal vectors
DFN.SetC.meannormalx=sind(DFN.SetC.dipdirmean).*cosd(90-DFN.SetC.dipmean);
DFN.SetC.meannormaly=cosd(DFN.SetC.dipdirmean).*cosd(90-DFN.SetC.dipmean);
DFN.SetC.meannormalz=sind(90-DFN.SetC.dipmean);
DFN.SetC.normalx=normrnd(DFN.SetC.meannormalx,(cosd(90-DFN.SetC.orientstandard)),DFN.SetC.njoint,1);
DFN.SetC.normaly=normrnd(DFN.SetC.meannormaly,(cosd(90-DFN.SetC.orientstandard)),DFN.SetC.njoint,1);
DFN.SetC.dip=normrnd(DFN.SetC.dipmean,DFN.SetC.orientstandard,DFN.SetC.njoint,1);
DFN.SetC.normalz=cosd(DFN.SetC.dip);
radii=exprnd(DFN.SetC.rmean,DFN.SetC.njoint,1);
radii_t=(randsample(radii(radii>0.1),DFN.SetC.njoint,1));  %Truncation of joints smaller than 20 cm in diameter
meanRT=mean(radii_t);
meanR=DFN.SetC.rmean;

while meanRT > (meanR*1.025)
radii=radii.*0.99;
    radii_t=(randsample(radii(radii>0.1),DFN.SetC.njoint,1));
meanRT=mean(radii_t);    
end
DFN.SetC.radius=radii_t;





clear radii radii_t meanR meanRT
u=numbersets;
for z=1:u
for p=1:length(DFN.SetC.normalx)
    if DFN.SetC.normalx(p,:)>1
        DFN.SetC.normalx(p,:)=1-(DFN.SetC.normalx(p,:)-1);
    end
    if DFN.SetC.normalx(p,:)<-1
        DFN.SetC.normalx(p,:)=-1+(-1-DFN.SetC.normalx(p,:));
    end
    if DFN.SetC.normaly(p,:)>1
        DFN.SetC.normaly(p,:)=1-(DFN.SetC.normaly(p,:)-1);
    end
    if DFN.SetC.normaly(p,:)<-1
        DFN.SetC.normaly(p,:)=-1+(-1-DFN.SetC.normaly(p,:));
    end
        if DFN.SetC.normalz(p,:)>1
        DFN.SetC.normalz(p,:)=-1+(DFN.SetC.normalz(p,:)-1);
       
        end
    if DFN.SetC.normalz(p,:)<-1
        DFN.SetC.normalz(p,:)=1-(DFN.SetC.normalz(p,:)+1);

    end
	
	end
end


z=1;
for z=1:u

DFN.SetC.TEMP.lowerx=DFN.SetC.normalx.*sign(DFN.SetC.normalz);
DFN.SetC.TEMP.lowery=DFN.SetC.normaly.*sign(DFN.SetC.normalz);
DFN.SetC.TEMP.lowerz=DFN.SetC.normalz.*sign(DFN.SetC.normalz);



if DFN.SetC.TEMP.lowerx==0 & DFN.SetC.TEMP.lowery==0
    DFN.SetC.TEMP.poleazimuth=0
elseif DFN.SetC.TEMP.lowerx<0 & DFN.SetC.TEMP.lowery>=0
    DFN.SetC.TEMP.poleazimuth=450-(atan2d(DFN.SetC.TEMP.lowery,DFN.SetC.TEMP.lowerx));
else
    DFN.SetC.TEMP.poleazimuth=90-(atan2d(DFN.SetC.TEMP.lowery,DFN.SetC.TEMP.lowerx));
end

DFN.SetC.TEMP.poleplunge=90-acosd((DFN.SetC.TEMP.lowerz));
DFN.SetC.dipdir=DFN.SetC.TEMP.poleazimuth;


% Simulation of the Radius distribution
% uniform distribution
DFN.SetC.dip=90-DFN.SetC.TEMP.poleplunge;
end
% convert dip dir/dip angle of joint in normal vectors
u=numbersets;
for z=1:u
for p=1:length(DFN.SetC.dip)
    if DFN.SetC.dipdir(p,:)>=360
        DFN.SetC.dipdir(p,:)=DFN.SetC.dipdir(p,:)-360;
    end
    if DFN.SetC.dipdir(p,:)<0
        DFN.SetC.dipdir(p,:)=360-abs(DFN.SetC.dipdir(p,:));
    end
end
end
DFN.SetC.normalx=sind(DFN.SetC.dipdir).*cosd(90-DFN.SetC.dip);
DFN.SetC.normaly=cosd(DFN.SetC.dipdir).*cosd(90-DFN.SetC.dip);
DFN.SetC.normalz=sind(90-DFN.SetC.dip);

    circles3(:,1)=DFN.SetC.xcenter;
    circles3(:,2)=DFN.SetC.ycenter;
	circles3(:,3)=DFN.SetC.zcenter;
    circles3(:,4)=DFN.SetC.radius;
    circles3(:,5)=DFN.SetC.normalx;
    circles3(:,6)=DFN.SetC.normaly;
    circles3(:,7)=DFN.SetC.normalz;
else
end
if numbersets==4 || numbersets==5    
    DFN.SetD.dipdirmean=str2double(get (handles.editDFN_dipdirD,'String'));
DFN.SetD.dipmean=str2double(get (handles.editDFN_dipD,'String'));
DFN.SetD.orientstandard=str2double(get (handles.editDFN_stdD,'String'));
DFN.SetD.rmean=(str2double(get (handles.editDFN_diaD,'String')))/2;
DFN.SetD.onedimfreq=1/(str2double(get (handles.editDFN_spaceD,'String')));

%Calculation of number of joints
mean_radius=DFN.SetD.rmean;
onedimfreq=DFN.SetD.onedimfreq;       
% tolerance value within which the simulated 1D frequency should fall.
tolval = 1e-3;
% no. of imaginary scanlines for which 1D frequency is calculated. the larger the number
% the longer the simulation will take.
nsims = 100;
% two ends of search space for number of circles. 
maxguess=(onedimfreq*2.5)/((mean_radius^2)*pi)*xup*yup*zup;  %usually approx. 5 times larger than the final number of circles
ncircles = [1 maxguess];
freq = 0;

while abs(onedimfreq - freq) > tolval && ncircles(1)~=ncircles(2)
    nc = round(mean(ncircles));
    pproc = xup*rand(nc, 1) + 1i*zup*rand(nc, 1);
    circle_r = exprnd(mean_radius, nc, 1);
    
    lines_loc = xup*rand(1,nsims) + 1i*zup*rand(1, nsims);
    dist_mat = abs(repmat(pproc, 1, nsims) - repmat(lines_loc, nc, 1));
    
    intersect_circles = sum((dist_mat-repmat(circle_r, 1, nsims))<0);
    freq = mean(intersect_circles/zup);
    
    if freq>onedimfreq
        ncircles(end) = nc;
    else
        ncircles(1) = nc;
    end
end
DFN.SetD.njoint=ncircles(1);


% 3D joint point generation (Poisson process)
DFN.SetD.xcenter=xlo+(xup-xlo)*rand(DFN.SetD.njoint,1);
DFN.SetD.ycenter=ylo+(yup-ylo)*rand(DFN.SetD.njoint,1);
DFN.SetD.zcenter=zlo+(zup-ylo)*rand(DFN.SetD.njoint,1);

%Normal vectors
DFN.SetD.meannormalx=sind(DFN.SetD.dipdirmean).*cosd(90-DFN.SetD.dipmean);
DFN.SetD.meannormaly=cosd(DFN.SetD.dipdirmean).*cosd(90-DFN.SetD.dipmean);
DFN.SetD.meannormalz=sind(90-DFN.SetD.dipmean);
DFN.SetD.normalx=normrnd(DFN.SetD.meannormalx,(cosd(90-DFN.SetD.orientstandard)),DFN.SetD.njoint,1);
DFN.SetD.normaly=normrnd(DFN.SetD.meannormaly,(cosd(90-DFN.SetD.orientstandard)),DFN.SetD.njoint,1);
DFN.SetD.dip=normrnd(DFN.SetD.dipmean,DFN.SetD.orientstandard,DFN.SetD.njoint,1);
DFN.SetD.normalz=cosd(DFN.SetD.dip);
radii=exprnd(DFN.SetD.rmean,DFN.SetD.njoint,1);
radii_t=(randsample(radii(radii>0.1),DFN.SetD.njoint,1));  %Truncation of joints smaller than 20 cm in diameter
meanRT=mean(radii_t);
meanR=DFN.SetD.rmean;

while meanRT > (meanR*1.025)
radii=radii.*0.99;
    radii_t=(randsample(radii(radii>0.1),DFN.SetD.njoint,1));
meanRT=mean(radii_t);    
end
DFN.SetD.radius=radii_t;





clear radii radii_t meanR meanRT
u=numbersets;
for z=1:u
for p=1:length(DFN.SetD.normalx)
    if DFN.SetD.normalx(p,:)>1
        DFN.SetD.normalx(p,:)=1-(DFN.SetD.normalx(p,:)-1);
    end
    if DFN.SetD.normalx(p,:)<-1
        DFN.SetD.normalx(p,:)=-1+(-1-DFN.SetD.normalx(p,:));
    end
    if DFN.SetD.normaly(p,:)>1
        DFN.SetD.normaly(p,:)=1-(DFN.SetD.normaly(p,:)-1);
    end
    if DFN.SetD.normaly(p,:)<-1
        DFN.SetD.normaly(p,:)=-1+(-1-DFN.SetD.normaly(p,:));
    end
        if DFN.SetD.normalz(p,:)>1
        DFN.SetD.normalz(p,:)=-1+(DFN.SetD.normalz(p,:)-1);
       
        end
    if DFN.SetD.normalz(p,:)<-1
        DFN.SetD.normalz(p,:)=1-(DFN.SetD.normalz(p,:)+1);

    end
	
	end
end


z=1;
for z=1:u

DFN.SetD.TEMP.lowerx=DFN.SetD.normalx.*sign(DFN.SetD.normalz);
DFN.SetD.TEMP.lowery=DFN.SetD.normaly.*sign(DFN.SetD.normalz);
DFN.SetD.TEMP.lowerz=DFN.SetD.normalz.*sign(DFN.SetD.normalz);



if DFN.SetD.TEMP.lowerx==0 & DFN.SetD.TEMP.lowery==0
    DFN.SetD.TEMP.poleazimuth=0
elseif DFN.SetD.TEMP.lowerx<0 & DFN.SetD.TEMP.lowery>=0
    DFN.SetD.TEMP.poleazimuth=450-(atan2d(DFN.SetD.TEMP.lowery,DFN.SetD.TEMP.lowerx));
else
    DFN.SetD.TEMP.poleazimuth=90-(atan2d(DFN.SetD.TEMP.lowery,DFN.SetD.TEMP.lowerx));
end

DFN.SetD.TEMP.poleplunge=90-acosd((DFN.SetD.TEMP.lowerz));
DFN.SetD.dipdir=DFN.SetD.TEMP.poleazimuth;


% Simulation of the Radius distribution
% uniform distribution
DFN.SetD.dip=90-DFN.SetD.TEMP.poleplunge;
end
% convert dip dir/dip angle of joint in normal vectors
u=numbersets;
for z=1:u
for p=1:length(DFN.SetD.dip)
    if DFN.SetD.dipdir(p,:)>=360
        DFN.SetD.dipdir(p,:)=DFN.SetD.dipdir(p,:)-360;
    end
    if DFN.SetD.dipdir(p,:)<0
        DFN.SetD.dipdir(p,:)=360-abs(DFN.SetD.dipdir(p,:));
    end
end
end
DFN.SetD.normalx=sind(DFN.SetD.dipdir).*cosd(90-DFN.SetD.dip);
DFN.SetD.normaly=cosd(DFN.SetD.dipdir).*cosd(90-DFN.SetD.dip);
DFN.SetD.normalz=sind(90-DFN.SetD.dip);

    circles4(:,1)=DFN.SetD.xcenter;
    circles4(:,2)=DFN.SetD.ycenter;
	circles4(:,3)=DFN.SetD.zcenter;
    circles4(:,4)=DFN.SetD.radius;
    circles4(:,5)=DFN.SetD.normalx;
    circles4(:,6)=DFN.SetD.normaly;
    circles4(:,7)=DFN.SetD.normalz;
else
end
if numbersets==5 
    DFN.SetE.dipdirmean=str2double(get (handles.editDFN_dipdirE,'String'));
DFN.SetE.dipmean=str2double(get (handles.editDFN_dipE,'String'));
DFN.SetE.orientstandard=str2double(get (handles.editDFN_stdE,'String'));
DFN.SetE.rmean=(str2double(get (handles.editDFN_diaE,'String')))/2;
DFN.SetE.onedimfreq=1/(str2double(get (handles.editDFN_spaceE,'String')));

%Calculation of number of joints
mean_radius=DFN.SetE.rmean;
onedimfreq=DFN.SetE.onedimfreq;       
% tolerance value within which the simulated 1D frequency should fall.
tolval = 1e-3;
% no. of imaginary scanlines for which 1D frequency is calculated. the larger the number
% the longer the simulation will take.
nsims = 100;
% two ends of search space for number of circles. 
maxguess=(onedimfreq*2.5)/((mean_radius^2)*pi)*xup*yup*zup;  %usually approx. 5 times larger than the final number of circles
ncircles = [1 maxguess];
freq = 0;

while abs(onedimfreq - freq) > tolval && ncircles(1)~=ncircles(2)
    nc = round(mean(ncircles));
    pproc = xup*rand(nc, 1) + 1i*zup*rand(nc, 1);
    circle_r = exprnd(mean_radius, nc, 1);
    
    lines_loc = xup*rand(1,nsims) + 1i*zup*rand(1, nsims);
    dist_mat = abs(repmat(pproc, 1, nsims) - repmat(lines_loc, nc, 1));
    
    intersect_circles = sum((dist_mat-repmat(circle_r, 1, nsims))<0);
    freq = mean(intersect_circles/zup);
    
    if freq>onedimfreq
        ncircles(end) = nc;
    else
        ncircles(1) = nc;
    end
end
DFN.SetE.njoint=ncircles(1);


% 3D joint point generation (Poisson process)
DFN.SetE.xcenter=xlo+(xup-xlo)*rand(DFN.SetE.njoint,1);
DFN.SetE.ycenter=ylo+(yup-ylo)*rand(DFN.SetE.njoint,1);
DFN.SetE.zcenter=zlo+(zup-ylo)*rand(DFN.SetE.njoint,1);

%Normal vectors
DFN.SetE.meannormalx=sind(DFN.SetE.dipdirmean).*cosd(90-DFN.SetE.dipmean);
DFN.SetE.meannormaly=cosd(DFN.SetE.dipdirmean).*cosd(90-DFN.SetE.dipmean);
DFN.SetE.meannormalz=sind(90-DFN.SetE.dipmean);
DFN.SetE.normalx=normrnd(DFN.SetE.meannormalx,(cosd(90-DFN.SetE.orientstandard)),DFN.SetE.njoint,1);
DFN.SetE.normaly=normrnd(DFN.SetE.meannormaly,(cosd(90-DFN.SetE.orientstandard)),DFN.SetE.njoint,1);
DFN.SetE.dip=normrnd(DFN.SetE.dipmean,DFN.SetE.orientstandard,DFN.SetE.njoint,1);
DFN.SetE.normalz=cosd(DFN.SetE.dip);
radii=exprnd(DFN.SetE.rmean,DFN.SetE.njoint,1);
radii_t=(randsample(radii(radii>0.1),DFN.SetE.njoint,1));  %Truncation of joints smaller than 20 cm in diameter
meanRT=mean(radii_t);
meanR=DFN.SetE.rmean;

while meanRT > (meanR*1.025)
radii=radii.*0.99;
    radii_t=(randsample(radii(radii>0.1),DFN.SetE.njoint,1));
meanRT=mean(radii_t);    
end
DFN.SetE.radius=radii_t;





clear radii radii_t meanR meanRT
u=numbersets;
for z=1:u
for p=1:length(DFN.SetE.normalx)
    if DFN.SetE.normalx(p,:)>1
        DFN.SetE.normalx(p,:)=1-(DFN.SetE.normalx(p,:)-1);
    end
    if DFN.SetE.normalx(p,:)<-1
        DFN.SetE.normalx(p,:)=-1+(-1-DFN.SetE.normalx(p,:));
    end
    if DFN.SetE.normaly(p,:)>1
        DFN.SetE.normaly(p,:)=1-(DFN.SetE.normaly(p,:)-1);
    end
    if DFN.SetE.normaly(p,:)<-1
        DFN.SetE.normaly(p,:)=-1+(-1-DFN.SetE.normaly(p,:));
    end
        if DFN.SetE.normalz(p,:)>1
        DFN.SetE.normalz(p,:)=-1+(DFN.SetE.normalz(p,:)-1);
       
        end
    if DFN.SetE.normalz(p,:)<-1
        DFN.SetE.normalz(p,:)=1-(DFN.SetE.normalz(p,:)+1);

    end
	
	end
end


z=1;
for z=1:u

DFN.SetE.TEMP.lowerx=DFN.SetE.normalx.*sign(DFN.SetE.normalz);
DFN.SetE.TEMP.lowery=DFN.SetE.normaly.*sign(DFN.SetE.normalz);
DFN.SetE.TEMP.lowerz=DFN.SetE.normalz.*sign(DFN.SetE.normalz);



if DFN.SetE.TEMP.lowerx==0 & DFN.SetE.TEMP.lowery==0
    DFN.SetE.TEMP.poleazimuth=0
elseif DFN.SetE.TEMP.lowerx<0 & DFN.SetE.TEMP.lowery>=0
    DFN.SetE.TEMP.poleazimuth=450-(atan2d(DFN.SetE.TEMP.lowery,DFN.SetE.TEMP.lowerx));
else
    DFN.SetE.TEMP.poleazimuth=90-(atan2d(DFN.SetE.TEMP.lowery,DFN.SetE.TEMP.lowerx));
end

DFN.SetE.TEMP.poleplunge=90-acosd((DFN.SetE.TEMP.lowerz));
DFN.SetE.dipdir=DFN.SetE.TEMP.poleazimuth;


% Simulation of the Radius distribution
% uniform distribution
DFN.SetE.dip=90-DFN.SetE.TEMP.poleplunge;
end
% convert dip dir/dip angle of joint in normal vectors
u=numbersets;
for z=1:u
for p=1:length(DFN.SetE.dip)
    if DFN.SetE.dipdir(p,:)>=360
        DFN.SetE.dipdir(p,:)=DFN.SetE.dipdir(p,:)-360;
    end
    if DFN.SetE.dipdir(p,:)<0
        DFN.SetE.dipdir(p,:)=360-abs(DFN.SetE.dipdir(p,:));
    end
end
end
DFN.SetE.normalx=sind(DFN.SetE.dipdir).*cosd(90-DFN.SetE.dip);
DFN.SetE.normaly=cosd(DFN.SetE.dipdir).*cosd(90-DFN.SetE.dip);
DFN.SetE.normalz=sind(90-DFN.SetE.dip);

    circles5(:,1)=DFN.SetE.xcenter;
    circles5(:,2)=DFN.SetE.ycenter;
	circles5(:,3)=DFN.SetE.zcenter;
    circles5(:,4)=DFN.SetE.radius;
    circles5(:,5)=DFN.SetE.normalx;
    circles5(:,6)=DFN.SetE.normaly;
    circles5(:,7)=DFN.SetE.normalz;
else
end  
  
    
if numbersets==1    
circles=vertcat(circles1); 
elseif numbersets==2
circles=vertcat(circles1,circles2); 
elseif numbersets==3
circles=vertcat(circles1,circles2,circles3); 
elseif numbersets==4
circles=vertcat(circles1,circles2,circles3,circles4); 
elseif numbersets==5
circles=vertcat(circles1,circles2,circles3,circles4,circles5); 
end

circlesboundxcenter=[dimensions/2, dimensions/2, 0, dimensions, dimensions/2, dimensions/2];
circlesboundycenter=[dimensions/2, dimensions/2, dimensions/2, dimensions/2, 0, dimensions];
circlesboundzcenter=[0, dimensions, dimensions/2, dimensions/2, dimensions/2, dimensions/2];
circlesboundradius=[dimensions*5, dimensions*5, dimensions*5, dimensions*5, dimensions*5, dimensions*5];
circlesboundnormalx=[0, 0, 1, 1, 0, 0];
circlesboundnormaly=[0,0,0,0,1,1];
circlesboundnormalz=[1,1,0,0,0,0];
boundary=vertcat(circlesboundxcenter,circlesboundycenter,circlesboundzcenter,circlesboundradius,circlesboundnormalx,circlesboundnormaly,circlesboundnormalz);




circles=vertcat(boundary',circles);
numberofcircles=length(circles);
xofcircles=circles(:,1);
yofcircles=circles(:,2);
zofcircles=circles(:,3);
radiiofcircles=circles(:,4);
xnormalvectors=circles(:,5);
ynormalvectors=circles(:,6);
znormalvectors=circles(:,7);
	
	
n_ = numberofcircles;
R1 = radiiofcircles(:);
C1 = [xofcircles(:) yofcircles(:) zofcircles(:)];
N1 = [xnormalvectors(:) ynormalvectors(:) znormalvectors(:)];
N1 = N1 ./ repmat(sqrt(dot(N1,N1,2)),1,3);
number=size(circles,1);
[R,idx] = sort(R1,'ascend');
C = C1(idx,:);
N = N1(idx,:);
dlmwrite ('n.txt', N', 'delimiter', ' ');
dlmwrite ('r.txt', R, 'delimiter', ' ');
dlmwrite ('c.txt', C', 'delimiter', ' ');
clear Volumes fid
geomfile=[dimensions dimensions dimensions;0 0 0;dimrast dimrast dimrast; number number number]';
dlmwrite ('geom.txt', geomfile', 'delimiter', '\t');

command=sprintf('cells.exe %d 2.4', number);

status = dos(command);
set(handles.textDFN_calc,'Visible','off');
set(handles.text296,'Visible','on');
set(handles.text297,'Visible','on');

fid = fopen('result.txt');
Volumes = fscanf(fid, '%g', Inf);
fclose(fid);    












%% MAJOR CALLBACKS (THE HAVE INTERACTION WITH GUI CHANGE
%===========================================================

% --- Executes on button press in pushbutton_filelocation.
function pushbutton_filelocation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_filelocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Loads *.csv into the GUI
[filename pathname]=uigetfile({'*.csv*'},'Select file');
fullpathname=strcat(pathname, filename);
set(handles.edit_filelocation, 'String', fullpathname);

function pushbutton_general_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_general (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.uipanel_dfn,'Visible','off');
set(handles.panel_general,'Visible','on');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.panel_rmc,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');

function pushbutton_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','on');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.panel_rmc,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');


% --- Executes on button press in pushbutton_savefig.
function pushbutton_savefig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_savefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GUI_fig_children=get(gcf,'children');
Fig_Axes=findobj(GUI_fig_children,'type','Axes');
fig=figure;ax=axes;clf;
new_handle=copyobj(Fig_Axes,fig);
set(gca,'ActivePositionProperty','outerposition')
set(gca,'Units','normalized')
set(gca,'OuterPosition',[0 0 1 1])
set(gca,'position',[0.1300 0.1100 0.7750 0.8150])

set(handles.popupmenu_clustersets,'Visible','on');
set(handles.pushbutton_calculate,'Visible','on');



function edit_anglebc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_anglebc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_anglebc as text
%        str2double(get(hObject,'String')) returns contents of edit_anglebc as a double

set(axesblockvol,'Tag','axes_blockvol');

% --- Executes on button press in pushbutton_merge.
function pushbutton_merge_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','on');
set(handles.panel_rmc,'Visible','off');
set(handles.uipanel_dfn,'Visible','off');

function pushbutton_dfn_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dfn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_dfn,'Visible','on');
set(handles.panel_general,'Visible','off');
set(handles.panel_clustering,'Visible','off');
set(handles.panel_statistics,'Visible','off');
set(handles.panel_volume,'Visible','off');
set(handles.uipanel_merge,'Visible','off');
set(handles.panel_rmc,'Visible','off');



% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filenamesl1 pathnamesl1]=uigetfile({'*.mat*'},'Select file');
fullpathnamesl1=strcat(pathnamesl1, filenamesl1);
set(handles.edit_data1, 'String', fullpathnamesl1);

% --- Executes on button press in pushbutton_locsl2.
function pushbutton_locsl2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_locsl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filenamesl2 pathnamesl2]=uigetfile({'*.mat*'},'Select file');
fullpathnamesl2=strcat(pathnamesl2, filenamesl2);
set(handles.edit_data2, 'String', fullpathnamesl2);



% --- Executes on button press in text187.
function text187_Callback(hObject, eventdata, handles)
% hObject    handle to text187 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
jlgsi=str2double(get(handles.edit_rmc_jl,'String'));
jrgsi=str2double(get(handles.edit_rmc_jr,'String'));
jagsi=str2double(get(handles.edit_rmc_ja,'String'));


jcgsi=jlgsi*jrgsi/jagsi;

set(handles.edit_rmc_jc, 'String', jcgsi);










%% Minor Callbacks (the don't have any major calculation but you can modify according to use
%============================================================================================

function edit_filelocation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filelocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_filelocation as text
%        str2double(get(hObject,'String')) returns contents of edit_filelocation as a double


function edit_sltrend_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sltrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_sltrend as text
%        str2double(get(hObject,'String')) returns contents of edit_sltrend as a double
% --- Executes during object creation, after setting all properties.


function edit_slplunge_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slplunge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_slplunge as text
%        str2double(get(hObject,'String')) returns contents of edit_slplunge as a double



function edit_sllength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sllength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_sllength as text
%        str2double(get(hObject,'String')) returns contents of edit_sllength as a double


% --- Executes when entered data in editable cell(s) in preview_table.
function preview_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to preview_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_newscanline.
function pushbutton_newscanline_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_newscanline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double



% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double




% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double



function edit_answer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_answer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_answer as text
%        str2double(get(hObject,'String')) returns contents of edit_answer as a double




function edit_setclusters_Callback(hObject, eventdata, handles)
% hObject    handle to edit_setclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_setclusters as text
%        str2double(get(hObject,'String')) returns contents of edit_setclusters as a double





function edit_clusterokanswer_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clusterokanswer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_clusterokanswer as text
%        str2double(get(hObject,'String')) returns contents of edit_clusterokanswer as a double



function edit_clusterok_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clusterok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_clusterok as text
%        str2double(get(hObject,'String')) returns contents of edit_clusterok as a double




function edit_statusmanual_Callback(hObject, eventdata, handles)
% hObject    handle to edit_statusmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_statusmanual as text
%        str2double(get(hObject,'String')) returns contents of edit_statusmanual as a double





function edit_finish_Callback(hObject, eventdata, handles)
% hObject    handle to text_clustermeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of text_clustermeans as text
%        str2double(get(hObject,'String')) returns contents of text_clustermeans as a double



% --- Executes on selection change in popupmenu_clustersets.
function popupmenu_clustersets_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_clustersets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_clustersets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_clustersets



% --- Executes on selection change in popupmenu_changeclusters.
function popupmenu_changeclusters_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_changeclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_changeclusters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_changeclusters




function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in check_box_scan_spacing.
function check_box_scan_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_spacing


% --- Executes on button press in check_box__scan_tracelength.
function check_box__scan_tracelength_Callback(hObject, eventdata, handles)
% hObject    handle to check_box__scan_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box__scan_tracelength


% --- Executes on button press in check_box_scan_jc.
function check_box_scan_jc_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_jc


% --- Executes on button press in check_box_scan_phi.
function check_box_scan_phi_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_phi


% --- Executes on button press in check_box_scan_jointsizefactor.
function check_box_scan_jointsizefactor_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_jointsizefactor


% --- Executes on button press in check_box_scan_roughness.
function check_box_scan_roughness_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_roughness


% --- Executes on button press in check_box_scan_alteration.
function check_box_scan_alteration_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_scan_alteration


% --- Executes on button press in check_box_sets_spacing.
function check_box_sets_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_spacing


% --- Executes on button press in check_box_sets_tracelengths.
function check_box_sets_tracelengths_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_tracelengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_tracelengths


% --- Executes on button press in check_box_sets_jc.
function check_box_sets_jc_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_jc


% --- Executes on button press in check_box_sets_phi.
function check_box_sets_phi_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_phi


% --- Executes on button press in check_box_sets_jointsizefactor.
function check_box_sets_jointsizefactor_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_jointsizefactor


% --- Executes on button press in check_box_sets_roughness.
function check_box_sets_roughness_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_roughness


% --- Executes on button press in check_box_sets_alteration.
function check_box_sets_alteration_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_box_sets_alteration


% --- Executes on button press in check_hist_scan_spacing.
function check_hist_scan_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_spacing


% --- Executes on button press in check_hist_scan_tracelength.
function check_hist_scan_tracelength_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_tracelength


% --- Executes on button press in check_hist_scan_jc.
function check_hist_scan_jc_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_jc


% --- Executes on button press in check_hist_scan_phi.
function check_hist_scan_phi_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_phi


% --- Executes on button press in check_hist_scan_jointsizefactor.
function check_hist_scan_jointsizefactor_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_jointsizefactor


% --- Executes on button press in ckeck_hist_scan_roughness.
function ckeck_hist_scan_roughness_Callback(hObject, eventdata, handles)
% hObject    handle to ckeck_hist_scan_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of ckeck_hist_scan_roughness


% --- Executes on button press in check_hist_scan_alteration.
function check_hist_scan_alteration_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_alteration

% --- Executes on button press in check_hist_sets_spacing.
function check_hist_sets_spacing_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_spacing


% --- Executes on button press in check_hist_sets_tracelength.
function check_hist_sets_tracelength_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_tracelength


% --- Executes on button press in check_hist_sets_jc.
function check_hist_sets_jc_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_jc


% --- Executes on button press in check_hist_sets_phi.
function check_hist_sets_phi_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_phi


% --- Executes on button press in check_hist_sets_jointsizefactor.
function check_hist_sets_jointsizefactor_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_jointsizefactor


% --- Executes on button press in check_hist_sets_roughness.
function check_hist_sets_roughness_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_roughness


% --- Executes on button press in check_hist_sets_alteration.
function check_hist_sets_alteration_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_alteration


% --- Executes on button press in check_rqd_hist.
function check_rqd_hist_Callback(hObject, eventdata, handles)
% hObject    handle to check_rqd_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rqd_hist


% --- Executes on button press in check_rqd_distance.
function check_rqd_distance_Callback(hObject, eventdata, handles)
% hObject    handle to check_rqd_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rqd_distance


% --- Executes on button press in check_rose_jc.
function check_rose_jc_Callback(hObject, eventdata, handles)
% hObject    handle to check_rose_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rose_jc


% --- Executes on button press in check_rose_phi.
function check_rose_phi_Callback(hObject, eventdata, handles)
% hObject    handle to check_rose_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rose_phi


% --- Executes on button press in check_rose_waviness.
function check_rose_waviness_Callback(hObject, eventdata, handles)
% hObject    handle to check_rose_waviness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rose_waviness


% --- Executes on button press in check_rose_roughness.
function check_rose_roughness_Callback(hObject, eventdata, handles)
% hObject    handle to check_rose_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rose_roughness


% --- Executes on button press in check_rose_alteration.
function check_rose_alteration_Callback(hObject, eventdata, handles)
% hObject    handle to check_rose_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of check_rose_alteration




function edit_angleab_Callback(hObject, eventdata, handles)
% hObject    handle to edit_angleab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_angleab as text
%        str2double(get(hObject,'String')) returns contents of edit_angleab as a double




function edit_angleac_Callback(hObject, eventdata, handles)
% hObject    handle to edit_angleac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_angleac as text
%        str2double(get(hObject,'String')) returns contents of edit_angleac as a double







function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double

function edit_spacing_seta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_spacing_seta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_spacing_seta as text
%        str2double(get(hObject,'String')) returns contents of edit_spacing_seta as a double





function edit_spacing_setb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_spacing_setb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_spacing_setb as text
%        str2double(get(hObject,'String')) returns contents of edit_spacing_setb as a double





function edit_spacing_setc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_spacing_setc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_spacing_setc as text
%        str2double(get(hObject,'String')) returns contents of edit_spacing_setc as a double


function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double






function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double




function edit_blockshapefactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_blockshapefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_blockshapefactor as text
%        str2double(get(hObject,'String')) returns contents of edit_blockshapefactor as a double





% --- Executes on button press in radio_3negexp.
function radio_3negexp_Callback(hObject, eventdata, handles)
% hObject    handle to radio_3negexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_3negexp


% --- Executes on button press in check_3uniform.
function check_3uniform_Callback(hObject, eventdata, handles)
% hObject    handle to check_3uniform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_3uniform


% --- Executes on button press in check_2negexp1uni.
function check_2negexp1uni_Callback(hObject, eventdata, handles)
% hObject    handle to check_2negexp1uni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_2negexp1uni


% --- Executes on button press in check_3lognorm.
function check_3lognorm_Callback(hObject, eventdata, handles)
% hObject    handle to check_3lognorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_3lognorm


% --- Executes on button press in check_3fractal.
function check_3fractal_Callback(hObject, eventdata, handles)
% hObject    handle to check_3fractal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_3fractal



function edit_persistence_Callback(hObject, eventdata, handles)
% hObject    handle to edit_persistence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_persistence as text
%        str2double(get(hObject,'String')) returns contents of edit_persistence as a double




% --- Executes on button press in pushbutton_calcvolpalm.
function pushbutton_calcvolpalm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_calcvolpalm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




 
function edit_data1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_data1 as text
%        str2double(get(hObject,'String')) returns contents of edit_data1 as a double




function edit_data2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_data2 as text
%        str2double(get(hObject,'String')) returns contents of edit_data2 as a double



function edit_data3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_data3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_data3 as text
%        str2double(get(hObject,'String')) returns contents of edit_data3 as a double


% --- Executes on button press in pushbutton_locsl3.
function pushbutton_locsl3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_locsl3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function edit_rmc_smr_rmrb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_rmrb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_smr_rmrb as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_smr_rmrb as a double






function edit_rmc_smr_dipdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_smr_dipdir as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_smr_dipdir as a double





function edit_rmc_ucs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_ucs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_ucs as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_ucs as a double




function edit_rmc_rmi_jc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmi_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmi_jc as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmi_jc as a double




function edit_rmc_rmi_vb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_vb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_vb as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_vb as a double




function edit_rmc_jc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jc as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jc as a double




function edit_rmc_vb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_vb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_vb as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_vb as a double



function edit_rmc_rqd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rqd as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rqd as a double




function edit_rmc_jn_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jn as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jn as a double




function edit_rmc_jc_q_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jc_q as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jc_q as a double



function edit_rmc_rmr_ucs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_ucs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmr_ucs as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmr_ucs as a double






function edit_rmc_rmr_rqd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rqd as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rqd as a double







function edit_rmc_rmr_jc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmr_jc as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmr_jc as a double




function edit_rmc_rmr_wc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmr_wc as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmr_wc as a double




function edit_rmc_rmr_strike_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_strike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmr_strike as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmr_strike as a double






% --- Executes on selection change in popup_rmc_rmr_wc.
function popup_rmc_rmr_wc_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rmc_rmr_wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rmc_rmr_wc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rmc_rmr_wc





function edit_rmc_rmr_dip_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rmr_dip as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rmr_dip as a double




% --- Executes on selection change in popup_rmr_orientation.
function popup_rmr_orientation_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rmr_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rmr_orientation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rmr_orientation





function edit_rmc_q_jw_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_q_jw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_q_jw as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_q_jw as a double




function edit_rmc_q_srf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_q_srf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_q_srf as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_q_srf as a double




% --- Executes on selection change in popup_smr_f4.
function popup_smr_f4_Callback(hObject, eventdata, handles)
% hObject    handle to popup_smr_f4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_smr_f4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_smr_f4






function edit_rmc_smr_dip_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_smr_dip as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_smr_dip as a double





function edit_rmc_smr_dipdirjoint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipdirjoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_smr_dipdirjoint as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_smr_dipdirjoint as a double





function edit_rmc_smr_dipjoint_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipjoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_smr_dipjoint as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_smr_dipjoint as a double






% --- Executes on selection change in popup_rmr_infilling.
function popup_rmr_infilling_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rmr_infilling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rmr_infilling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rmr_infilling




% --- Executes on selection change in popup_rmc_water.
function popup_rmc_water_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rmc_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rmc_water contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rmc_water



% --- Executes on selection change in popup_rmc_q_srf.
function popup_rmc_q_srf_Callback(hObject, eventdata, handles)
% hObject    handle to popup_rmc_q_srf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_rmc_q_srf contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_rmc_q_srf




function edit_rmc_jc_gsi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc_gsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jc_gsi as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jc_gsi as a double





function edit_rmc_rqd_gsi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd_gsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_rqd_gsi as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_rqd_gsi as a double





function edit_rmc_ja_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_ja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_ja as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_ja as a double





function edit_rmc_jw_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jc as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jc as a double






function edit_rmc_jl_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jl as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jl as a double




% --- Executes during object creation, after setting all properties.


function edit_rmc_blockshape_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_blockshape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_blockshape as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_blockshape as a double




function edit_rmc_waterscanline_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_waterscanline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_waterscanline as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_waterscanline as a double




function edit_height_Callback(hObject, eventdata, handles)
% hObject    handle to edit_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_height as text
%        str2double(get(hObject,'String')) returns contents of edit_height as a double



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over popup_sets_a.
function popup_sets_a_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to popup_sets_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_rmc_jr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rmc_jr as text
%        str2double(get(hObject,'String')) returns contents of edit_rmc_jr as a double


% --- Executes on button press in check_box_scan_spacing.
function checkbox61_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_spacing


% --- Executes on button press in check_box_scan_phi.
function checkbox62_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_phi


% --- Executes on button press in check_box_scan_jointsizefactor.
function checkbox63_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_jointsizefactor


% --- Executes on button press in check_box_scan_roughness.
function checkbox64_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_roughness


% --- Executes on button press in check_box_scan_alteration.
function checkbox65_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_alteration


% --- Executes on button press in check_box_sets_spacing.
function checkbox66_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_spacing


% --- Executes on button press in check_box_sets_tracelengths.
function checkbox67_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_tracelengths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_tracelengths


% --- Executes on button press in check_box_sets_phi.
function checkbox68_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_phi


% --- Executes on button press in check_box_sets_jointsizefactor.
function checkbox69_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_jointsizefactor


% --- Executes on button press in check_box_sets_roughness.
function checkbox70_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_roughness


% --- Executes on button press in check_box_sets_alteration.
function checkbox71_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_alteration


% --- Executes on button press in check_box__scan_tracelength.
function checkbox72_Callback(hObject, eventdata, handles)
% hObject    handle to check_box__scan_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box__scan_tracelength


% --- Executes on button press in pushbutton_drawgraphs.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_drawgraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in check_hist_scan_spacing.
function checkbox73_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_spacing


% --- Executes on button press in check_hist_scan_tracelength.
function checkbox74_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_tracelength


% --- Executes on button press in check_hist_scan_phi.
function checkbox75_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_phi


% --- Executes on button press in check_hist_scan_jointsizefactor.
function checkbox76_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_jointsizefactor


% --- Executes on button press in ckeck_hist_scan_roughness.
function checkbox77_Callback(hObject, eventdata, handles)
% hObject    handle to ckeck_hist_scan_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ckeck_hist_scan_roughness


% --- Executes on button press in check_hist_scan_alteration.
function checkbox78_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_scan_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_scan_alteration


% --- Executes on button press in check_hist_sets_spacing.
function checkbox79_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_spacing


% --- Executes on button press in check_hist_sets_tracelength.
function checkbox80_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_tracelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_tracelength


% --- Executes on button press in check_hist_sets_phi.
function checkbox81_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_phi


% --- Executes on button press in check_hist_sets_jointsizefactor.
function checkbox82_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_jointsizefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_jointsizefactor


% --- Executes on button press in check_hist_sets_roughness.
function checkbox83_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_roughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_roughness


% --- Executes on button press in check_hist_sets_alteration.
function checkbox84_Callback(hObject, eventdata, handles)
% hObject    handle to check_hist_sets_alteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_hist_sets_alteration


% --- Executes on button press in check_rqd_hist.
function checkbox85_Callback(hObject, eventdata, handles)
% hObject    handle to check_rqd_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_rqd_hist


% --- Executes on button press in check_rqd_distance.
function checkbox86_Callback(hObject, eventdata, handles)
% hObject    handle to check_rqd_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_rqd_distance




% --- Executes on button press in check_box_scan_term.
function check_box_scan_term_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_scan_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_scan_term


% --- Executes on button press in check_box_sets_term.
function check_box_sets_term_Callback(hObject, eventdata, handles)
% hObject    handle to check_box_sets_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_box_sets_term


function editDFN_dipdirD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipdirD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipdirD as a double



function editDFN_traceD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_traceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_traceD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_traceD as a double



function editDFN_traceE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_traceE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_traceE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_traceE as a double




function editDFN_stdA_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_stdA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_stdA as text
%        str2double(get(hObject,'String')) returns contents of editDFN_stdA as a double



function editDFN_stdB_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_stdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_stdB as text
%        str2double(get(hObject,'String')) returns contents of editDFN_stdB as a double





function editDFN_stdC_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_stdC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_stdC as text
%        str2double(get(hObject,'String')) returns contents of editDFN_stdC as a double





function editDFN_stdD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_stdD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_stdD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_stdD as a double





function editDFN_stdE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_stdE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_stdE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_stdE as a double





function editDFN_size_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_size as text
%        str2double(get(hObject,'String')) returns contents of editDFN_size as a double




% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton52.
function pushbutton52_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu16.
function popupmenu16_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu16




% --- Executes on selection change in popupmenu17.
function popupmenu17_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu17 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu17




function editDFN_dipdirA_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipdirA as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipdirA as a double



function editDFN_dipdirB_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipdirB as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipdirB as a double



function editDFN_dipdirC_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipdirC as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipdirC as a double



function editDFN_dipdirE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipdirE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipdirE as a double



function editDFN_dipA_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipA as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipA as a double



function editDFN_dipB_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipB as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipB as a double



function editDFN_dipC_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipC as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipC as a double



function editDFN_dipD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipD as a double



function editDFN_dipE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_dipE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_dipE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_dipE as a double



function editDFN_spaceA_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_spaceA as text
%        str2double(get(hObject,'String')) returns contents of editDFN_spaceA as a double





function editDFN_spaceB_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_spaceB as text
%        str2double(get(hObject,'String')) returns contents of editDFN_spaceB as a double




function editDFN_spaceC_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_spaceC as text
%        str2double(get(hObject,'String')) returns contents of editDFN_spaceC as a double




function editDFN_spaceD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_spaceD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_spaceD as a double





function editDFN_spaceE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_spaceE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_spaceE as a double





function editDFN_diaA_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_diaA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_diaA as text
%        str2double(get(hObject,'String')) returns contents of editDFN_diaA as a double




function editDFN_diaB_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_diaB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_diaB as text
%        str2double(get(hObject,'String')) returns contents of editDFN_diaB as a double





function editDFN_diaC_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_diaC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_diaC as text
%        str2double(get(hObject,'String')) returns contents of editDFN_diaC as a double






function editDFN_diaD_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_diaD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_diaD as text
%        str2double(get(hObject,'String')) returns contents of editDFN_diaD as a double





function editDFN_diaE_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_diaE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_diaE as text
%        str2double(get(hObject,'String')) returns contents of editDFN_diaE as a double




function editDFN_raster_Callback(hObject, eventdata, handles)
% hObject    handle to editDFN_raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDFN_raster as text
%        str2double(get(hObject,'String')) returns contents of editDFN_raster as a double


%% ================================================================================================================================================================================
% FUNCTIONS CALLED UPON THE CREATION OF gui ELEMENTS 
% ==================================================================================================================================================================================
% --- Executes during object creation, after setting all properties.
function edit_filelocation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filelocation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%2

function edit_sltrend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sltrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%2

function edit_slplunge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slplunge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%4


function edit_sllength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slplunge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%5

function edit_rmc_jr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sllength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%6

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%7

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%8

% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%9


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%10


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%11

% --- Executes during object creation, after setting all properties.
function edit_answer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_answer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%12

% --- Executes during object creation, after setting all properties.
function edit_setclusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_setclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%13

% --- Executes during object creation, after setting all properties.
function edit_clusterokanswer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clusterokanswer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%14

% --- Executes during object creation, after setting all properties.
function edit_clusterok_CreateFcn(hObject, eventdata, ~)
% hObject    handle to edit_clusterok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%15

% --- Executes during object creation, after setting all properties.
function edit_statusmanual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_statusmanual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%16

% --- Executes during object creation, after setting all properties.
function text_clustermeans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_clustermeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%17


% --- Executes during object creation, after setting all properties.
function popupmenu_clustersets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_clustersets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%18


% --- Executes during object creation, after setting all properties.
function popupmenu_changeclusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_changeclusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%19

% --- Executes during object creation, after setting all properties.
function popup_sets_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sets_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%20

% --- Executes during object creation, after setting all properties.
function popup_sets_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sets_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%21

% --- Executes during object creation, after setting all properties.
function popup_sets_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_sets_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%22


% --- Executes during object creation, after setting all properties.
function edit_angleab_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_angleab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%23

% --- Executes during object creation, after setting all properties.
function edit_angleac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_angleac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%24

% --- Executes during object creation, after setting all properties.
function edit_anglebc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_anglebc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%25

% --- Executes during object creation, after setting all properties.
function edit_spacing_seta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_spacing_seta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%26


% --- Executes during object creation, after setting all properties.
function edit_spacing_setb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_spacing_setb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%27

% --- Executes during object creation, after setting all properties.
function edit_spacing_setc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_spacing_setc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%28

% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%29

% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%30

% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%31


% --- Executes during object creation, after setting all properties.
function edit_blockshapefactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_blockshapefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%32

% --- Executes when uipanel34 is resized.
function uipanel34_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uipanel34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%33

% --- Executes during object creation, after setting all properties.
function edit_persistence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_persistence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%34


% --- Executes during object creation, after setting all properties.
function edit_data1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_data1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%35


% --- Executes during object creation, after setting all properties.
function edit_data2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_data2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%36

% --- Executes during object creation, after setting all properties.
function edit_data3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_data3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%37

% --- Executes during object creation, after setting all properties.
function edit_rmc_smr_rmrb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_rmrb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%38

% --- Executes during object creation, after setting all properties.
function edit_rmc_smr_dipdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%39

% --- Executes during object creation, after setting all properties.
function edit_rmc_ucs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_ucs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%40

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmi_jc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmi_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%41


% --- Executes during object creation, after setting all properties.
function edit_rmc_rmi_vb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_vb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%42

% --- Executes during object creation, after setting all properties.
function edit_rmc_jc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%43

% --- Executes during object creation, after setting all properties.
function edit_rmc_vb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_vb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%44

% --- Executes during object creation, after setting all properties.
function edit_rmc_rqd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%45


% --- Executes during object creation, after setting all properties.
function edit_rmc_jn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%46

% --- Executes during object creation, after setting all properties.
function edit_rmc_jc_q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc_q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%47

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_ucs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_ucs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%48

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_rqd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%49

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_jc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%50

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_wc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%51

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_strike_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_strike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%52

% --- Executes during object creation, after setting all properties.
function popup_rmc_rmr_wc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rmc_rmr_wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%53

% --- Executes during object creation, after setting all properties.
function edit_rmc_rmr_dip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rmr_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%54

% --- Executes during object creation, after setting all properties.
function popup_rmr_orientation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rmr_orientation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%55

% --- Executes during object creation, after setting all properties.
function edit_rmc_q_jw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_q_jw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%56

% --- Executes during object creation, after setting all properties.
function edit_rmc_q_srf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_q_srf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%57

% --- Executes during object creation, after setting all properties.
function popup_smr_f4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_smr_f4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%58

% --- Executes during object creation, after setting all properties.
function edit_rmc_smr_dip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%59

% --- Executes during object creation, after setting all properties.
function edit_rmc_smr_dipdirjoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipdirjoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%60

% --- Executes during object creation, after setting all properties.
function edit_rmc_smr_dipjoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_smr_dipjoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%61


% --- Executes during object creation, after setting all properties.
function text_betafromspacings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_betafromspacings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%62

% --- Executes during object creation, after setting all properties.
function popup_rmr_infilling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rmr_infilling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%63

% --- Executes during object creation, after setting all properties.
function popup_rmc_water_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rmc_water (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%64


% --- Executes during object creation, after setting all properties.
function popup_rmc_q_srf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_rmc_q_srf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%65

% --- Executes during object creation, after setting all properties.
function edit_rmc_jc_gsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc_gsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%66

% --- Executes during object creation, after setting all properties.
function edit_rmc_rqd_gsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_rqd_gsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%67

% --- Executes during object creation, after setting all properties.
function edit_rmc_ja_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_ja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%68

% --- Executes during object creation, after setting all properties.
function edit_rmc_jw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%69

% --- Executes during object creation, after setting all properties.
function edit_rmc_jl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_jl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%70

% --- Executes during object creation, after setting all properties.
function edit_rmc_blockshape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_blockshape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%71

% --- Executes during object creation, after setting all properties.
function edit_rmc_waterscanline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rmc_waterscanline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%72

% --- Executes during object creation, after setting all properties.
function edit_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%73

% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%74


% --- Executes during object creation, after setting all properties.
function editDFN_dipdirA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%75

% --- Executes during object creation, after setting all properties.
function editDFN_dipdirB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%76

% --- Executes during object creation, after setting all properties.
function editDFN_dipdirC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%77

% --- Executes during object creation, after setting all properties.
function edit84_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%78


% --- Executes during object creation, after setting all properties.
function edit85_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%79


% --- Executes during object creation, after setting all properties.
function edit86_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%80


%81


%82


%83


%84


%85


%86


% --- Executes during object creation, after setting all properties.
function editDFN_dipdirD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%87


% --- Executes during object creation, after setting all properties.
function editDFN_dipdirE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipdirE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%88




% --- Executes during object creation, after setting all properties.
function editDFN_dipA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%89


% --- Executes during object creation, after setting all properties.
function editDFN_dipB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%90

% --- Executes during object creation, after setting all properties.
function editDFN_dipC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%91


% --- Executes during object creation, after setting all properties.
function editDFN_dipD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%92


% --- Executes during object creation, after setting all properties.
function editDFN_dipE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_dipE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%93

% --- Executes during object creation, after setting all properties.
function editDFN_spacA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spacA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%94

% --- Executes during object creation, after setting all properties.
function editDFN_spacB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spacB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%95

% --- Executes during object creation, after setting all properties.
function editDFN_spacC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spacC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%96

% --- Executes during object creation, after setting all properties.
function editDFN_spacD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spacD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%97


% --- Executes during object creation, after setting all properties.
function editDFN_spacE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spacE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%98


% --- Executes during object creation, after setting all properties.
function editDFN_traceA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_traceA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%99


% --- Executes during object creation, after setting all properties.
function editDFN_traceB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_traceB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%100


% --- Executes during object creation, after setting all properties.
function editDFN_traceC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_traceC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%101


% --- Executes during object creation, after setting all properties.
function editDFN_traceD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_traceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%102

% --- Executes during object creation, after setting all properties.
function editDFN_traceE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_traceE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%103

% --- Executes during object creation, after setting all properties.
function editDFN_stdA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_stdA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%104

% --- Executes during object creation, after setting all properties.
function editDFN_stdB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_stdB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%105

% --- Executes during object creation, after setting all properties.
function editDFN_stdC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_stdC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%106

% --- Executes during object creation, after setting all properties.
function editDFN_stdD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_stdD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%107

% --- Executes during object creation, after setting all properties.
function editDFN_stdE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_stdE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%108

% --- Executes during object creation, after setting all properties.
function editDFN_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%109

% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%110

% --- Executes during object creation, after setting all properties.
function popupmenu17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%111

% --- Executes during object creation, after setting all properties.
function popupDFN_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupDFN_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%112

% --- Executes during object creation, after setting all properties.
function editDFN_spaceA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%113

% --- Executes during object creation, after setting all properties.
function editDFN_spaceB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%114

% --- Executes during object creation, after setting all properties.
function editDFN_spaceC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%115

% --- Executes during object creation, after setting all properties.
function editDFN_spaceD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%116

% --- Executes during object creation, after setting all properties.
function editDFN_spaceE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_spaceE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%117


% --- Executes during object creation, after setting all properties.
function editDFN_diaA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_diaA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%118

% --- Executes during object creation, after setting all properties.
function editDFN_diaB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_diaB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%119

% --- Executes during object creation, after setting all properties.
function editDFN_diaC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_diaC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%120

% --- Executes during object creation, after setting all properties.
function editDFN_diaD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_diaD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%121

% --- Executes during object creation, after setting all properties.
function editDFN_diaE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_diaE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%122


% --- Executes during object creation, after setting all properties.
function editDFN_raster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDFN_raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%123


% --- Executes during object creation, after setting all properties.
function pushbutton_clusterok_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_clusterok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




