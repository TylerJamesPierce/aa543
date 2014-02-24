%----------------------------------------------------------------------
%
% Euler1D is a code that allows the investigation of several explicit
% and implicit algorithms for the solution of the 1D Euler equations:
%
%                  dq/dt + df/dx = 0
%
% This code was developed by Uri Shumlak at the University of Washington.
% Questions can be emailed to shumlak@u.washington.edu.
%
% Last modified January 2014.
%
%----------------------------------------------------------------------

function varargout = plottype(varargin)
% PLOTTYPE M-file for plottype.fig
%      PLOTTYPE, by itself, creates a new PLOTTYPE or raises the existing
%      singleton*.
%
%      H = PLOTTYPE returns the handle to a new PLOTTYPE or the handle to
%      the existing singleton*.
%
%      PLOTTYPE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOTTYPE.M with the given input arguments.
%
%      PLOTTYPE('Property','Value',...) creates a new PLOTTYPE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plottype_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plottype_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plottype

% Last Modified by GUIDE v2.5 08-Jan-2014 11:30:30

% Begin initialization code - DO NOT EDIT
  gui_Singleton = 1;
  gui_State = struct('gui_Name',       mfilename, ...
                     'gui_Singleton',  gui_Singleton, ...
                     'gui_OpeningFcn', @plottype_OpeningFcn, ...
                     'gui_OutputFcn',  @plottype_OutputFcn, ...
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


% --- Executes just before plottype is made visible.
function plottype_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plottype (see VARARGIN)
% Choose default command line output for plottype
  handles.output = hObject;
  store=[];
  prop='Density';
  geo='Nozzle';
  algorithm='MacCormack';
  assignin('base','Algorithm',algorithm);
  assignin('base','Property',prop);
  assignin('base','savedata',store);
  assignin('base','Type',geo);
% Update handles structure
  guidata(hObject, handles);

% UIWAIT makes plottype wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = plottype_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
  varargout{1} = handles.output;

function Courant_Callback(hObject, eventdata, handles)
% hObject    handle to Courant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Courant as text
%        str2double(get(hObject,'String')) returns contents of Courant as a double

% --- Executes during object creation, after setting all properties.
function Courant_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Courant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

function tfinal_Callback(hObject, eventdata, handles)
% hObject    handle to tfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tfinal as text
%        str2double(get(hObject,'String')) returns contents of tfinal as a double

% --- Executes during object creation, after setting all properties.
function tfinal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfinal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

function imax_Callback(hObject, eventdata, handles)
% hObject    handle to imax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of imax as text
%        str2double(get(hObject,'String')) returns contents of imax as a double

% --- Executes during object creation, after setting all properties.
function imax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double

% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes on selection change in Algorithm.
function Algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to Algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns Algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Algorithm
  val = get(hObject,'Value');
  string_list = get(hObject,'String');
  selected_string = string_list{val};
  assignin('base','Algorithm',selected_string);
  guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes on selection change in wavetype.
function plottype_Callback(hObject, eventdata, handles)
% hObject    handle to wavetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns wavetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wavetype
  val = get(hObject,'Value');
  string_list = get(hObject,'String');
  selected_string = string_list{val};
  assignin('base','Wavetype',selected_string);
  guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function plottype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  guidata(hObject, handles);
  courant = str2num(get(handles.Courant,'String'))
  theta = str2num(get(handles.theta,'String'))
  imax = str2num(get(handles.imax,'String'))
  idim = imax;
  dt = 1;
  nsteps = str2num(get(handles.miter,'String'));;

  shock=evalin('base','Type');
  epse=str2num(get(handles.epse,'String'));
  epsi=str2num(get(handles.epsi,'String'));
  algorithm=evalin('base','Algorithm'); 
  prop = evalin('base','Property');
  if strcmp(prop, 'Density')==1
    set(handles.Plot_label,'String',sprintf('rho(x)'));
  elseif strcmp(prop, 'Pressure')==1
    set(handles.Plot_label,'String',sprintf('p(x)'));
  elseif strcmp(prop, 'Energy')==1
    set(handles.Plot_label,'String',sprintf('e(x)'));
  end
  gamma=1.4;
  gamma1=0.4;
  extrap='False';
  x=zeros(idim,1);
  q=zeros(idim,3);
  qbar=zeros(idim,3);
  flux=zeros(idim,3);
  s=zeros(idim,3);
  area=zeros(idim,1);
  dq=zeros(idim,3);
  uu=zeros(idim,1);
  dens=zeros(idim,1);
  amom=zeros(idim,1);
  ener=zeros(idim,1);
  press=zeros(idim,1);
  l2norm=[];

  t=0;
  tout=0;
  iter=1;
  dx=10/(imax-1);

  for i=1:imax
    x(i,1)=10*(i-1)/(imax-1);
    s(i,1)=0.0;
    s(i,2)=0.2776/(cosh(0.8*(x(i,1)-4.0))^2);
    s(i,3)=0;
    area(i,1)=1.398 + .347 * tanh(0.8*(x(i,1)-4.0));
  end

  p=1/gamma;

  if strcmp(shock,'Shock Tube')==1
    for i=1:int32(imax/2)
        s(i,2)=0;
        area(i,1)=1.0;
        q(i,1)=1.0;
        q(i,2)=0;
        q(i,3)=p/gamma1;
    end
    for i=int32(imax/2):imax
        s(i,2)=0;
        area(i,1)=1.0;
        q(i,1)=4.0;
        q(i,2)=0;
        q(i,3)=4.0*p/gamma1;
    end
  elseif strcmp(shock, 'Nozzle')==1
    q(:,1)=1.0*area(:,1);
    q(:,2)=1.26*area(:,1);
    q(:,3)=((p/gamma1)+.7938)*area(:,1);
  end
  for i=1:imax
    dens(i,1)=q(i,1)/area(i,1);
    amom(i,1)=q(i,2)/area(i,1);
    ener(i,1)=q(i,3)/area(i,1);
    press(i,1)=gamma1*(ener(i,1)-.5*(amom(i,1)*amom(i,1)/dens(i,1)));
  end
  if strcmp(prop,'Density')==1
    handles.plot = plot(x(1:imax),dens(1:imax))
  elseif strcmp(prop,'Pressure')==1
    handles.plot = plot(x(1:imax),press(1:imax))
  elseif strcmp(prop,'Energy')==1
    handles.plot = plot(x(1:imax),ener(1:imax))
  end
  while iter < nsteps
    for i=1:imax
        dens(i,1)=q(i,1)/area(i,1);
        amom(i,1)=q(i,2)/area(i,1);
        uu(i,1)=amom(i,1)/dens(i,1);
        ener(i,1)=q(i,3)/area(i,1);
        press(i,1)=gamma1*(ener(i,1)-.5*(amom(i,1)*amom(i,1)/dens(i,1)));
    end
    al2n=0.0;
    dt=1000.0;
    for i=1:imax
        u=q(i,2)/q(i,1);
        c2=gamma*gamma1*(q(i,3)-.5*q(1,1)*u*u)/q(i,1);
        c=sqrt(abs(c2));
        eig=c+abs(u);
        cfl=(dt/dx)*eig;
        dtt=courant*dx/eig;
        dt=min(dt,dtt);
    end
    t=t+dt;
    if strcmp(algorithm,'MacCormack')==1
      % Insert MacCormack Algorithm here.
      [flux]=calfx(q,flux,1,imax,gamma1);
      for i=1:imax-1
        for nx=1:3
          qbar(i,nx) = q(i,nx) - (dt/dx) * (flux(i+1,nx) - flux(i,nx));
        end
      end
      [flux]=calfx(qbar,flux,1,imax-1,gamma1);
      for i=2:imax
        for nx=1:3
          q(i,nx) = 0.5*(q(i-1,nx) + qbar(i-1,nx)) - (dt/dx)*0.5*(flux(i,nx)-flux(i-1,nx));
        end
      end
      for i=2:imax
        rhoa = q(i,1);
        rhoua = q(i,2);
        ea = q(i,3);
        u = rhoua / rhoa;
        p = gamma1 * (ea-(.5 * rhoua * u)) / area(i,1);
        dq(i,2) = dq(i,2) + (dt * p * s(i,2));
      end
    elseif strcmp(algorithm,'Lax Wendroff')==1
       [flux]=calfx(q,flux,1,imax,gamma1);
       for i=1:imax-1
           for nx=1:3
              qbar(i,nx) = 0.5*(q(i,nx) + q(i+1,nx)) - (dt/dx)*0.5*(flux(i+1,nx)-flux(i,nx));
           end
       end
       [flux]=calfx(qbar,flux,1,imax-1,gamma1);
       for i=2:imax-1
            for nx = 1:3
                 dq(i,nx) = -(dt/dx) * (flux(i,nx) - flux(i-1,nx));
            end
       end
       for i = 2:imax
            rhoa = q(i,1) + dq(i,1);
            rhoua = q(i,2) + dq(i,2);
            ea = q(i,3) + dq(i,3);
            u = rhoua / rhoa;
            p = gamma1 * (ea - (.5 * rhoua * u)) / area(i,1);
            dq(i,2) = dq(i,2) + (dt * p * s(i,2));
       end
    elseif strcmp(algorithm,'Roe')==1
        sw=zeros(idim,3);
        dalfa=zeros(idim,3);
        g=zeros(idim,3);
        a=zeros(idim,1);
        u=zeros(idim,1);
        enth=zeros(idim,1);
        qtemp=zeros(idim,3);
        arg=zeros(idim,3);
        gami = gamma - 1.0;
        diff = 1.0E-4;
        for i=1:imax-1
            rojinv=1/q(i,1);
            rojpinv=1/q(i+1,1);
            uj = q(i,2)*rojinv;
            uj1= q(i+1,2)*rojpinv;
            pj = gami*(q(i,3)-.5*q(i,1)*uj*uj);
            pj1= gami*(q(i+1,3)-.5*q(i+1,1)*uj1*uj1);
            dro= q(i+1,1)-q(i,1);
            dm = q(i+1,2)-q(i,2);
            de = q(i+1,3)-q(i,3);
            u(i)=.5*(uj+uj1);
            enth(i)=.5*((q(i,3)+pj)*rojinv+(q(i+1,3)+pj1)*rojpinv);
            a(i)=sqrt(gami*(enth(i)-.5*u(i)*u(i)));
            temp1=u(i)*dm-de-dro*(0.5*u(i)*u(i));
            dalfa(i,1) = dro + gami/(a(i)*a(i))*temp1;
            dalfa(i,2) = 0.5*(dm-u(i)*dro)/a(i);
            dalfa(i,3) =-dalfa(i,2)-0.5*gami*temp1/(a(i)*a(i));
            dalfa(i,2) = dalfa(i,2)-0.5*gami*temp1/(a(i)*a(i));
            rocuj1=q(i,2);
            rocuj2=q(i+1,2);
            qtemp(i,1)= rocuj1+rocuj2;
            qtemp(i,2)= (uj*rocuj1+pj) + (uj1*rocuj2+pj1);
            qtemp(i,3)= uj*(q(i,3)+pj) + uj1*(q(i+1,3)+pj1);
        end
        for i=1:imax-1
             e1=abs(u(i));
             e2=abs(u(i)+a(i));
             e3=abs(u(i)-a(i));
             if e1 < diff
                e1=.5*(e1*e1+diff*diff)/diff;
             end
             if e2 < diff 
                 e2=.5*(e2*e2+diff*diff)/diff;
             end
             if e3 < diff 
                 e3=.5*(e3*e3+diff*diff)/diff;
             end
             sw(i,1)= - e1*dalfa(i,1);
             sw(i,2)= - e2*dalfa(i,2);
             sw(i,3)= - e3*dalfa(i,3);
        end
        for i=1:imax-1
            temp1= sw(i,2)-sw(i,3);
            temp2= sw(i,1)+sw(i,2)+sw(i,3);
            ftemp1=temp2;
            ftemp2=u(i)*temp2 + a(i)*temp1;
            ftemp3= enth(i)*temp2 + u(i)*a(i)*temp1-a(i)*a(i)/gami*sw(i,1);
            qtemp(i,1)=.5*(qtemp(i,1)+ftemp1);
            qtemp(i,2)=.5*(qtemp(i,2)+ftemp2);
            qtemp(i,3)=.5*(qtemp(i,3)+ftemp3);
        end
        for i=2:imax-1
            dq(i,1)= - dt/dx*(qtemp(i,1)-qtemp(i-1,1));
            dq(i,2)= - dt/dx*(qtemp(i,2)-qtemp(i-1,2));
            dq(i,3)= - dt/dx*(qtemp(i,3)-qtemp(i-1,3));
            rhoa = q(i,1) + dq(i,1);
            rhoua = q(i,2) + dq(i,2);
            ea = q(i,3) + dq(i,3);
            us = rhoua / rhoa;
            p = gamma1*(ea - .5*rhoua*us)/area(i);
            dq(i,2) = dq(i,2) + dt*p*s(i,2);
        end
    elseif strcmp(algorithm,'Harten-Yee')==1
        sw=zeros(idim,3);
        dalfa=zeros(idim,3); 
        g=zeros(idim,3); 
        a=zeros(idim,1);
        u=zeros(idim,1); 
        enth=zeros(idim,1);
        qtemp=zeros(idim,3);
        arg=zeros(idim,3);

        gami = gamma - 1.0;
        diff = 1.0E-4;
% Roe's scheme - First Order
%           ----- evaluate primitive variables at nodal points -------
        for i=1:imax-1
            rojinv=1.0/q(i,1);
            rojpinv=1.0/q(i+1,1);
            uj = q(i,2)*rojinv;
            uj1= q(i+1,2)*rojpinv;
            pj = gami*(q(i,3)-.5*q(i,1)*uj*uj);
            pj1= gami*(q(i+1,3)-.5*q(i+1,1)*uj1*uj1);
%           ----- evaluate quantities at half points -----------------
            dro= q(i+1,1)-q(i,1);
            dm = q(i+1,2)-q(i,2);
            de = q(i+1,3)-q(i,3);
            u(i)=.5*(uj+uj1);
            enth(i)=.5*((q(i,3)+pj)*rojinv+(q(i+1,3)+pj1)*rojpinv);
            a(i)=sqrt(gami*(enth(i)-.5*u(i)*u(i)));
%           ----- get alfa terms -------------------------------------
            temp1=u(i)*dm-de-dro*(0.5*u(i)*u(i));
            dalfa(i,1) = dro + gami/(a(i)*a(i))*temp1;
            dalfa(i,2) = 0.5*(dm-u(i)*dro)/a(i);
            dalfa(i,3) =-dalfa(i,2)-0.5*gami*temp1/(a(i)*a(i));
            dalfa(i,2) = dalfa(i,2)-0.5*gami*temp1/(a(i)*a(i));
%           ----- get twice the real flux term ------------------------
            rocuj1=q(i,2);
            rocuj2=q(i+1,2);
            qtemp(i,1)= rocuj1+rocuj2;
            qtemp(i,2)= (uj*rocuj1+pj) + (uj1*rocuj2+pj1);
            qtemp(i,3)= uj*(q(i,3)+pj) + uj1*(q(i+1,3)+pj1);
        end
%        ----- Second order TVD dissipation for internal nodes -------
%            ----- calculate sigma and put them in arg(i) ------------
        for i=1:imax-1
             e1=abs(u(i));
             e2=abs(u(i)+a(i));
             e3=abs(u(i)-a(i));
             if e1 < diff
                 arg(i,1)=.25*(e1*e1+diff*diff)/diff;
             else
                 arg(i,1)=.5*e1;
             end
             if e2 < diff
                 arg(i,2)=.25*(e2*e2+diff*diff)/diff;
             else
                 arg(i,2)=.5*e2;
             end
             if e3 < diff
                 arg(i,3)=.25*(e3*e3+diff*diff)/diff;
             else
                 arg(i,3)=.5*e3;
             end
        end
        for i=1:imax-1
              if dalfa(i,1) <= 0
                  sw(i,1)=-1;
              elseif dalfa(i,1) > 0
                sw(i,1)=1.0;
              end
              if dalfa(i,2) <= 0
               sw(i,2)=-1.0;
              elseif dalfa(i,2) > 0
               sw(i,2)=1.0;
              end
              if dalfa(i,3) <= 0
               sw(i,3)=-1.0;
              elseif dalfa(i,3) > 0
               sw(i,3)=1.0;
              end
        end
         g(1,1)=0.0;
         g(1,2)=0.0;
         g(1,3)=0.0;
         g(imax,1)=0.0;
         g(imax,2)=0.0;
         g(imax,3)=0.0;
         for i=2:imax-1
             g(i,1)=sw(i,1)*max(0.0,min(arg(i,1)*abs(dalfa(i,1)),sw(i,1)*arg(i-1,1)*dalfa(i-1,1)));
             g(i,2)=sw(i,2)*max(0.0,min(arg(i,2)*abs(dalfa(i,2)),sw(i,2)*arg(i-1,2)*dalfa(i-1,2)));
             g(i,3)=sw(i,3)*max(0.0,min(arg(i,3)*abs(dalfa(i,3)),sw(i,3)*arg(i-1,3)*dalfa(i-1,3)));
         end
%        ---- store dgamma into arg(i,m) -----------------------------
         for i=1:imax-1
             if dalfa(i,1)==0
                arg(i,1)=0.0;
             else
                arg(i,1)=(g(i+1,1)-g(i,1))/dalfa(i,1);
             end
             if dalfa(i,2)==0
                arg(i,2)=0.0;
             else
                arg(i,2)=(g(i+1,2)-g(i,2))/dalfa(i,2);
             end
             if dalfa(i,3)==0
                arg(i,3)=0.0;
             else
                arg(i,3)=(g(i+1,3)-g(i,3))/dalfa(i,3);
             end
         end
%        ---- store phi(m) into sw(i,m) -------------------------------
         for i=1:imax-1
             e1=abs(u(i)+arg(i,1));
             e2=abs(u(i)+a(i)+arg(i,2));
             e3=abs(u(i)-a(i)+arg(i,3));
             if e1 < diff
                 e1=.5*(e1*e1+diff*diff)/diff;
             end
             if e2 < diff 
                 e2=.5*(e2*e2+diff*diff)/diff;
             end
             if e3 < diff 
                 e3=.5*(e3*e3+diff*diff)/diff;
             end
             sw(i,1)=g(i,1)+g(i+1,1) - e1*dalfa(i,1);
             sw(i,2)=g(i,2)+g(i+1,2) - e2*dalfa(i,2);
             sw(i,3)=g(i,3)+g(i+1,3) - e3*dalfa(i,3);
         end
%        ----- correct the flux using TVD dissipation term --------
         for i=1:imax-1
              temp1= sw(i,2)-sw(i,3);
              temp2= sw(i,1)+sw(i,2)+sw(i,3);
              ftemp1=temp2;
              ftemp2=u(i)*temp2 + a(i)*temp1;
              ftemp3= enth(i)*temp2 + u(i)*a(i)*temp1-a(i)*a(i)/gami*sw(i,1);
              qtemp(i,1)=.5*(qtemp(i,1)+ftemp1);
              qtemp(i,2)=.5*(qtemp(i,2)+ftemp2);
              qtemp(i,3)=.5*(qtemp(i,3)+ftemp3);
         end
%        ------ store residual in i-direction   --------------------	 
         for i=2:imax-1
            dq(i,1)= - dt/dx*(qtemp(i,1)-qtemp(i-1,1));
            dq(i,2)= - dt/dx*(qtemp(i,2)-qtemp(i-1,2));
            dq(i,3)= - dt/dx*(qtemp(i,3)-qtemp(i-1,3));
            rhoa = q(i,1) + dq(i,1);
            rhoua = q(i,2) + dq(i,2);
            ea = q(i,3) + dq(i,3);
            us = rhoua / rhoa;
            p = gamma1 * (ea - .5 * rhoua * us) / area(i);
            dq(i,2) = dq(i,2) + dt * p * s(i,2);
         end
    elseif strcmp(algorithm,'Theta')==1
      % Insert Theta-Method Algorithm here.
    end
    al2n=0.0;
    for i=2:imax-1
        q(i,1)=q(i,1) + dq(i,1);
        q(i,2)=q(i,2) + dq(i,2);
        q(i,3)=q(i,3) + dq(i,3);
        al2n=al2n + dq(i,1)*dq(i,1);
    end
    l2norm(iter)=sqrt(al2n);
    if strcmp(shock,'Nozzle')==1
        [q]=bc(q,area,imax,extrap,gamma,gamma1);
    end
    if t >= tout
        tout=tout+dt;
    end
    if strcmp(prop,'Density')==1
    plot(x(1:imax),dens(1:imax))
    elseif strcmp(prop,'Pressure')==1
    plot(x(1:imax),press(1:imax))
    elseif strcmp(prop,'Energy')==1
    plot(x(1:imax),ener(1:imax))
    end
    pause(0.025)
    set(handles.Iter,'String',sprintf('Iterations: %g',iter));
    iter=iter+1;
  %end
    data=horzcat(x, dens, press, ener);
    store = mat2str(data);
    assignin('base','savedata',store);
  end

% --- Executes on button press in datasave.
function datasave_Callback(hObject, eventdata, handles)
% hObject    handle to datasave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  data=evalin('base','savedata');
  [a b]=size(data)
  if a==0 && b==0
    error('No data to save')
  end
  data=eval(data);
  dlmwrite('q_output_data.dat',data,'precision','%14.6e')

function miter_Callback(hObject, eventdata, handles)
% hObject    handle to miter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of miter as text
%        str2double(get(hObject,'String')) returns contents of miter as a double

% --- Executes during object creation, after setting all properties.
function miter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to miter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes on selection change in Boundary.
function geotype_Callback(hObject, eventdata, handles)
% hObject    handle to Boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns Boundary contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Boundary
  val = get(hObject,'Value');
  string_list = get(hObject,'String');
  selected_string = string_list{val};
  assignin('base','Type',selected_string);
  guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function geotype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Boundary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5

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

% --- Executes on selection change in property.
function property_Callback(hObject, eventdata, handles)
% hObject    handle to property (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = get(hObject,'String') returns property contents as cell array
%        contents{get(hObject,'Value')} returns selected item from property
  val = get(hObject,'Value');
  string_list = get(hObject,'String');
  selected_string = string_list{val};
  assignin('base','Property',selected_string);
  guidata(hObject, handles);
  if strcmp(selected_string, 'Density')==1
    set(handles.Plot_label,'String',sprintf('Density vs. X plot'));
  elseif strcmp(selected_string, 'Pressure')==1
    set(handles.Plot_label,'String',sprintf('Pressure vs. X plot'));
  elseif strcmp(selected_string, 'Energy')==1
    set(handles.Plot_label,'String',sprintf('Energy vs. X plot'));
  end

% --- Executes during object creation, after setting all properties.
function property_CreateFcn(hObject, eventdata, handles)
% hObject    handle to property (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

function epse_Callback(hObject, eventdata, handles)
% hObject    handle to epse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of epse as text
%        str2double(get(hObject,'String')) returns contents of epse as a double

% --- Executes during object creation, after setting all properties.
function epse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

function epsi_Callback(hObject, eventdata, handles)
% hObject    handle to epsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of epsi as text
%        str2double(get(hObject,'String')) returns contents of epsi as a double

% --- Executes during object creation, after setting all properties.
function epsi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
  end

% --- Executes during object creation, after setting all properties.
function plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate plot

% --- Executes on button press in Pause.
function Pause_Callback(hObject, eventdata, handles)
% hObject    handle to Pause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  pause(5);
