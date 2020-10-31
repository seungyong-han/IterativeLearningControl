function [A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib(ex);
%
%  COMPleib: COnstrained Matrix-optimization Problem library 
%
%   Version: 1.1
%      Date: 09.03.2005 (last change)
%    Author: Friedemann Leibfritz
%
%
%     References: 
%     [1] COMPl{_e}ib: COnstrained Matrix-optimization Problem library --
%         a collection of test examples for nonlinear semidefinite
%         programs, control system design and related problems
%         F. Leibfritz, Tech.-Report, 2003
%
%     [2] COMPl{_e}ib 1.0 -- User manual and quick reference
%         F. Leibfritz and W. Lipinski, Tech.-Report, 2004
%
%     [3] Description of the benchmark examples in COMPl{_e}ib 1.0
%         F. Leibfritz and W. Lipinski, Tech.-Report, 2003
%
%
%  Input: ex -- string variable which refers to an example in COMPlib,
%                i.e. 'AC1'. For more details, i.e., see [2], [3]
%
%  Output: data matrices
%            A -- real (nx times nx) matrix
%            B -- real (nx times nu) matrix
%            C -- real (ny times nx) matrix
%           B1 -- real (nx times nw) matrix
%           C1 -- real (nz times nx) matrix
%          D11 -- real (nz times nw) matrix
%          D12 -- real (nz times nu) matrix
%          D21 -- real (ny times nw) matrix
%           nx -- dimension of the 'state'
%           nu -- dimension of the 'control'
%           ny -- dimension of the 'measurement'
%           nw -- dimension of the 'noise'
%           nz -- dimension of the 'regulated output'
%
%  Note: From the output data of COMPleib, one can define a control
%        system of the form
%
%        (d/dt)x(t) =  A*x(t)+ B1*w(t)+ B*u(t);  (x - state, u - control, w- noise)
%              z(t) = C1*x(t)+D11*w(t)+D12*u(t); (z - reg. output)
%              y(t) =  C*x(t)+D21*w(t).          (y - measurement)
%
%        Depending on the control design goals, from the COMPleib data
%        it is possible to form several constraint matrix optimization
%        problems, i.e.
%          --- nonlinear semidefinite programs (NSDPs)
%          --- or, equivalently, bilinear matrix inequality (BMI) problems 
%          --- linear semidefinite programs.
%        For more details, see [1], [2].
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[];B1=[];B=[];C1=[];C=[];D11=[];D12=[];D21=[];

%------------------------------------------------------------------
% (AC1)/(AC2): Y. S. Hung and A. G. J. MacFarlane, "Multivariable
%              feedback: A quasi--classical approach", Springer-Verlag,
%              "Lecture Notes in Control and Information Sciences",
%              1982  
%------------------------------------------------------------------
if(strmatch('AC1',ex,'exact'))
 nx=5;nu=3;ny=3; 
 A=[0,0,1.132,0,-1; 0,-0.0538,-0.1712,0,0.0705;
    0,0,0,1,0; 0,0.0485,0,-0.8556,-1.013; 0,-0.2909,0,1.0532,-0.6859];
 B=[0,0,0;-0.12,1,0;0,0,0;4.419,0,-1.665;1.575,0,-0.0732];
 C=[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0];
 B1=[3.593e-2, 0, 1.672e-2; 0, 9.89e-3, 0; 0, -7.548e-2, 0; 
     0, 0, 5.635e-2;  1.45e-3,0,6.743e-2];
 C1=(1/sqrt(2))*C(2:3,:); 
 D12=(1/sqrt(2))*[1,0,0;0,1,0];
 [nx,nw]=size(B1); [nz,nx]=size(C1); [nz,nu]=size(D12); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (AC2): like (AC1) with changed C1, D11, D12
%------------------------------------------------------------------
if(strmatch('AC2',ex,'exact'))
 nx=5;nu=3;ny=3; 
 A=[0,0,1.132,0,-1; 0,-0.0538,-0.1712,0,0.0705;
    0,0,0,1,0; 0,0.0485,0,-0.8556,-1.013; 0,-0.2909,0,1.0532,-0.6859];
 B=[0,0,0;-0.12,1,0;0,0,0;4.419,0,-1.665;1.575,0,-0.0732];
 C=[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0];
 B1=[3.593e-2, 0, 1.672e-2; 0, 9.89e-3, 0; 0, -7.548e-2, 0; 
     0, 0, 5.635e-2;  1.45e-3,0,6.743e-2]; 
 C1=(1/sqrt(2))*[0,1,0,0,0;0,0,1,0,0;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0];
 D12=(1/sqrt(2))*[0,0,0;0,0,0;1,0,0;0,1,0;0,0,1];
 [nx,nw]=size(B1); [nz,nx]=size(C1); [nz,nu]=size(D12); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
%(AC3): L-1011 aircraft in cruise flight conditions
%       C. Edwards and S. K. Spurgeon,
%       On the development of discontinuous observers",
%       IJOC, Vol. 59, Nr. 5, pp. 1211-1229, 1994
%------------------------------------------------------------------
if(strmatch('AC3',ex,'exact'))
 nx=5;nu=2;ny=4;
 A=[0 0 1 0 0;0 -0.154 -0.0042 1.54 0;0 0.249 -1 -5.2 0;
    0.0386 -0.996 -0.0003 -0.117 0;0 0.5 0 0 -0.5];
 B=[0 0;-0.744 -0.032;0.337 -1.12;0.02 0;0 0];
 C=[0 1  0 0 -1;0 0 1 0 0;0 0 0 1 0;1 0 0 0 0];
end

%------------------------------------------------------------------
% (AC4): Missile autopilot;  %%ehemals (MA1)
%        B. Fares, P. Apkarian and D. Noll,
%        "An Augmented Lagrangian Method for a Class of LMI-Constrained
%        Problems in Robust Control Theory",
%        IJOC, Vol. 74, Nr. 4, pp. 348-360
%------------------------------------------------------------------
if(strmatch('AC4',ex,'exact')) 
 A=[-0.876 1 -0.1209 0; 8.9117 0 -130.75 0; 0 0 -150 0; -1 0 0 -0.05];
 B=[0;0;150;0];
 C=[-1 0 0 0;0 -1 0 0];
 B1=[0 0;0 0; 0 0; 0 1];
 C1=[-0.25 0 0 3.487;0 0 -3 0];
 D12=[0; 3]; D11=[0 0.25;0 0]; D21=[0 1;0.01 0];
end

%------------------------------------------------------------------
%(AC5): Boeing B-747 aircraft
%       T. Ishihara, H.-J. Guo and H. Takeda,
%       "A Design of Discrete-Time Integral Controllers with
%       Computation Delays via Loop Transfer Recovery",
%       AUTO, Vol. 28, Nr. 3, pp. 599-603, 1992
%------------------------------------------------------------------
if(strmatch('AC5',ex,'exact'))
 nx=4;nu=2;ny=2;
 A=[0.9801  0.0003 -0.0980  0.0038;
   -0.3868  0.9071  0.0471 -0.0008;
    0.1591 -0.0015  0.9691  0.0003;
   -0.0198  0.0958  0.0021  1      ];
B=[-0.0001  0.0058; 0.0296  0.0153; 0.0012 -0.0908; 0.0015  0.0008];
C=[1 0 0 0;0 0 0 1];
end

%------------------------------------------------------------------
% (AC6): THE AIRCRAFT L- 1011 MODEL ???
%------------------------------------------------------------------
if(strmatch('AC6',ex,'exact'))
 nx=7;nu=2;ny=4;
 A=[0 0 1 0 0 0 0;
    0 -0.154 -0.0042 1.54 0 -0.744 -0.032;
    0 0.249 -1 -5.2 0 0.337 -1.12;
    0.0386 -0.996 -0.0003 -2.117 0 0.02 0;
    0 0.5 0 0 -4 0 0;
    0 0 0 0 0 -20 0;
    0 0 0 0 0 0 -25];
 B=[0 0;0 0;0 0;0 0;0 0;20 0;0 25];
 C=[0 -0.154 -0.0042 1.54 0 -0.744 -0.032;
    0 0.249 -1 -5.2 0 0.337 -1.12;
    1 0 0 0 0 0 0;
    0 0 0 0 1 0 0];
end

%------------------------------------------------------------------
% (AC7): Transport Aircraft model (Boing flight condition VMIN)
%        D. Gangsaas, K. R. Bruce, J. D. Blight and U.-L. Ly,
%        "Application of Modern Synthesis to Aircraft Control:
%        Three Case Studies", TOAC, Vol.31, Nr.11, pp.995-1014, 1986
%        Case study III 2)
%------------------------------------------------------------------
if(strmatch('AC7',ex,'exact'))
 nx=9;nu=1;ny=2;
 A=[-0.06254, 0.01888,      0,-0.56141,-0.02751,  0, 0.06254,-0.00123, 0;
     0.01089,-0.99280, 0.99795, 0.00097,-0.07057,  0,-0.01089, 0.06449, 0;
     0.07743, 1.67540,-1.31111,-0.00030,-4.25030,  0,-0.07743,-0.10883, 0;
           0,       0,       1,       0,       0,  0,       0,       0, 0;
           0,       0,       0,       0,-20.0000, 20,       0,       0, 0;
           0,       0,       0,       0,       0,-30,       0,       0, 0;
           0,       0,       0,       0,       0,  0,-0.88206,       0, 0;
           0,       0,       0,       0,       0,  0,       0,-0.88206, 0.00882;
           0,       0,       0,       0,       0,  0,       0,-0.00882,-0.88206];
 B=[0;0;0;0;0;30;0;0;0];
 C=[-0.00519,0.47604, 0.00098, -0.00031, 0.03378, 0,0.00519,-0.03086, 0;
            0,      0,       1,        0,       0, 0,      0,       0, 0];
 B1=[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0;
     1.3282,0,0,0; 0,1.62671,0,0;0,-68.75283,0,0];
 C1=(1/sqrt(2))*C(1:1,:);
 D12=(1/sqrt(2))*[1];
 D21=[0 0 0 0; 0 0 0 1];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (AC8): Transport Aircraft model (Boing flight condition CRUISE)
%        (see (AC7)!)
%        Case study II, p.1001/1012
%------------------------------------------------------------------
if(strmatch('AC8',ex,'exact'))
 nx=9;nu=1;ny=5;
 A=[-0.01365,0.17800, 0.00017,-0.56100,-0.03726,    0, 0.01365,-0.01311,       0;
   -0.01516,-0.75200, 1.00100, 0.00127,-0.06311,    0, 0.01516, 0.05536,       0;
    0.00107, 0.07896,-0.87250,       0,-3.39900,    0,-0.00107,-0.00581,       0;
          0,       0,       1,       0,       0,    0,       0,       0,       0;
          0,       0,       0,       0,-20.0000,10.72,       0,       0,       0;
          0,       0,       0,       0,       0,  -50,       0,       0,       0;
          0,       0,       0,       0,       0,    0,-0.44470,       0,       0;
          0,       0,       0,       0,       0,    0,       0,-0.44470, 0.00440;
          0,       0,       0,       0,       0,    0,       0,-0.00440,-0.44470];
 B=[0;0;0;0;0;50;0;0;0];
 C=[0.00646,0.3203,-0.03358,     0,-0.10320, 0,-0.00646,-0.02358, 0;
          1,     0,       0,     0,       0, 0,      -1,       0, 0;
   -0.01365,0.1780, 0.00017,-0.561,-0.03726, 0, 0.01365,-0.01311, 0;
          0,-13.58,       0, 13.58,       0, 0,       0,       0, 0;
          0,     0,       1,     0,       0, 0,       0,       0, 0];
 C1=(1/2)*C(1:2,:);
 B1=[0,     0,0,1, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0,50,0,0,0,0,0;
0.9431,     0,0,0, 0,0,0,0,0,0;
     0, 1.155,0,0, 0,0,0,0,0,0;
     0,-48.82,0,0, 0,0,0,0,0,0];
  D12=(1/2)*[1;1];
  D21=[0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 1 0];
  [nx,nw]=size(B1); [nz,nx]=size(C1);
  D11=zeros(nz,nw);    
end

%------------------------------------------------------------------
% (AC9): Transport Aircraft model (Boing flight condition CRUISE); %%ehemals (AC12)
%        (see (AC7)!)     
%------------------------------------------------------------------
if(strmatch('AC9',ex,'exact'))
 nx=10;
 A=[-0.01365,0.17800, 0.00017,-0.56100,-0.03726,    0, 0.01365,-0.01311,       0 -1;
   -0.01516,-0.75200, 1.00100, 0.00127,-0.06311,    0, 0.01516, 0.05536,       0  0 ;
    0.00107, 0.07896,-0.87250,       0,-3.39900,    0,-0.00107,-0.00581,       0  0;
          0,       0,       1,       0,       0,    0,       0,       0,       0  0;
          0,       0,       0,       0,-20.0000,10.72,       0,       0,       0  0;
          0,       0,       0,       0,       0,  -50,       0,       0,       0  0;
          0,       0,       0,       0,       0,    0,-0.44470,       0,       0  0;
          0,       0,       0,       0,       0,    0,       0,-0.44470, 0.00440  0;
          0,       0,       0,       0,       0,    0,       0,-0.00440,-0.44470  0;
	  0,       0,       0,       0,       0,    0,       0,       0        0  0];
 B=[0;0;0;0;0;50;0;0;0;0];
 C=[0.00646,0.3203,-0.03358,     0,-0.10320, 0,-0.00646,-0.02358, 0 0;
          1,     0,       0,     0,       0, 0,      -1,       0,  0 0;
   -0.01365,0.1780, 0.00017,-0.561,-0.03726, 0, 0.01365,-0.01311,  0 0;
          0,-13.58,       0, 13.58,       0, 0,       0,       0,  0 0;
          0,     0,       1,     0,       0, 0,       0,       0,  0 0];
 C1=C(1:2,:);
 B1=[0,     0,0,1, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0, 0,0,0,0,0,0;
     0,     0,0,0,50,0,0,0,0,0;
0.9431,     0,0,0, 0,0,0,0,0,0;
     0, 1.155,0,0, 0,0,0,0,0,0;
     0,-48.82,0,0, 0,0,0,0,0,0;
     0,     0,1,0, 0,0,0,0,0,0];
  B=[B B1(:,1:3)];   
  [nz,nx]=size(C1);  [nx,nu]=size(B); 
  D12=ones(nz,nu);
  D21=[0 0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 0 0 0 1;
      0 0 0 0 0 0 0 0 1 0];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (AC10): B767 aircraft at a flutter condition
%         E. J. Davison, "Benchmark Problems for Control System Design",
%         "Report of the IFAC Theory Comittee", 1990 
%------------------------------------------------------------------
if(strmatch('AC10',ex,'exact'))
 nx=55;nu=2;ny=2;
 load ac10; % .mat file containing A, B, C, C1, B1, D12, D21
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (AC11): CCV-type aircraft; 
%         A. T. Alexandridis and P. N. Paraskevopoulos, "A New Approach
%         to Eigenstructure Assignment by Output Feedback",
%         TOAC, Vol. 41, Nr. 7, pp. 1046-1050, 1996; Example 2 
%------------------------------------------------------------------
if(strmatch('AC11',ex,'exact'))
 A=[-1.3410  0.9933 0 -0.1689 -0.2518;
    43.2230 -0.8693 0 -17.251 -1.5766;
    1.34100  0.0067 0  0.1689  0.2518;
    0 0 0 -20 0;0 0 0 0 -20];
 B=[0 0;0 0;0 0;20 0;0 20];
 C=[0 0 1 0 0; 47.76 -0.268 0 -4.56 4.45;
    0 0 0 1 0; 0 0 0 0 1];
end

%------------------------------------------------------------------
% (AC12): ASTOVL Aircraft    %%ehemals (AC9)
%         S. Toffner-Clausen, "System Identification and Robust Control:
%         A Case Study Approach", Springer-Verlag,
%         "Advances in Industrial Conrol", 1996; p. 274         
%------------------------------------------------------------------ 
if(strmatch('AC12',ex,'exact'))
 nx=4;nu=3;ny=4;
 A=[-0.0017, 0.0413,-5.3257,-9.7565;
    -0.0721,-0.3393,49.5146,-1.0097;
    -0.0008, 0.0138,-0.2032, 0.0009;
          0,      0,      1,      0];
 B=[ 0.2086,-0.0005,-0.0271;
     -0.0005, 0.2046, 0.0139;
     -0.0047, 0.0023, 0.1226;
          0,      0,      0];
 C=[      0,      0,57.2958,      0;
           0,      0,      0,57.2958;
      0.1045,-0.9945, 0.1375,51.5791;
     -0.0002, 0.0045,      0,      0];
 B1=[ 0.033,     0,     0;
         0, 0.048,-0.002;
    -0.064,     0, 0.340;
         0,     0, 0.006];    
 C1=(1/sqrt(2))*[1, 0, 0, 0];
 D12=(1/sqrt(2))*[0,0,1];
 D21=[0,0,0;0,0,0;0,0,0;0.0212,0,0]; % (ny,nw) 
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (AC13): ASTOVL Aircraft,  LARGE model
%         save ASTOVL_aug A B1 B2 C1 C2 D11 D12 D21 D22 Ag Bg Cg Dg
%            in ASTOVL.m /export/home/leibfr/Lipinski/matlab
%------------------------------------------------------------------
if(strmatch('AC13',ex,'exact'))
 nx=28;nu=3;ny=4;
 load ac13_14; %% ehemals "ASTOVL_aug.mat"
 % .mat file containing Ag, Bg, Cg -> A, B, C
 A=Ag; B=Bg; C=Cg; B1=[]; C1=[]; D11=[]; D12=[]; D21=[]; 
end

%------------------------------------------------------------------
% (AC14): ASTOVL Aircraft,  LARGE model
%         save ASTOVL_aug A B1 B2 C1 C2 D11 D12 D21 D22 Ag Bg Cg Dg
%            in ASTOVL.m /export/home/leibfr/Lipinski/matlab
%------------------------------------------------------------------ 
if(strmatch('AC14',ex,'exact'))
 nx=40;nu=3;ny=4;
 load ac13_14; %% ehemals "ASTOVL_aug.mat"
 % .mat file containing A, B, C, C1, B1, D12, D21
 A=A; B=B2; C=C2;	  
end

%--------------------------------------------------------------------
% (AC15) Mach 2.7 flight condition of a supersonic transport aircraft  %%ehemals(NN2) 	 
%        "Computation of Optimal Output Feedback Gains for Linear
%        Multivariable Systems", TOAC, Vol. 19, pp. 257-258, 1974
%-------------------------------------------------------------------- 
if(strmatch('AC15',ex,'exact'))
 A=[-.037 .0123 .00055 -1;0 0 1 0;-6.37 0 -.23 .0618;1.25 0 .016 -.0457];
 B=[.00084 .000236;0 0;.08 .804;-.0862 -.0665];
 C=[0 1 0 0;0 0 1 0;0 0 0 1];
 B1=eye(size(A)); 
 C1=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
 D12=[0,0;0,0;0,0;0,0;1,0;0,1];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%--------------------------------------------------------------------
% (AC16)  like (AC15) with changed C and D21  %%ehemals(NN3) 	 
%-------------------------------------------------------------------- 
if(strmatch('AC16',ex,'exact'))
 A=[-.037 .0123 .00055 -1;0 0 1 0;-6.37 0 -.23 .0618;1.25 0 .016 -.0457];
 B=[.00084 .000236;0 0;.08 .804;-.0862 -.0665];
 C=eye(size(A));
 B1=eye(size(A)); 
 C1=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
 D12=[0,0;0,0;0,0;0,0;1,0;0,1];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (AC17): Leteral axis dynamic for a L-1011 aircraft
%         A. R. Galimidi and B. R. Bramish
%        "The constrained Lyapunov problem and its application
%          to robust output feedback stabilization"
%         TOAC Vol. 31,Nr. 5, pp.410-419, 1986
%------------------------------------------------------------------
if(strmatch('AC17',ex,'exact'))
 nx=4;nu=1;ny=2;
 A=[-2.98 , 0.93  , 0     ,-0.034 ;
    -0.99 ,-0.21  , 0.035 ,-0.0011;
     0    , 0     , 0     , 1     ;
     0.39 ,-5.555 , 0     ,-1.89   ];
 B=[-0.032; 0; 0;-1.6];
 C=[ 0 , 0 , 1 , 0;
     0 , 0 , 0 , 1 ];
end

%------------------------------------------------------------------
% (AC18): B767 aircraft at a flutter condition 
%         Davison (see (AC10)!)
%         reduced order system generated by
%            /export/home/leibfr/bsp37/bsp37bal.m 
%------------------------------------------------------------------
if(strmatch('AC18',ex,'exact'))
 n=10;nu=2;ny=2;
 load ac18; % .mat file containing A, B, C, C1, B1, D12, D21
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (HE1): Longitudinal motion of a VTOL helicopter
%        S. N. Singh and A. A. R. Coelho,
%        "Nonlinear control of mismatched uncertain linear systems
%        and application to control of aircraft",
%        Journal of Dynamic Systems, Measurement and Control, Vol. 106,
%        pp. 203-210, 1984
%------------------------------------------------------------------
if(strmatch('HE1',ex,'exact'))
 nx=4;nu=2;ny=1;
 A=[-0.0366 0.0271 0.0188 -0.4555;0.0482 -1.01 0.0024 -4.0208;
    0.1002 0.3681 -0.7070 1.42;0 0 1 0];
 B=[0.4422 0.1761;3.5446 -7.5922;-5.52 4.49;0 0];
 C=[0 1 0 0];
 B1=[4.678e-2, 0; 4.572e-2, 9.88e-3; 4.369e-2, 1.11e-3; -2.179e-2, 0];
 C1=(1/sqrt(2))*[2,0,0,0;0,1,0,0];
 D12=(1/sqrt(2))*[1, 0; 0, 1];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (HE2): AH-64  HELICOPTER at 130 knots
%        Ph. M. Fitzsimons, "Reducing the computation required to solve
%        a standard minimax problem", AUTO, Vol.31, pp.1885-1887, 1995
%------------------------------------------------------------------
if(strmatch('HE2',ex,'exact'))
 nx=4;nu=2;ny=2;
 A=[-0.0649 0.0787 0.1705 -0.5616;0.0386 -0.939 4.2277 0.0198;
     0.1121 -0.4254 -0.7968 0;0 0 1 0];
 B=[-0.9454 0.5313;-8.6476 -10.769;19.0824 -2.8959;0 0];
 C=[1 0 0 0;0 0 0 1];
end

%------------------------------------------------------------------
% (HE3): Bell 201A-1 helicopter
%        D.-W. Gu, P. Hr. Petkov and M. M. Konstantinov,
%        "H_inf and H_2 Optimization Toolbox in SLICOT",
%        SLICOT Working Note 1999-12,
%        available via ftp: wgs.esat.kuleuven.ac.be/ 
%                           pub/WGS/REPORTS/SLWN1999-12.ps.Z
%------------------------------------------------------------------
if(strmatch('HE3',ex,'exact'))
 nx=8;nu=4;ny=6;
 A=[-0.0046 , 0.038  , 0.3259 ,-0.0045 ,-0.402  ,-0.073  ,-9.81 , 0   ;
    -0.1978 ,-0.5667 , 0.357  ,-0.0378 ,-0.2149 , 0.5683 , 0    , 0   ;
     0.0039 ,-0.0029 ,-0.2947 , 0.007  , 0.2266 , 0.0148 , 0    , 0   ;
     0.0133 ,-0.0014 ,-0.4076 ,-0.0654 ,-0.4093 , 0.2674 , 0    , 9.81;
     0.0127 ,-0.01   ,-0.8152 ,-0.0397 ,-0.821  , 0.1442 , 0    , 0   ;
    -0.0285 ,-0.0232 , 0.1064 , 0.0709 ,-0.2786 ,-0.7396 , 0    , 0   ;
     0      , 0      , 1      , 0      , 0      , 0      , 0    , 0   ;
     0      , 0      , 0      , 0      , 1      , 0      , 0    , 0   ];
 B=[ 0.0676 , 0.1221 ,-0.0001 ,-0.0016;
    -1.1151 , 0.1055 , 0.0039 , 0.0035;
     0.0062 ,-0.0682 , 0.001  ,-0.0035;
    -0.017  , 0.0049 , 0.1067 , 0.1692;
    -0.0129 , 0.0106 , 0.2227 , 0.143 ;
     0.139  , 0.0059 , 0.0326 ,-0.407 ;
     0      , 0      , 0      , 0     ; 
     0      , 0      , 0      , 0      ];
 C=[ 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0;
     0 , 0 , 0 , 0 , 0 , 0 , 1 , 0;
     0 , 0 , 0 , 0 , 0 , 0 , 0 , 1;
     0 , 0 , 0 , 0 , 0 , 1 , 0 , 0;
     0 , 0 , 1 , 0 , 0 , 0 , 0 , 0;
     0 , 0 , 0 , 0 , 1 , 0 , 0 , 0];
 C1=[C;zeros(4,8)];
 D12=[zeros(6,4);eye(4)];
 B1=B(:,1);
 D21=0.1*[0; 1; 0; 0; 0.5; 0];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (HE4): Helicopter control
%        "Multivariable feedback control: Analysis and design"
%        S. Skogestad and I. Postlethwaite
%        John Wiley and Sons, 1996, Section 12.2
%        Note: Matlab files 
%              http://www.nt.ntnu.no/users/skoge/book/matlab.html
%        stored in  /export/home/leibfr/Lipinski/matlab/..
%                                ..Examples_Multi_Feedback_Control/matlab_m/
%        F. Leibfritz, 16.09.2003
%        Data matrices unscaled in Sec12_2.m (in directory above)
%        cf. page 472
%------------------------------------------------------------------
if(strmatch('HE4',ex,'exact'))
 nx=8;nu=4;ny=6;
 a01 = [          0                  0                  0   0.99857378005981;
                  0                  0   1.00000000000000  -0.00318221934140;
                  0                  0 -11.57049560546880  -2.54463768005371;
                  0                  0   0.43935656547546  -1.99818229675293;
                  0                  0  -2.04089546203613  -0.45899915695190;
 -32.10360717773440                  0  -0.50335502624512   2.29785919189453;
   0.10216116905212  32.05783081054690  -2.34721755981445  -0.50361156463623;
  -1.91097259521484   1.71382904052734  -0.00400543212891  -0.05741119384766];

 a02 = [0.05338427424431             0                  0                  0;
   0.05952465534210                  0                  0                  0;
  -0.06360262632370   0.10678052902222  -0.09491866827011   0.00710757449269;
                  0   0.01665188372135   0.01846204698086  -0.00118747074157;
  -0.73502779006958   0.01925575733185  -0.00459562242031   0.00212036073208;
                  0  -0.02121581137180  -0.02116791903973   0.01581159234047;
   0.83494758605957   0.02122657001019  -0.03787973523140   0.00035400385968;
                  0   0.01398963481188  -0.00090675335377  -0.29051351547241];
 
 A=[a01 a02];
 B1=eye(8);
 B=[               0                  0                  0                  0;
                   0                  0                  0                  0;
    0.12433505058289   0.08278584480286  -2.75247764587402  -0.01788876950741;
   -0.03635892271996   0.47509527206421   0.01429074257612                  0;
    0.30449151992798   0.01495801657438  -0.49651837348938  -0.20674192905426;
    0.28773546218872  -0.54450607299805  -0.01637935638428                  0;
   -0.01907348632812   0.01636743545532  -0.54453611373901   0.23484230041504;
   -4.82063293457031  -0.00038146972656                  0                 0];
 C1=[eye(8);zeros(4,8)];
 D11=zeros(12,8);
 D12=[zeros(8,4);eye(4)];
 C=[0 0 0 0 0 0.0595 0.05329 -0.9968; 1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0; 0 0 0 -0.05348 1 0 0 0; 0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0];
 D21=zeros(6,8);
end

%------------------------------------------------------------------
% (HE5): A variation of the system above with eight state, two
%        measurement and four control variables. The matrices A and B
%        are the same as in (HE4).
%------------------------------------------------------------------
if(strmatch('HE5',ex,'exact'))
 nx=8; nu=4; ny=2; nz=4; nw=3;
 a01 = [          0                  0                  0   0.99857378005981;
                  0                  0   1.00000000000000  -0.00318221934140;
                  0                  0 -11.57049560546880  -2.54463768005371;
                  0                  0   0.43935656547546  -1.99818229675293;
                  0                  0  -2.04089546203613  -0.45899915695190;
 -32.10360717773440                  0  -0.50335502624512   2.29785919189453;
   0.10216116905212  32.05783081054690  -2.34721755981445  -0.50361156463623;
  -1.91097259521484   1.71382904052734  -0.00400543212891  -0.05741119384766];

 a02 = [0.05338427424431             0                  0                  0;
   0.05952465534210                  0                  0                  0;
  -0.06360262632370   0.10678052902222  -0.09491866827011   0.00710757449269;
                  0   0.01665188372135   0.01846204698086  -0.00118747074157;
  -0.73502779006958   0.01925575733185  -0.00459562242031   0.00212036073208;
                  0  -0.02121581137180  -0.02116791903973   0.01581159234047;
   0.83494758605957   0.02122657001019  -0.03787973523140   0.00035400385968;
                  0   0.01398963481188  -0.00090675335377  -0.29051351547241];
 
 A=[a01 a02];
 B1=[0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 0 0 0; 0 0 0];
 B=[               0                  0                  0                  0;
                   0                  0                  0                  0;
    0.12433505058289   0.08278584480286  -2.75247764587402  -0.01788876950741;
   -0.03635892271996   0.47509527206421   0.01429074257612                  0;
    0.30449151992798   0.01495801657438  -0.49651837348938  -0.20674192905426;
    0.28773546218872  -0.54450607299805  -0.01637935638428                  0;
   -0.01907348632812   0.01636743545532  -0.54453611373901   0.23484230041504;
   -4.82063293457031  -0.00038146972656                  0                 0];
 C1=[ 0        0         0         0         0    0.0595   0.05329  -0.9968;
      1.0      0         0         0         0         0         0        0;
      0      1.0         0         0         0         0         0        0;
      0        0         0  -0.05348       1.0         0         0        0];
 D11=zeros(4,3);
 D12=eye(nu);
 C =[   0        0       1.0         0         0         0         0       0;
        0        0         0       1.0         0         0         0       0];  
 D21=[0.01 0 0; 0 0.01 0];
end
%------------------------------------------------------------------
% (HE6): Helicopter control
%        "Multivariable feedback control: Analysis and design"
%        S. Skogestad and I. Postlethwaite
%        John Wiley and Sons, 1996, Section 12.2.3
%        Note: Matlab files 
%              http://www.nt.ntnu.no/users/skoge/book/matlab.html
%        stored in  /export/home/leibfr/Lipinski/matlab/..
%                                ..Examples_Multi_Feedback_Control/matlab_m/
%        F. Leibfritz, 29.10.2003
%        Data matrices H_inf mixed-sensitivity design
%        generated by Sec12_2.m (in directory above) on Laptop
%        cf. page 474, 475
%        save Heli_Sec12_2_3_Hinf A B1 B2 C1 C2 D11 D12 D21 D22
%------------------------------------------------------------------
if(strmatch('HE6',ex,'exact'))
 nx=20; nu=4; ny=6; nw=6; nz=16; 
 load he6  % ehemals: "Heli_Sec12_2_3_Hinf.mat"
 % .mat file containing A, B1, B2, C1, C2, D11, D12, D21, D22
 B=B2; C=C2;
end

%------------------------------------------------------------------
% (HE7): Helicopter control
%        "Multivariable feedback control: Analysis and design"
%        S. Skogestad and I. Postlethwaite
%        John Wiley and Sons, 1996, Section 12.2.4
%        difference: nw=9 --> B1, D11, D21,  rest like in (HE6)
%------------------------------------------------------------------
if(strmatch('HE7',ex,'exact'))
 nx=20; nu=4; ny=6; nw=9; nz=16; 
 load he7  % ehemals: "Heli_Sec12_2_4_dist_rejec.mat"
 % .mat file containing A, B1, B2, C1, C2, D11, D12, D21, D22
 B=B2; C=C2;
end

%------------------------------------------------------------------
% (JE1): Multivariable servomechanism problem (J-100 jet engine);
%        E. J. Davison and W. Gesing, "The systematic design of control
%        systems for for the multivariable servomechanism problem",
%        Nat. Eng. Consortium Inc., Chicago,
%        "Alternatives for Linear Multivariable Control", 1978 
%------------------------------------------------------------------
if(strmatch('JE1',ex,'exact')) 
 nx=30; nu=3; ny=5;
 load je1;  % .mat file containing A, B, C
 C1=sqrt(0.5)*[C; zeros(nu,nx)]; D12=sqrt(0.5)*[zeros(ny,nu); eye(nu)];
 B1=eye(nx);
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (JE2): Aero engine control
%        "Multivariable feedback control: Analysis and design"
%        S. Skogestad and I. Postlethwaite
%        John Wiley and Sons, 1996, Section 12.3.3
%        Note: Matlab files 
%              http://www.nt.ntnu.no/users/skoge/book/matlab.html
%        stored in  /export/home/leibfr/Lipinski/matlab/..
%                         ..Examples_Multi_Feedback_Control/matlab_m/
%        F. Leibfritz, 16.09.2003
%        Data matrices generated by Sec12_33.m (in directory above) on Laptop
%        save Aero_Engine a b c d A_Hinf B1 B2 C1 C2 D11 D12 D21 D22
%        a b c d -- data set G5 ==> [a,b,c,d]=unpck(W2GWpbar);
%        F. Leibfritz, 07.03.2005
%        renamed to A_je2=a; B_je2=b; C_je2=c; D_je2=d 
%------------------------------------------------------------------
if(strmatch('JE2',ex,'exact'))
 nx=21; nu=3; ny=3;
 load je2_3 
  % .mat file containing A_je2, B_je2, C_je2 --> A, B, C
 A=A_je2; B=B_je2; C=C_je2; B1=[]; C1=[]; D11=[]; D12=[]; D21=[];
end

%------------------------------------------------------------------
% (JE3): Aero engine control
%        "Multivariable feedback control: Analysis and design"
%        S. Skogestad and I. Postlethwaite
%        John Wiley and Sons, 1996, Section 12.3.3
%        Note: Matlab files 
%              http://www.nt.ntnu.no/users/skoge/book/matlab.html
%        stored in  /export/home/leibfr/Lipinski/matlab/..
%                         ..Examples_Multi_Feedback_Control/matlab_m/
%        F. Leibfritz, 16.09.2003
%        Data matrices generated by Sec12_33.m (in directory above) on Laptop
%        save Aero_Engine a b c d A_Hinf B1 B2 C1 C2 D11 D12 D21 D22
%        A_Hinf B1 B2 C1 C2 D11 D12 D21 D22 -- data set G5 for H_inf design
%              ==> A_Hinf := A in Sec12_33.m
%------------------------------------------------------------------
if(strmatch('JE3',ex,'exact'))
 nx=24; nu=3; ny=6; nw=6; nz=9;
 load je2_3 %%ehemals: "Aero_Engine.mat"
 % .mat file containing A_Hinf, B2, C2 (--> A, B, C) and C1, B1, D11, D12, D21
 A=A_Hinf; B=B2; C=C2; 
end

%------------------------------------------------------------------
% (REA1): The Chemical Reactor Example    %%ehemals (CHR2)
%         Y. S. Hung and A. G. J. MacFarlane, "Multivariable feedback:
%         A  quasi-classical approach", Springer-Verlag,
%         "Lecture Notes in Control and Information Sciences", 1982
%------------------------------------------------------------------
if(strmatch('REA1',ex,'exact')) 
 nx=4;nu=2;ny=2;
 A=[1.38 -0.2077 6.715 -5.676;-0.5814 -4.29 0 0.675;1.067 4.273 -6.654 5.893;
    0.048 4.273 1.343 -2.104];
 B=[0 0;5.679 0;1.136 -3.146;1.136 0];
 C=[1 0 1 -1;0 1 0 0;0 0 1 -1];
end

%------------------------------------------------------------------
% (REA2): Obtained from (REA1) by leaving out the last row of the
%         matrix C                                   %%ehemals (CHR1)
%------------------------------------------------------------------
if(strmatch('REA2',ex,'exact')) 
 nx=4;nu=2;ny=3;
 A=[1.4 -0.208 6.715 -5.676;-0.581 -4.29 0 .675;1.067 4.273 -6.654 5.893;
    .048 4.273 1.343 -2.104];
 B=[0 0;5.679 0;1.136 -3.146;1.136 0];
 C=[1 0 1 -1;0 1 0 0];
end

%------------------------------------------------------------------
% (REA3): Nuclear reactor model, L. F. Miller, R. G. Cochran, J. W. Howze
%         "Computation of Optimal Output Feedback Gains for Linear
%         Multivariable Systems", TOAC, Vol. 19, pp. 257--258, 1974
%------------------------------------------------------------------
if(strmatch('REA3',ex,'exact')) 
 nx=4;nu=2;ny=3;
 A=[-0.4044 , 0     , 0     , 0.4044, 0     , 0     , 0     , 0  , 0     , 0    , 0    , 0    ;
     0      ,-0.4044, 0     , 0     , 0.4044, 0     , 0     , 0  , 0     , 0    , 0    , 0    ;
     0      , 0     ,-0.4044, 0     , 0     , 0.4044, 0     , 0  , 0     , 0    , 0    , 0    ;
     0.01818, 0     , 0     ,-0.5363, 0     , 0     , 0.4045, 0  , 0     , 0    , 0    , 0    ;
     0      , 0.0818, 0     , 0.4545,-0.5363, 0     , 0     , 0  , 0     , 0    , 0    , 0    ;
     0      , 0     , 0.0818, 0     , 0.4545,-0.5363, 0     , 0  , 0     , 0    , 0    , 0    ;
     0      , 0     , 0     , 0     , 0.15  , 0     ,-0.15  , 0  , 0     , 0    , 0    , 0    ;
     0      , 0     , 0     , 0     , 0     , 0     , 0     , 0  , 0     , 0    , 0    , 0    ;
     0      ,-7.5   , 0     , 0     , 75    , 0     , 0     , 600,-74.995, 0.033, 0.346, 0.621;
     0      , 0     , 0     , 0     , 0     , 0     , 0     , 0  , 2.475 ,-0.033, 0    , 0    ;
     0      , 0     , 0     , 0     , 0     , 0     , 0     , 0  , 25.95 , 0    ,-0.346, 0    ;
     0      , 0     , 0     , 0     , 0     , 0     , 0     , 0  , 46.57 , 0    , 0    ,-0.621];
 B=[0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
 C=[0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ];
end

%------------------------------------------------------------------
% (REA4): Chemical reactor model by P. M. Maekilae, "Parametric LQ
%         Control", IJOC, Vol. 41, Nr. 6, pp. 1413-1428, 1985 
%         (discrete modell)
%------------------------------------------------------------------
if(strmatch('REA4',ex,'exact'))
 nx=8;nu=1;ny=1;
 A=[0.5623 ,-0.01642, 0.01287,-0.0161 , 0.02094,-0.02988, 0.0183 , 0.008743;
    0.102  , 0.6114 ,-0.02468, 0.02468,-0.03005, 0.04195,-0.02559, 0.03889 ;
    0.1361 , 0.2523 , 0.641  ,-0.03404, 0.03292,-0.04296, 0.02588, 0.08467 ;
    0.09951, 0.2859 , 0.3476 , 0.6457 ,-0.03249, 0.03316,-0.01913, 0.1103  ;
   -0.04794, 0.08708, 0.3297 , 0.3102 , 0.6201 ,-0.03015, 0.01547, 0.08457 ;
   -0.1373 ,-0.1224 , 0.1705 , 0.3106 , 0.191  , 0.5815 ,-0.01274, 0.05394 ;
   -0.1497 ,-0.1692 , 0.1165 , 0.2962 , 0.1979 , 0.07631, 0.5242 , 0.04702 ;
    0      , 0      , 0      , 0      , 0      , 0      , 0      , 0.6065  ];
 B=[-0.1774;-0.2156;-0.2194;-0.09543; 0.0579; 0.09303; 0.08962; 0];
 C=[-0.0049 , 0.0049 ,-0.006 , 0.01 , 0.0263 , 0.3416 , 0.6759 , 0];
 B1=[0; 0; 0; 0; 0; 0; 0; 1];
 C1=[-0.0465 ,-0.1135 ,-0.1909 ,-0.2619 ,-0.2634 ,-0.1422 ,-0.0002 , 0.1856];
 D12=[0.1001];
 D11=[0];
 D21=[1];
end

%------------------------------------------------------------------
% (DIS1): Decentralized Interconnected System;
%         H. Singh, R. H. Brown and D. S. Naidu, "Unified approach to
%         linear quadrtic regulator with time-scale property",
%         Optimal Control Applications and Methods, Vol.22, pp.1-16, 2001      
%------------------------------------------------------------------
if(strmatch('DIS1',ex,'exact'))
 A11=[0.144 -0.058 0.056 0.042; -0.506 -0.236 -0.020 -0.012;
     0 0 -0.278 0.291;0 0 0 -0.330];
 A12=[0.120 2.1454 0 0.080;-0.06 -0.909 1.093 -0.040; 0 0 0 0.580; 0 0 0 0]; 
 A21=[0 0 0.303 0.029;-0.154 0.133 -0.006 -0.004;
     -0.345 0.304 -0.018 -0.014; 0 0 0 0.247];  
 A22=[-1.67 0 0 0.092;-0.014 -1.688 0.236 0.013;-0.032 -0.611 -1.824 -0.024;
     0 0 0 -1.978];
 A=[A11 A12;A21 A22];
 Bb1=[-0.076 0.02;0.588 -0.006;0 0.152; 0 1.45];
 Bb2=[0 0.012;0.162 -0.002;0.414 -0.008;0 0.248]; 
 B=[Bb1, zeros(4,2); zeros(4,2) Bb2];
 Ch=eye(8);
 C=Ch(1:4,:);
 C1=[C; zeros(4,8)];
 B1=[1;0;1;0;1;0;1;0];
 D11=zeros(8,1);
 D12=[zeros(4,4);eye(4)];
 D21=zeros(4,1);  
end

%------------------------------------------------------------------
% (DIS2): Decentralized system with 2 control stations
%         W. Q. Liu and V. Sreeram, "New Algorithm for Computing LQ
%         Suboptimal Output Feedback Gains of Decentralized Control
%         Systems", JOTA, Vol. 93   
%------------------------------------------------------------------
if(strmatch('DIS2',ex,'exact'))
 A=[-4 2 1;3 -2 5;-7 0 3];
 Bb1=[1;1;0]; Bb2=[0;0;1];
 B=[Bb1 Bb2];
 Cc1=[0 1 0]; Cc2=[0 0 1];
 C=[Cc1;Cc2];
end

%------------------------------------------------------------------
% (DIS3): M. Saif and Y. Guan,"Decentralized State Estimation in
%         Large-Scale Interconnected Dynamical Systems",
%         AUTO, Vol. 28, Nr. 1, pp. 215-219 
%------------------------------------------------------------------
if(strmatch('DIS3',ex,'exact'))
 A=[-1.0 0.0 0.0 0.0 0.0 0.0;-1.0 1.0 1.0 0.0 0.0 0.0;
    1.0 -2.0 -1.0 -1.0 1.0 1.0;0.0 0.0 0.0 -1.0 0.0 0.0;
   -8.0 1.0 -1.0 -1.0 -2.0 0.0;4.0 -0.5 0.5 0.0 0.0 -4.0];
 B=[0 1 0 0; 1 0 0 0; 1 1 0 0; 0 0 0 -1; 0 0 1 0; 0 0 0 1]; 
 C=[0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];		   
end

%------------------------------------------------------------------
% (DIS4): H. T. Toivonen and P. M. Maekilae, "A descent Anderson-
%         Moore algorithm for optimal decentralized control",
%         AUTO, Vol. 21, Nr. 6, pp.743-744, 1985 
%------------------------------------------------------------------
if(strmatch('DIS4',ex,'exact'))
 A=[0 1 0.5 1   0.6  0  ;-2 -3 1 0 0 1; 0  2 0.5 1   1   0.5;
    1 3 0   0.5 0   -0.5; 0  1 1 0 1 0;-3 -4 0   0.5 0.5 0  ];
 B=[1 0 0 0;1 0 0 0;0 3 0 0;0 0 4 0;0 0 0 2;0 0 0 3]; 
 C=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 1];	   
end

%------------------------------------------------------------------
% (DIS5): M. C. de Oliveira, J. F. Camino and R. E. Skelton,
%         A Convexifying Algorithm of Structured Linear Controllers
%         Tech. Report, FAPESP and CAPES, Brazil
%         (discrete model)
%------------------------------------------------------------------
if(strmatch('DIS5',ex,'exact'))
 A=[0.8189  0.0863 0.09    0.0813;
    0.2524  1.0033 0.0313  0.2004;
   -0.0545  0.0102 0.7901 -0.258 ;
   -0.1918 -0.1034 0.1602  0.8604];
 B=[0.0045  0.0044;
    0.1001  0.01  ;
    0.0003 -0.0136;
   -0.0051  0.0936];
 C=[1 0 0 0; 0 0 1 0];
 B1=[0.0953 0 0;
     0.0145 0 0;
     0.0862 0 0;
    -0.0011 0 0];
 C1=[1 0 -1 0; 0 0 0 0;0 0 0 0];  
 D11=zeros(3,3);
 D12=[0 0; 1 0; 0 1];
 D21=[0 1 0; 0 0 1];	   
end 

%------------------------------------------------------------------
% (TG1): Turbo-Generator
%        Y. S. Hung and A. G. J. MacFarlane, "Multivariable feedback:
%        A  quasi-classical approach", Springer-Verlag,
%        "Lecture Notes in Control and Information Sciences", 1982
%        p. 117/167
%------------------------------------------------------------------
if(strmatch('TG1',ex,'exact')) 
 nx=10;nu=2;ny=2;
 A=[0 1 0 0 0 0 0 0 0 0;
    0 -0.11323 -.98109 -11.847 -11.847 -63.08 -34.339 -34.339 -27.645 0;
    324.121 -1.1755 -29.101 .12722 2.83448 -967.73 -678.14 -678.14 0 -129.29;
    -127.3 0.46167 11.4294 -1.0379 13.1237 380.079 266.341 266.341 0 1054.85;
    -186.05 .67475 16.7045 0.86092 -17.068 555.502 389.268 389.268 0 -874.92;
    341.917 1.09173 1052.75 756.465 756.465 -29.774 0.16507 3.27626 0 0;
   -30.748 -.09817 -94.674 -68.029 -68.029 2.67753 -2.6558 4.88497 0 0;
   -302.36 -.96543 -930.96 -668.95 -668.95 26.3292 2.42028 -9.5603 0 0;
   0 0 0 0 0 0 0 0 -1.6667 0;
   0 0 0 0 0 0 0 0 0 -10];
 B=[0 0;0 0;0 0;0 0;0 0;0 0;0 0;0 0; 1.6667 0;0 10];
 C=[1 0 0 0 0 0 0 0 0 0;-.49134 0 -.63203 0 0 -.20743 0 0 0 0];
end

%------------------------------------------------------------------
% (AGS): Automobile Gas Turbine
%        Y. S. Hung and A. G. J. MacFarlane, "Multivariable feedback:
%        A  quasi-classical approach", Springer-Verlag,
%        "Lecture Notes in Control and Information Sciences", 1982
%        p. 27/163
%------------------------------------------------------------------
if(strmatch('AGS',ex,'exact'))
 nx=12;nu=2;ny=2;
 A=[0 1 0 0 0 0 0 0 0 0 0 0;
   -0.202 -1.15 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 -2.36 -13.6 -12.8 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 -1.62 -9.4 -9.15 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 -188 -111.6 -116.4 -20.8];
 B=[0 0;1.0439 4.1486;0 0;0 0;-1.794 2.6775;0 0;0 0;1.0439 4.1486;0 0;0 0;
    0 0;-1.794 2.6775];
 C=[0.264 0.806 -1.42 -15 0 0 0 0 0 0 0 0;
    0 0 0 0 0 4.9 2.12 1.95 9.35 25.8 7.14 0];
end

%------------------------------------------------------------------
% (WEC1): Wind energy conversion system
%         "Dynamic modelling and robust control of a wind energy
%         conversion system"
%         Maarten Steinbuch, 1989, PHD-Thesis University of Delft
%         Appendix A.5: Linear models 
%         operation point v=12 m/s (v : wind speed)
%------------------------------------------------------------------
if(strmatch('WEC1',ex,'exact'))
 A=zeros(10,10); 
 A(1,1)=-5.0e0; 
 A(2,3)=1.0e0;
 A(3,1:4)=[-5.5005e0 -1.4791e3 -3.2812e0 -1.7889e-2];
 A(3,7:10)=[1.6968e2 3.6137e1 3.6137e1 1.4483e2];
 A(4,2:3)=[1.4164e3 3.1250e0];
 A(4,7:10)=[-1.6968e2 -3.6137e1 -3.6137e1 -1.4483e2];
 A(5,4:5)=[9.5493e-2 -1.0e1];
 A(6,6)=-1.0e1;
 A(7,4)=7.8416e0;
 A(7,6:10)=[1.1552e-1 -1.2571e3 1.0151e3 1.0111e3 4.9909e2];
 A(8,4)=4.6042e0;
 A(8,6:10)=[ 2.0960e0 -6.9313e2 5.5933e2 6.3131e2 3.0618e2];
 A(9,4)=5.7968e0;
 A(9,6:10)=[-1.8671e0 -9.7681e2 7.8851e2 7.0825e2 3.5508e2];
 A(10,4)=-2.8663e0;
 A(10,6:10)=[-4.7856e-2 4.1358e2 -3.4335e2 -3.4163e2 -2.1245e2];
 
 B=zeros(10,3);
 B(1,1)=5.0e0; B(6,2)=1.0e1;
 B(7:10,3)=[-3.0565e2;-1.6627e2;-2.3988e2; 9.6020e1];
 
 C=zeros(4,10);
 C(1,5)=1; C(2,7)=1;
 C(3,3)=4.5455e-2; C(3,4)=C(3,3);
 C(4,2)=1.2249e1; C(4,3)=2.7025e-2;
end

%------------------------------------------------------------------
% (WEC2): like (WEC1) at an operation point of v=16m/s (v: wind speed)
%------------------------------------------------------------------
if(strmatch('WEC2',ex,'exact'))
 A=zeros(10,10); 
 A(1,1)=-5.0e0; 
 A(2,3)=1.0e0;
 A(3,1:4)=[-3.3690e1 -1.4791e3 -3.3531e0 -8.9802e-2];
 A(3,7:10)=[1.6968e2 3.6137e1 3.6137e1 1.4483e2];
 A(4,2:3)=[1.4164e3 3.1250e0];
 A(4,7:10)=[-1.6968e2 -3.6137e1 -3.6137e1 -1.4483e2];
 A(5,4:5)=[9.5493e-2 -1.0e1];
 A(6,6)=-1.0e1;
 A(7,4)=7.8416e0;
 A(7,6:10)=[1.1552e-1 -1.2571e3 1.0151e3 1.0111e3 4.9909e2];
 A(8,4)=4.6042e0;
 A(8,6:10)=[ 2.0960e0 -6.9313e2 5.5933e2 6.3131e2 3.0618e2];
 A(9,4)=5.7968e0;
 A(9,6:10)=[-1.8671e0 -9.7681e2 7.8851e2 7.0825e2 3.5508e2];
 A(10,4)=-2.8663e0;
 A(10,6:10)=[-4.7856e-2 4.1358e2 -3.4335e2 -3.4163e2 -2.1245e2];
 
 B=zeros(10,3);
 B(1,1)=5.0e0; B(6,2)=1.0e1;
 B(7:10,3)=[-3.0565e2;-1.6627e2;-2.3988e2; 9.6020e1];
 
 C=zeros(4,10);
 C(1,5)=1; C(2,7)=1;
 C(3,3)=4.5455e-2; C(3,4)=C(3,3);
 C(4,2)=1.2249e1; C(4,3)=2.7025e-2;
end

%------------------------------------------------------------------
% (WEC3): like (WEC1) at an operation point of v=20m/s (v: wind speed)
%------------------------------------------------------------------
if(strmatch('WEC3',ex,'exact'))
 A=zeros(10,10); 
 A(1,1)=-5.0e0; 
 A(2,3)=1.0e0;
 A(3,1:4)=[-7.0878e1 -1.4791e3 -3.4321e0 -1.6877e-1];
 A(3,7:10)=[1.6968e2 3.6137e1 3.6137e1 1.4483e2];
 A(4,2:3)=[1.4164e3 3.1250e0];
 A(4,7:10)=[-1.6968e2 -3.6137e1 -3.6137e1 -1.4483e2];
 A(5,4:5)=[9.5493e-2 -1.0e1];
 A(6,6)=-1.0e1;
 A(7,4)=7.8416e0;
 A(7,6:10)=[1.1552e-1 -1.2571e3 1.0151e3 1.0111e3 4.9909e2];
 A(8,4)=4.6042e0;
 A(8,6:10)=[ 2.0960e0 -6.9313e2 5.5933e2 6.3131e2 3.0618e2];
 A(9,4)=5.7968e0;
 A(9,6:10)=[-1.8671e0 -9.7681e2 7.8851e2 7.0825e2 3.5508e2];
 A(10,4)=-2.8663e0;
 A(10,6:10)=[-4.7856e-2 4.1358e2 -3.4335e2 -3.4163e2 -2.1245e2];
 
 B=zeros(10,3);
 B(1,1)=5.0e0; B(6,2)=1.0e1;
 B(7:10,3)=[-3.0565e2;-1.6627e2;-2.3988e2; 9.6020e1];
 
 C=zeros(4,10);
 C(1,5)=1; C(2,7)=1;
 C(3,3)=4.5455e-2; C(3,4)=C(3,3);
 C(4,2)=1.2249e1; C(4,3)=2.7025e-2;
end

%------------------------------------------------------------------
% (HF1): Heat flow in a thin rod (1D model)
%        A. S. Hodel, K. P. Poolla and B. Tension, "Numerical Solution
%        of the Lyapunpv Equation by Approximate Power Iteration",
%        Linear Algebra Appl., Vol. 236, pp. 205-230, 1996
%------------------------------------------------------------------
if(strmatch('HF1',ex,'exact'))
 nx=130;nu=1;ny=2;
 h=1/(nx+1);
 ad=-2*ones(nx,1); ad(1)=-1.0; an=ones(nx,1);
 Ah=spdiags([an ad an],-1:1,nx,nx); 
 A=(1/h)*full(Ah);
 B=zeros(nx,1); B(nx,1)=1/h;
 C=zeros(2,nx); C(1,1)=1; C(2,nx/2)=1;
 B1=B; C1=C;
 D12=[0;1];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (BDT1): Binary distillation tower with pressure variation
%         E. J. Davison, "Benchmark Problems for Control System Design",
%         "Report of the IFAC Theory Comittee", 1990 
%------------------------------------------------------------------
if(strmatch('BDT1',ex,'exact'))
 nx=11;nu=3;ny=3;
 A=[-0.014 0.0043 0 0 0 0 0 0 0 0 0;
    0.0095 -0.0138 0.0046 0 0 0 0 0 0 0 0.0005;
    0 0.0095 -0.0141 0.0063 0 0 0 0 0 0 0.0002;
    0 0 0.0095 -0.0158 0.011 0 0 0 0 0 0;
    0 0 0 0.0095 -0.0312 0.015 0 0 0 0 0;
    0 0 0 0 0.0202 -0.0352 0.022 0 0 0 0;
    0 0 0 0 0 0.0202 -0.0422 0.028 0 0 0;
    0 0 0 0 0 0 0.0202 -0.0482 0.037 0 0.0002;
    0 0 0 0 0 0 0 0.0202 -0.0572 0.042 0.0005;
    0 0 0 0 0 0 0 0 0.0202 -0.0483 0.0005;
    0.0255 0 0 0 0 0 0 0 0 0.0255 -0.0185];
 B=[0 0 0; 0.000005 -0.00004 0.0025; 0.000002 -0.00002 0.005;
    0.000001 -0.00001 0.005; 0 0 0.005; 0 0 0.005; -0.000005 0.00001 0.005;
   -0.00001 0.00003 0.005; -0.00004 0.000005 0.0025; -0.00002 0.000002 0.0025;
    0.00046 0.00046 0];
 C=[0 0 0 0 0 0 0 0 0 1 0; 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1];
 B1=zeros(11,1); B1(5,1)=1.0e-2; C1=[C ; zeros(3,nx)];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D12=[zeros(nz-nu,nu); eye(nu)];
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (BDT2): Binary distillation tower
%         Test and generate the linearized 82-state distillation column
%         "Multivariable feedback control: Analysis and design"
%         S. Skogestad and I. Postlethwaite
%         John Wiley and Sons, 1996, Section 12.4
%         Note: Matlab files 
%            http://www.nt.ntnu.no/users/skoge/book/matlab.html
%         stored in  /export/home/leibfr/Lipinski/matlab/..
%                         ..Examples_Multi_Feedback_Control/matlab_m/cola
%         F. Leibfritz, 29.10.2003
%         save Dist_Column_82 A B B1 C C1 D12
%------------------------------------------------------------------
if(strmatch('BDT2',ex,'exact'))
 nx=82; nu=4; ny=4; nz=4; nw=2;
 load bdt2 %-- generated by test_cola_lin.m (in directory above)
           % ehemals: "Dist_Column_82.mat"
 A=full(A); B=B; C=C;
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (MFP): Moored Floating Platform  
%        C. Scherer, P. Gahinet and M. Chilali, "Multiobjective Output-
%        Feedback Control via LMI Optimization",
%        TOAC, Vol. 42, Nr. 7, pp. 896-911, 1997    
%------------------------------------------------------------------
if(strmatch('MFP',ex,'exact'))
 nx=4;nu=3;ny=2;
 A=[0 0 1 0;0 0 0 1;-0.101 -0.1681 -0.04564 -0.01075;
    0.06082 -2.1407 -0.05578 -0.1273];
 B=[0 0 0;0 0 0;0.1179 0.1441 0.1476;0.1441 1.7057 -0.7557];
 C=[1 0 0 0;0 1 0 0];
end

%------------------------------------------------------------------
% (UWV): Control surface servo for an underwater vehicle,
%         E. J. Davison, "Benchmark Problems for Control System Design",
%         "Report of the IFAC Theory Comittee", 1990; p.32
%------------------------------------------------------------------
if(strmatch('UWV',ex,'exact'))
 nx=8;nu=2;ny=2;nw=2;nz=1;
 A=[  0   , 850 , 0    , 0   , 0    , 0  , 0   , 0  ;
     -850 ,-120 ,-4100 , 0   , 0    , 0  , 0   , 0  ;
      33  , 0   ,-33   , 0   ,-700  , 0  , 0   , 0  ;
      0   , 0   , 0    , 0   , 1400 , 0  , 0   , 0  ;
      0   , 0   , 1600 ,-450 ,-110  , 0  , 0   , 0  ;
      0   , 0   , 0    , 81  , 0    ,-1  , 0   ,-900;
      0   , 0   , 0    , 0   , 0    , 0  , 0   , 110;
      0   , 0   , 0    , 0   , 0    , 12 ,-1.1 ,-22 ];
 B=[  0 , 0; 4.6 , 99000; 0 , 0; 0 , 0;
      0 , 0; 0   , 0    ; 0 , 0; 0 , 0];
 C=[  0 , 0 , 0 , 0 , 0 , 1 , 0 , 0;
      0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ];
 B1=[0 , 0; 9900 , 0; 0 , 0; 0 , 0 ;
     0 , 0; 0    , 0; 0 , 0; 0 , 99];
 C1=[0 , 0 , 0 , 0 , 0 , 0 , 1 , 0];
 D12=[1,0];
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (IH): Intelligent highway: model of position and velocity control
%       for a string of high-speed vehicles;                   %%ehemals (IH1)
%       J. Abels and P. Benner, "CAREX - A Collection of Benchmark Examples
%       for Continuous-Time Algebraic Riccati Equations (Version 2.0)",
%       SLICOT Working Note 1999-14, Ex. 3.1
%       available via ftp: wgs.esat.kuleuven.ac.be/
%                          pub/WGS/REPORTS/SLWN1999-14.ps.Z",
%------------------------------------------------------------------
if(strmatch('IH',ex,'exact'))
  nx=21; nu=11; ny=10;
  load ih;
  C1=sqrt(5)*[C;zeros(1,nx)]; D12=sqrt(0.5)*eye(nu);
  B1=eye(nx);
  [nx,nw]=size(B1); [nz,nx]=size(C1);
  D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (CSE1): Coupled spring experiment,  l=10  (2nd order system)
%         J. Abels and P. Benner, "CAREX - A Collection of Benchmark Examples
%         for Continuous-Time Algebraic Riccati Equations (Version 2.0)",
%         SLICOT Working Note 1999-14, Ex. 4.3
%         available via ftp: wgs.esat.kuleuven.ac.be/
%                          pub/WGS/REPORTS/SLWN1999-14.ps.Z",
%------------------------------------------------------------------
if(strmatch('CSE1',ex,'exact'))
 l=10; 
 mu=4.0e0; delta=mu; kappa=1.0e0;
 nx=2*l; nu=2; ny=l;
 M=mu*eye(l); L=delta*eye(l); N=eye(l); P=eye(l);
 e=ones(l,1); 
 K=full(spdiags([-e 2*e -e], -1:1,l,l)); K(1,1)=1; K(l,l)=1;
 K=kappa*K;
 D=zeros(l,nu); D(1,1)=1; D(l,2)=-1;
  
 A = [zeros(l,l) eye(l); -inv(M)*K -inv(M)*L];
 B = [zeros(l,nu); inv(M)*D];
 B1= -0.02*B(:,2);
 C =[N P];
 C1=[C;zeros(nu,nx)];
 D12=[zeros(ny,nu);eye(nu)];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (CSE2): Coupled spring experiment,   l=30   (2nd order system)
%         J. Abels and P. Benner, "CAREX - A Collection of Benchmark Examples
%         for Continuous-Time Algebraic Riccati Equations (Version 2.0)",
%         SLICOT Working Note 1999-14, Ex. 4.3
%         available via ftp: wgs.esat.kuleuven.ac.be/
%                          pub/WGS/REPORTS/SLWN1999-14.ps.Z",
%------------------------------------------------------------------
if(strmatch('CSE2',ex,'exact'))
 l=30; 
 mu=4.0e0; delta=mu; kappa=1.0e0;
 nx=2*l; nu=2; ny=l;
 M=mu*eye(l); L=delta*eye(l); N=eye(l); P=eye(l);
 e=ones(l,1); 
 K=full(spdiags([-e 2*e -e], -1:1,l,l)); K(1,1)=1; K(l,l)=1;
 K=kappa*K;
 D=zeros(l,nu); D(1,1)=1; D(l,2)=-1;
  
 A = [zeros(l,l) eye(l); -inv(M)*K -inv(M)*L];
 B = [zeros(l,nu); inv(M)*D];
 B1= -0.02*B(:,2);
 C =[N P];
 C1=[C;zeros(nu,nx)];
 D12=[zeros(ny,nu);eye(nu)];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (EB1): Euler-Bernoulli Beam; J. C. Geromel and P. B. Gapski
%        "Synthesis of positive real H2 controllers"
%	 TOAC, Vol. 42, Nr. 7, pp. 988-992, 1997
%        low damping --> xi=1e-2
%        s_max=5 --> nx=10 
%------------------------------------------------------------------
if(strmatch('EB1',ex,'exact'))
 nx=10;nu=1;ny=1;
 O=zeros(2,2); xi=0.01;
 s=1; w=s^2; A1=[0, 1; -w^2 -2*xi*w];
 s=2; w=s^2; A2=[0, 1; -w^2 -2*xi*w];
 s=3; w=s^2; A3=[0, 1; -w^2 -2*xi*w];
 s=4; w=s^2; A4=[0, 1; -w^2 -2*xi*w];
 s=5; w=s^2; A5=[0, 1; -w^2 -2*xi*w];
 A=[A1 O O O O; O A2 O O O; O O A3 O O; O O O A4 O; O O O O A5];
 B=[0;0.9877;0;-0.309;0;-0.891;0;0.5878;0;0.7071];
 C=B';
 B1=[0 0;0.9877 0;0 0;-0.309 0;0 0;-0.891 0;0 0;0.5878 0;0 0;0.7071 0];
 C1=[0 0.809 0 -0.9511 0 0.309 0 0.5878 0 -1;
     0 0     0  0      0 0     0 0      0  0];
 D12=[0;1.9];
 D21=[0 1.9];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (EB2): Euler-Bernoulli Beam; like (EB1) with changed performance
%        criteria, (-> C1, D12)
%        low damping --> xi=1e-2
%        s_max=5 --> nx=10  
%------------------------------------------------------------------
if(strmatch('EB2',ex,'exact'))
 nx=10;nu=1;ny=1;
 O=zeros(2,2); xi=0.01;
 s=1; w=s^2; A1=[0, 1; -w^2 -2*xi*w];
 s=2; w=s^2; A2=[0, 1; -w^2 -2*xi*w];
 s=3; w=s^2; A3=[0, 1; -w^2 -2*xi*w];
 s=4; w=s^2; A4=[0, 1; -w^2 -2*xi*w];
 s=5; w=s^2; A5=[0, 1; -w^2 -2*xi*w];
 A=[A1 O O O O; O A2 O O O; O O A3 O O; O O O A4 O; O O O O A5];
 B=[0;0.9877;0;-0.309;0;-0.891;0;0.5878;0;0.7071];
 C=B';
 B1=[0 0;0.9877 0;0 0;-0.309 0;0 0;-0.891 0;0 0;0.5878 0;0 0;0.7071 0];
 C1=[0.809 0 -0.9511 0 0.309 0 0.5878 0 -1 0;
     0     0  0      0 0     0 0      0  0 0];
 D12=[0;0.5];
 D21=[0 1.9];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
end

%------------------------------------------------------------------
% (EB3): Euler-Bernoulli Beam; like (EB2)
%        very low damping --> xi=1e-7  
%        s_max=5 --> nx=10
%------------------------------------------------------------------
if(strmatch('EB3',ex,'exact'))
 xi=1e-7; %1e-2; 1e-5
 s_max=5; % 5; 10; 20; 80
 s=(1:s_max)'; %-- no. of modal coordinates q_r, r=1,..,s_max
 w=s.^2;       %-- w_r:=r^2
 p_a=1.7279e0; %-- location p_a of the point force actuator
 p_z=2.1991e0; %-- location for controlled output z
 i=1; Ah=[0, 1; -w(i)^2 -2*xi*w(i)]; A=Ah;
      Bh=[0; sin(s(i)*p_a)]; B=Bh;
      %C1h=[0, sin(s(i)*p_z)]; C1=C1h; %-- performance defined in terms of
                                      %-- velocity at p_z
      C1h=[sin(s(i)*p_z), 0]; C1=C1h; %-- performance defined in terms of
                                      %-- transverse deflection at p_z		      
 for i=2:s_max
   Ah=[0, 1; -w(i)^2 -2*xi*w(i)];
   A=blkdiag(A,Ah);
   Bh=[0; sin(s(i)*p_a)];
   B=[B;Bh];
   %C1h=[0, sin(s(i)*p_z)]; %-- performance defined in terms of
   %C1=[C1, C1h];           %-- velocity at p_z
   C1h=[sin(s(i)*p_z), 0]; %-- performance defined in terms of
   C1=[C1, C1h];           %-- transverse deflection at p_z
 end
 C=B';
 B1=[B, zeros(2*s_max,1)];
 D21=[0, 1.9];
 %-- performance defined in terms of velocity at p_z
 %C1=[C1; zeros(1,2*s_max)]; D12=[0;1.9];
 %-- performance defined in terms of transverse deflection at p_z
 C1=[C1; zeros(1,2*s_max)]; D12=[0;0.5];
 D11=zeros(2,2);
end

%------------------------------------------------------------------
% (EB4): Euler-Bernoulli Beam; like (EB2)
%        very low damping --> xi=1e-7
%        higher dimension:  s_max=10 --> nx=20
%------------------------------------------------------------------
if(strmatch('EB4',ex,'exact'))
 xi=1e-7; %1e-2; 1e-5
 s_max=10; % 5; 10; 20; 80 
 s=(1:s_max)'; %-- no. of modal coordinates q_r, r=1,..,s_max
 w=s.^2;       %-- w_r:=r^2
 p_a=1.7279e0; %-- location p_a of the point force actuator
 p_z=2.1991e0; %-- location for controlled output z
 i=1; Ah=[0, 1; -w(i)^2 -2*xi*w(i)]; A=Ah;
      Bh=[0; sin(s(i)*p_a)]; B=Bh;
      %C1h=[0, sin(s(i)*p_z)]; C1=C1h; %-- performance defined in terms of
                                      %-- velocity at p_z
      C1h=[sin(s(i)*p_z), 0]; C1=C1h; %-- performance defined in terms of
                                      %-- transverse deflection at p_z		      
 for i=2:s_max
   Ah=[0, 1; -w(i)^2 -2*xi*w(i)];
   A=blkdiag(A,Ah);
   Bh=[0; sin(s(i)*p_a)];
   B=[B;Bh];
   %C1h=[0, sin(s(i)*p_z)]; %-- performance defined in terms of
   %C1=[C1, C1h];           %-- velocity at p_z
   C1h=[sin(s(i)*p_z), 0]; %-- performance defined in terms of
   C1=[C1, C1h];           %-- transverse deflection at p_z
 end
 C=B';
 B1=[B, zeros(2*s_max,1)];
 D21=[0, 1.9];
 %-- performance defined in terms of velocity at p_z
 %C1=[C1; zeros(1,2*s_max)]; D12=[0;1.9];
 %-- performance defined in terms of transverse deflection at p_z
 C1=[C1; zeros(1,2*s_max)]; D12=[0;0.5];
 D11=zeros(2,2); 
end

%------------------------------------------------------------------
% (EB5): Euler-Bernoulli Beam; like (EB2)
%        very low damping --> xi=1e-7  
%        higher dimension:  s_max=20 --> nx=40
%------------------------------------------------------------------
if(strmatch('EB5',ex,'exact'))
 xi=1e-7; %1e-2; 1e-5
 s_max=20;%5; 10; 20; 80
 s=(1:s_max)'; %-- no. of modal coordinates q_r, r=1,..,s_max
 w=s.^2;       %-- w_r:=r^2
 p_a=1.7279e0; %-- location p_a of the point force actuator
 p_z=2.1991e0; %-- location for controlled output z
 i=1; Ah=[0, 1; -w(i)^2 -2*xi*w(i)]; A=Ah;
      Bh=[0; sin(s(i)*p_a)]; B=Bh;
      %C1h=[0, sin(s(i)*p_z)]; C1=C1h; %-- performance defined in terms of
                                      %-- velocity at p_z
      C1h=[sin(s(i)*p_z), 0]; C1=C1h; %-- performance defined in terms of
                                      %-- transverse deflection at p_z		      
 for i=2:s_max
   Ah=[0, 1; -w(i)^2 -2*xi*w(i)];
   A=blkdiag(A,Ah);
   Bh=[0; sin(s(i)*p_a)];
   B=[B;Bh];
   %C1h=[0, sin(s(i)*p_z)]; %-- performance defined in terms of
   %C1=[C1, C1h];           %-- velocity at p_z
   C1h=[sin(s(i)*p_z), 0]; %-- performance defined in terms of
   C1=[C1, C1h];           %-- transverse deflection at p_z
 end
 C=B';
 B1=[B, zeros(2*s_max,1)];
 D21=[0, 1.9];
 %-- performance defined in terms of velocity at p_z
 %C1=[C1; zeros(1,2*s_max)]; D12=[0;1.9];
 %-- performance defined in terms of transverse deflection at p_z
 C1=[C1; zeros(1,2*s_max)]; D12=[0;0.5];
 D11=zeros(2,2); 
end

%------------------------------------------------------------------
% (EB6): Euler-Bernoulli Beam; like (EB2)
%        very low damping --> xi=1e-7  
%        higher dimension:  s_max=80 --> nx=160
%------------------------------------------------------------------
if(strmatch('EB6',ex,'exact'))
 xi=1e-7; %1e-2; 1e-5
 s_max=80;%5; 10; 20; 80 
 s=(1:s_max)'; %-- no. of modal coordinates q_r, r=1,..,s_max
 w=s.^2;       %-- w_r:=r^2
 p_a=1.7279e0; %-- location p_a of the point force actuator
 p_z=2.1991e0; %-- location for controlled output z
 i=1; Ah=[0, 1; -w(i)^2 -2*xi*w(i)]; A=Ah;
      Bh=[0; sin(s(i)*p_a)]; B=Bh;
      %C1h=[0, sin(s(i)*p_z)]; C1=C1h; %-- performance defined in terms of
                                      %-- velocity at p_z
      C1h=[sin(s(i)*p_z), 0]; C1=C1h; %-- performance defined in terms of
                                      %-- transverse deflection at p_z		      
 for i=2:s_max
   Ah=[0, 1; -w(i)^2 -2*xi*w(i)];
   A=blkdiag(A,Ah);
   Bh=[0; sin(s(i)*p_a)];
   B=[B;Bh];
   %C1h=[0, sin(s(i)*p_z)]; %-- performance defined in terms of
   %C1=[C1, C1h];           %-- velocity at p_z
   C1h=[sin(s(i)*p_z), 0]; %-- performance defined in terms of
   C1=[C1, C1h];           %-- transverse deflection at p_z
 end
 C=B';
 B1=[B, zeros(2*s_max,1)];
 D21=[0, 1.9];
 %-- performance defined in terms of velocity at p_z
 %C1=[C1; zeros(1,2*s_max)]; D12=[0;1.9];
 %-- performance defined in terms of transverse deflection at p_z
 C1=[C1; zeros(1,2*s_max)]; D12=[0;0.5];
 D11=zeros(2,2); 
end

%------------------------------------------------------------------
% (PAS): Piezoelectric actuator system
%        B. M. Chen, "H_inf Control and Its Applications",
%        Springer-Verlag, "Lecture Notes in Control and Information Sciences",
%        Vol.235, 1998; p.283
%------------------------------------------------------------------
if(strmatch('PAS',ex,'exact'))
 nx=5;nu=1;ny=3;
 A=[         0,       1,         0,0,0;
    -274921.63,-73.2915,-274921.63,0,0;
             0,       0,   -0.9597,0,0;
             1,       0,         0,0,0;
             0,       0,         0,1,0];
 B=[0;0.12841;-3.39561e-7;0;0];
 C=[1,0,0,0,0; 0,0,0,1,0; 0,0,0,0,1];
 B1=[0,0;-274921.63,0;0,0;0,-1;0,0];
 C1=[0 0 0 0 1];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);
 [nx,nw]=size(B1); [nz,nx]=size(C1);  
 D11=zeros(nz,nw); D12=zeros(nz,nu); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (TF1): Terrain following model
%       Gershon, Shaked, Yaesh, Tech.-Rep. 2003 (Uni. Tel-Aviv)
%       "Static output feedback of state multiplicative systems with
%       application to terrain following"
% Note: This is not a classical SOF control design --> special ROC
% Q: Is the problem SOF stabilizable too ?
%------------------------------------------------------------------
if(strmatch('TF1',ex,'exact'))
 beta=1/20; eps2=-1e-5;
 A=[-1 0 0 0  0 0 0; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 0 0 0 0 0;
     0 0 0 1 -1 0 0; -0.088 0.0345 0 0 1 -0.0032 0;
     0 0 beta 0 0 0 eps2];
 B=[1 0; 0 0; 0 0; 0 0.09; 0 0; 0 0; 0 0];  
 C=[0 0 1 0 0 0 0; 0 1 0 0 0 0 0; 0 0 0 0 0 0 1;
    0 0 0 0 0 1 0]; 
 B1=[0; 0; 0; 0; 0; 0; -beta];
 C1=zeros(4,7); C1(1,7)=1.e0; C1(2,6)=2.23e0;
 D11=zeros(4,1);
 D12=[zeros(2,2); sqrt(3) 0; 0 sqrt(0.3)];
 D21=[0.04; 0; 0; 0];
end

%------------------------------------------------------------------
% (TF2): Like (TF1) with a different sensor matrix C.
%------------------------------------------------------------------
if(strmatch('TF2',ex,'exact'))
 beta=1/20; eps2=-1e-5;
 A=[-1 0 0 0  0 0 0; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 0 0 0 0 0;
     0 0 0 1 -1 0 0; -0.088 0.0345 0 0 1 -0.0032 0;
     0 0 beta 0 0 0 eps2];
 B=[1 0; 0 0; 0 0; 0 0.09; 0 0; 0 0; 0 0];  
 C=[0 0 1 0 0 0 0; 0 1 0 0 0 0 0; 0 0 0 1 0 0 0]; 
 B1=[0; 0; 0; 0; 0; 0; -beta];
 C1=zeros(4,7); C1(1,7)=1.e0; C1(2,6)=2.23e0;
 D11=zeros(4,1);
 D12=[zeros(2,2); sqrt(3) 0; 0 sqrt(0.3)];
 D21=[0.04; 0; 0];
end

%------------------------------------------------------------------
% (TF3): Another sensor matrix $C$ for the terrain following model.
%------------------------------------------------------------------
if(strmatch('TF3',ex,'exact'))
 beta=1/20; eps2=-1e-5;
 A=[-1 0 0 0  0 0 0; 1 0 0 0 0 0 0; 0 1 0 0 0 0 0; 0 0 0 0 0 0 0;
     0 0 0 1 -1 0 0; -0.088 0.0345 0 0 1 -0.0032 0;
     0 0 beta 0 0 0 eps2];
 B=[1 0; 0 0; 0 0; 0 0.09; 0 0; 0 0; 0 0];  
 C=[0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 0 0 1]; 
 B1=[0; 0; 0; 0; 0; 0; -beta];
 C1=zeros(4,7); C1(1,7)=1.e0; C1(2,6)=2.23e0;
 D11=zeros(4,1);
 D12=[zeros(2,2); sqrt(3) 0; 0 sqrt(0.3)];
 D21=[0.04; 0; 0];
end

%------------------------------------------------------------------
% (PSM): Power system model
%        A. Varga, "Model Reduction Routines for {SLICOT}",
%        NICONET Report 1999-8, p. 32
%        and
%        C. E. Fosha and O. I. Elgerd,"The megawatt-frequency control
%        problem: a new approach via optimal control theory",
%        IEEE Trans. on Power Apparatus and Systems,Vol.89,pp.563-571,1970
%------------------------------------------------------------------
if(strmatch('PSM',ex,'exact'))
  nx=7; nu=2; ny=3;
  A=[-0.04165, 0   , 4.92,-4.92, 0      , 0   , 0   ;
     -5.21   ,-12.5, 0   , 0   , 0      , 0   , 0   ;
      0      , 3.33,-3.33, 0   , 0      , 0   , 0   ;
      0.545  , 0   , 0   , 0   ,-0.545  , 0   , 0   ;
      0      , 0   , 0   , 4.92,-0.04165, 0   , 4.92;
      0      , 0   , 0   , 0   ,-5.21   ,-12.5, 0   ;
      0      , 0   , 0   , 0   , 0      , 3.33,-3.33];
  B=[-4.92 0; 0 0; 0 0; 0 0; 0 -4.92 ; 0 0; 0 0];
  C=[1 0 0 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 1 0 0];
  B1=[0 0; 12.5 0; 0 0; 0 0; 0 0 ; 0 12.5; 0 0];
  C1=[C; zeros(2,7)];
  D12=[zeros(3,2); eye(2)];
  [nz,nx]=size(C1); [nx,nw]=size(B1);
  D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (TL): Transmission Line
%       SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.2
%       W. Draijer, M. Steinbuch, O.H. Bosgra
%       and
%       J.-R. Li and J. White, "Efficient Model Reduction of Inter-
%       connect via Approximate System Gramians", IEEE,
%       0-7803-5832-5/99, 1999
%       Note: Ex'=Ax+Bu; y=Cx given
%------------------------------------------------------------------
if(strmatch('TL',ex,'exact'))
  nx=256; nu=2; ny=2;
  load tl; invE=inv(E);
  A=invE*A; B=invE*B; B1=invE*eye(nx);
  C1=eye(nx); [nx,nw]=size(B1); [nz,nx]=size(C1);
  D12=[zeros(nz-nu,nu);eye(nu)];  
  D11=zeros(nz,nw); D21=zeros(ny,nw);
end

%------------------------------------------------------------------
% (CDP): CD player
%        SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.3
%        W. Draijer, M. Steinbuch, O.H. Bosgra
%        and
%       "Adaptive control of the radial servo system of a compact disc player"
%        Automatica, Vol. 28, No. 3, pp. 455-462, 1992 
%------------------------------------------------------------------
if(strmatch('CDP',ex,'exact'))
  nx=120; nu=2; ny=2;
  load cdp;
  A=full(A);
  B1=1.0e-2*B; C1=[C; zeros(2,nx)]; D12=sqrt(2)*[zeros(2,2);eye(nu)];
  [nz,nx]=size(C1); [nx,nw]=size(B1);
  D21=[0 1; 0.03 2e-5]; D11=zeros(nz,nw);
end

%--------------------------------------------------------------------
% (NN1): L. F. Miller, R. G. Cochran and J. W. Howze, 
%        "Output feedback stabilization of a spectral radius functional", 
%        IJOC, Vol. 27, pp. 455-462, 1978 
%--------------------------------------------------------------------
if(strmatch('NN1',ex,'exact'))
 nx=3;nu=1;ny=2;
 A=[0 1 0;0 0 1;0 13 0];B=[0;0;1];C=[0 5 -1;-1 -1 0];
end

%--------------------------------------------------------------------
% (NN2): Classical example 
%        W. S. Levine and M. Athans, "On the determination of the optimal
%        constant output feedback gains for linear multivariable systems",
%        TOAC, Vol. 15, Nr. 8, pp. 44-48 
%--------------------------------------------------------------------
if(strmatch('NN2',ex,'exact'))
 A=[0,1;-1,0];
 B1=[1,0;0,1]; B=[0;1];
 C1=[1,0;0,0]; C=[0,1];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B); 
 [nz,nx]=size(C1); [ny,nx]=size(C); 
 D11=zeros(nz,nw); 
 D21=zeros(ny,nw);
 D12=[0;1];
end

%--------------------------------------------------------------------
% (NN3): C. W. Scherer, "An Efficient Solution to Multi--Objective
%        Control Problems with LMI Objectives",
%        Delft University of Technology, The Netherlands, 2000
%--------------------------------------------------------------------
if(strmatch('NN3',ex,'exact'))
 A=[0.5 1 1.5 1; -1 3 2.1 2; 1 -1 -0.6 1; -2 2 -1 1];
 B1=[0; 0; 1;0];
 B=[0;0;0;1];
 C1=[1 0 0 0];
 C=[0 0 0 1];
 D12=[0];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
 D21=zeros(ny,nw);
end

%--------------------------------------------------------------------
% (NN4): L. F. Miller, R. G. Cochran and J. W. Howze, 
%        "Output feedback stabilization of a spectral radius functional", 
%        IJOC, Vol. 27, pp. 455-462, 1978 
%--------------------------------------------------------------------
if(strmatch('NN4',ex,'exact'))
 nx=4;nu=2;ny=3;
 A=[0 1 0 0;0 -2.93 -4.75 -0.78;0.086 0 -0.11 -1;0 -0.042 2.59 -0.39];
 B=[0 0;0 -3.91;0.035 0;-2.53 0.31];
 C=[1 0 0 0;0 1 0 0;0 0 1 0];
end

%--------------------------------------------------------------------
% (NN5): Saturn V booster
%        L. F. Miller, R. G. Cochran and J. W. Howze, 
%        "Output feedback stabilization of a spectral radius functional", 
%        IJOC, Vol. 27, pp. 455-462, 1978
%--------------------------------------------------------------------
if(strmatch('NN5',ex,'exact'))
 nx=7;nu=1;ny=2;
 A=[0 1 0 0 0 0 0;0 0 .2 -.65 -.002 2.6 0;-.014 1 -.041 .0002 -.015 -.033 0;
    0 0 0 0 1 0 0;0 0 0 -45 -.13 255 0;0 0 0 0 0 0 1;0 0 0 0 0 -50 -10];
 B=[0;0;0;0;0;0;1];
 C=[1 0 0 0 0 0 0;0 1 0 0 0 0 0];
end

%--------------------------------------------------------------------
% (NN6): H. P. Horisberger and P. R. Belanger, "Solution of the Optimal
%        Constant Output Feedback Problem by Conjugate Gradients",
%        TOAC, Vol. 19, pp. 434-435, 1974   %%%ehemals (HB1)
%--------------------------------------------------------------------
if(strmatch('NN6',ex,'exact'))
 nx=9;nu=1;ny=4;
 A=[0,1,0,0,0,0,0,0,0;0,-20,-4.2,0,4.45,12.5,0,100,0;0,0,0,1,0,0,0,0,0;
    0,4.7,8.35,0,-1.1,0,0,0,0;0,0,0,0,-3.3,0,0,0,0;0,0,0,0,0,0,1,0,0;
    0,10.9,0,0,-2.55,-250,0,0,0;0,0,0,0,0,0,0,0,1;0,5.9,0,0,-1.39,0,0,-3700,0];
 B=[0;0;0;0;3.3;0;0;0;0];
 C=[1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;
     0,0,1,0,0,0.66,0,1.2,0;0,0,0,1,0,0,0.66,0,1.2];    
 B1=sqrt(0.1)*eye(nx); C1=sqrt(10)*eye(nx);
 [nz,nx]=size(C1); [nx,nw]=size(B1); 
 D12=10*[zeros(nz-nu,nu);eye(nu)];
 D11=zeros(nz,nw); D21=zeros(ny,nw);  
end

%--------------------------------------------------------------------
% (NN7): like (NN6) with changed B1, C1, D11, D12 and D21 %%%ehemals (HB2)
%--------------------------------------------------------------------
if(strmatch('NN7',ex,'exact'))
 nx=9;nu=1;ny=4;
 A=[0,1,0,0,0,0,0,0,0;0,-20,-4.2,0,4.45,12.5,0,100,0;0,0,0,1,0,0,0,0,0;
    0,4.7,8.35,0,-1.1,0,0,0,0;0,0,0,0,-3.3,0,0,0,0;0,0,0,0,0,0,1,0,0;
    0,10.9,0,0,-2.55,-250,0,0,0;0,0,0,0,0,0,0,0,1;0,5.9,0,0,-1.39,0,0,-3700,0];
 B=[0;0;0;0;3.3;0;0;0;0];
 C=[1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;
    0,0,1,0,0,0.66,0,1.2,0;0,0,0,1,0,0,0.66,0,1.2];    
 B1=[1.45e-1, 4.78e-1, 0, 0, 0; 0,0, 0,-1, 0; 5.23e-2,0, 1, 0, 0;
     0,0, 0, 0, 1; 0, 0, 0, 0, 0; 0, 5.98e-1, 0, 1, 0;
     1,0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 1, 0, 0]; 
 C1=[1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0];
 D12=[0;0;1];
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);  
end

%--------------------------------------------------------------------
% (NN8)
%--------------------------------------------------------------------
if(strmatch('NN8',ex,'exact'))
 nx=3;nu=2;ny=2;
 A=[-0.2 0.1 1;-0.05 0 0;0 0 -1];
 B=[0 1; 0 0.7;1 0];
 C=[1 0 0; 0 1 0];
end

%--------------------------------------------------------------------
% (NN9): B. M. Chen, "H_inf Control and Its Applications",
%        Springer-Verlag, "Lecture Notes in Control and Information Sciences",
%        Vol.235, 1998; p. 110
%--------------------------------------------------------------------
if(strmatch('NN9',ex,'exact'))
 nx=5; nu=3; ny=2;
 A=[1 1 1 0 1;0 1 0 0 1;0 1 1 0 1;1 1 1 1 1;1 1 1 1 0];
 B=[0 0 0;0 0 0;1 0 0;0 0 1;0 1 0];
 C=[0 -2 -3 -2 -1;1 2 3 2 1];
 B1=[5 1;0 0;0 0;2 3;1 4];
 C1=[0 0 1 0 0;0 0 0 0 1;0 1 0 0 0;0 0 1 0 0];
 D12=[1 0 0;0 0 0;0 0 0;0 0 0];
 D21=[1 0;0 0];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw);
end

%--------------------------------------------------------------------
% (NN10): X. A. Wang, "Grassmannian, Central Projection, and Output
%         Feedback Pole Assignment of Linear Systems", TOAC, Vol. 41, 
%         Nr. 6, pp. 786-794, 1996; Example 3.7
%--------------------------------------------------------------------
if(strmatch('NN10',ex,'exact'))
 A=[0,-1, 0, 0, 0, 0, 0, 1;
    1, 2, 0, 0, 1, 0, 0,-2;
    0,-1, 0, 0, 5, 0, 0, 0;
    0, 0, 1, 0,-7, 0, 0,-2;
    0,-1, 0, 1, 4, 0, 0, 2;
    0,-2, 0, 0, 2, 0, 0, 3;
    0, 0, 0, 0,-1, 1, 0,-2;
    0,-1, 0, 0, 1, 0, 1,-1];
 B=[ 0, 1 ,2;1, 0, 1;-1,-1,-3;1, 0, 1;
     0, 2, 4;2, 1, 5;-1, 1, 1;1,-1,-1];
 C=[0,1,0,0,-2,0,0,0;0,0,0,0, 1,0,0,0;0,0,0,0, 0,0,0,1];
 B1=zeros(8,3);
 C1=[1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw);
 D12=zeros(nz,nu); D21=zeros(ny,nw);            
end

%--------------------------------------------------------------------
% (NN11): P. Apkarian and H. D. Tuan, "Robust Conrol via Concave
%         Minimization, Local and Global Algorithms", TOAC, Vol. 45, 
%         Nr. 2, pp. 299-305, 2000 
%--------------------------------------------------------------------
if(strmatch('NN11',ex,'exact'))
 A1=[-101 -99.9 0 0 0 0 0 0;
      0 -101 0 0 0 0 0 0;
      0 0 -101 -99.9 0 0 0 0;
        0 0 0 -101 0 0 0 0;
          0 0 0 0 -1 0 0 0;
          0 0 0 0 0 -1 0 0;
          0 0 0 0 0 0 -1 0;
          0 0 0 0 0 0 0 -1;
      0 0 0 0 0 0 427.098 -46.8341;
      0 0 0 0 0 0 232.0719 120.4649;
      0 0 0 0 0 0 -764.2456 85.4154;
      0 0 0 0 0 0 166.827 -264.7739;
      0 0 0 0 0.3162 0 0 0;
      0 0 0 0 -0.125 0 0 0;
      0 0 0 0 0 0.3162 0 0;
      0 0 0 0 0 -0.125 0 0];
 A2=[0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
   -1 0 0.4271 -0.0468 0 0 0 0;
    0 -1 0.2321 0.1205 0 0 0 0;
    0 0 -1.7642 0.0854 0 0 0 0;
    0 0 0.1668 -1.2648 0 0 0 0;
    0 0 0 0 -1.1 -0.0759 0 0;
    0 0 0 0 0 -1 0 0;
    0 0 0 0 0 0 -1.1 -0.0759;
    0 0 0 0 0 0 0 -1];
 A=[A1 A2];
 B=[0 -9.995 0; 0.199 -9.995 0; 0.211 0 -9.995; -0.233 0 -9.995; 0 0 0; 0 0 0; 0 0 0; 0 0 0;
   0 2.7173 1.4274; 0 1.4274 2.8382; 0 -4.7909 -2.6032; 0 1.0261 -2.6393; 0.11 0 0; 0 0 0;
   0 0 0; 0.01 0 0];
 C=[0 0 0 0 0 0   0       0   0 0    0     0    -0.3162 0    0    0;
    0 0 0 0 0 0   0       0   0 0    0     0        0   0 -0.3162 0;
    0 0 0 0 0 0 1.5564 3.4834 0 0 0.0016 0.0035     0   0    0    0;
    0 0 0 0 0 0   0       0   0 0    0     0    -0.4743 0    0    0;
    0 0 0 0 0 0   0       0   0 0    0     0        0   0 -0.3479 0];
 B1=[0 -0.001 0;0 -0.001 0;0 0 -0.001;0 0 -0.001;
     0 0 0 ; 0 0 0; 0 0 0; 0 0 0;0.1787 0.0003 0.0001;
     -0.8364 0.0001 0.0003; 0.0818 -0.0005 -0.0003;
     0.3577 0.0001 -0.0003;0 -0.3162 0;0 0.125 0;-0.3162 0 0;0 0 0.125];
 C1=[0 0 0 0 0 0 1.5564 3.4834 0 0 0.0016 0.0035 0 0 0 0;
     0 0 0 0 0 0      0      0 0 0      0      0 -0.4743 0 0 0;
     0 0 0 0 0 0      0      0 0 0      0      0 0 0 -0.3479 0];         
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw);
 D12=zeros(nz,nu); D21=zeros(ny,nw);     
end

%--------------------------------------------------------------------
% (NN12): J. Rosenthal and X. A. Wang, "Output Feedback Pole Placement
%         with Dynamic Compensators", TOAC, Vol. 41, Nr. 6, 
%         pp.830-843, 1996, Example 3.14
%--------------------------------------------------------------------
if(strmatch('NN12',ex,'exact'))
  A=[0,0,0,0,0,0;1,0,0,0,0,-1;0,1,0,0,0,0;
     0,0,0,0,0,0;0,0,0,1,0,0;0,0,-1,0,1,0];
  B=-[1,3;0,0;0,-1;0,1;0,1;0,0];
  C=[0,0,1,0,0,0;0,0,0,0,0,1];  
end

%--------------------------------------------------------------------
% (NN13): D.-W. Gu, P. Hr. Petkov and M. M. Konstantinov,
%         "H_inf and H_2 Optimization Toolbox in SLICOT",
%         SLICOT Working note 1999-12, p. 15
%         available via ftp: wgs.esat.kuleuven.ac.be/
%                            pub/WGS/REPORTS/SLWN1999-12.ps.Z
%--------------------------------------------------------------------
if(strmatch('NN13',ex,'exact'))
 A=[-1 0 4 5 -3 -2;-2 4 -7 -2 0 3; -6 9 -5 0 2 -1;-8 4 7 -1 -3 0;
    2 5 8 -9 1 -4; 3 -5 8 0 2 -6];
 B=[1 0;-5 2; 7 -2; 1 -2;0 5;-6 -2];
 C=[9 -3 4 0 3 7;0 1 -2 1 -6 -2];
 B1=[-3 -4 -2;2 0 1;-5 -7 0;4 -6 1;-3 9 -8;1 -2 3];
 C1=[1 -1 2 -4 0 -3;-3 0 5 -1 1 1;-7 5 0 -8 2 -2];
 D11=[1 -2 -3;0 4 0;5 -3 -4];
 D12=[0 0;1 0;0 1]; D21=[0 1 0;0 0 1];  
end


%--------------------------------------------------------------------
% (NN14): D.-W. Gu, P. Hr. Petkov and M. M. Konstantinov,
%         "H_inf and H_2 Optimization Toolbox in SLICOT",
%         SLICOT Working note 1999-12, p. 19
%         same as NN13, diff. D11, D12, D21
%         available via ftp: wgs.esat.kuleuven.ac.be/
%                            pub/WGS/REPORTS/SLWN1999-12.ps.Z
%--------------------------------------------------------------------
if(strmatch('NN14',ex,'exact'))
  A=[-1 0 4 5 -3 -2;-2 4 -7 -2 0 3; -6 9 -5 0 2 -1;-8 4 7 -1 -3 0;
    2 5 8 -9 1 -4; 3 -5 8 0 2 -6];
 B=[1 0;-5 2; 7 -2; 1 -2;0 5;-6 -2];
 C=[9 -3 4 0 3 7;0 1 -2 1 -6 -2];
 B1=[-3 -4 -2;2 0 1;-5 -7 0;4 -6 1;-3 9 -8;1 -2 3];
 C1=[1 -1 2 -4 0 -3;-3 0 5 -1 1 1;-7 5 0 -8 2 -2];
 D11=zeros(3,3);
 D12=[-4 -1;1 0;0 1]; D21=[3 1 0;-2 0 1];
end

%--------------------------------------------------------------------
% (NN15): Space backpack model
%         P. L. D. Peres and J. C. Geromel,
%         "An Alternate Numerical Solution to the Linear Quadratic Problem",
%         TOAC, Vol. 39, Nr. 1, pp. 198-202, 1994
%--------------------------------------------------------------------
if(strmatch('NN15',ex,'exact'))
 A=[0,1,0;-79.285,-0.113,0;28.564,0.041,0];
 %B1=eye(size(A));
 B1=[0;0.041;-0.03];
 B=[0,0;0.041,-0.0047;-0.03,-0.0016];
 C1=[0,0,1;1,0,0;0,0,0;0,0,0];
 C=[0,0,1;1,0,0];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);  
 [nz,nx]=size(C1); [ny,nx]=size(C);    
 D11=zeros(nz,nw);    
 D21=zeros(ny,nw);
 D12=[0,0;0,0;0.1,0;0,0.1];
end

%------------------------------------------------------------------
% (NN16): Application for an large space structure 
%         A. J. Calise and D. D. Moerder, "Optimal Output Feedback
%         Design of Systems with Ill-conditioned Dynamics",
%         AUTO, Vol. 21, Nr. 3, pp. 271-276, 1985
%------------------------------------------------------------------
if(strmatch('NN16',ex,'exact'))
  nx=8;nu=4;ny=4;
  A=[ 0    , 1 , 0      , 0 , 0    , 0 , 0    , 0;
     -0.42 , 0 , 0      , 0 , 0    , 0 , 0    , 0;
      0    , 0 , 0      , 1 , 0    , 0 , 0    , 0;
      0    , 0 ,-0.1849 , 0 , 0    , 0 , 0    , 0;
      0    , 0 , 0      , 0 , 0    , 1 , 0    , 0;
      0    , 0 , 0      , 0 ,-4.41 , 0 , 0    , 0;
      0    , 0 , 0      , 0 , 0    , 0 , 0    , 1;
      0    , 0 , 0      , 0 , 0    , 0 ,-4.84 , 0];
  B=[ 0    , 0    , 0    , 0  ;
     -0.92 ,-1.4  , 0.92 ,-1.4;
      0    , 0    , 0    , 0  ;
      0.65 , 1.6  , 0.65 ,-1.6;
      0    , 0    , 0    , 0  ;
      1.4  ,-1    , 1.4  , 1  ;
      0    , 0    , 0    , 0  ;
      2    ,-0.8  ,-2    ,-0.8];
  C=[ 0 ,-1.8 , 0 , 1.3 , 0 , 2.9 , 0 , 4.1;
      0 ,-2.7 , 0 , 3.2 , 0 ,-2.1 , 0 ,-1.6;
      0 , 1.8 , 0 , 1.3 , 0 , 2.9 , 0 ,-4.1;
      0 ,-2.7 , 0 ,-3.2 , 0 , 2.1 , 0 ,-1.6];
  C1=zeros(4,8); C1(1,1)=0.065;C1(2,2)=0.065;
  D12=eye(4);
  B1=eye(nx);
  D11=zeros(4,8);
  D21=zeros(4,8);
end

%--------------------------------------------------------------------
% (NN17): Rank-deficient matrix D12, P. Gahinet and A. J. Laub,
%         "Numerically reliable computation of optimal performance in
%         singular H_inf control", SIOPT, Vol. 35, Nr. 5,
%         pp. 1690-1710, 1997
%--------------------------------------------------------------------
if(strmatch('NN17',ex,'exact'))
 A=[0,-1,2; 1,-2,3; 0,1,0];
 B1=[1; -1; 0];
 B=[1,0; 0,0; 0,-1];
 C1=[1,0,1; 1,0,1];
 C=[1 0 0];   
 D11=zeros(2,1);    
 D21=zeros(1,1); 
 D12=[0,1; 0,0];
end

%------------------------------------------------------------------
% (NN18): SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.5
%         W. Draijer, M. Steinbuch, O.H. Bosgra
%         and
%         T. Penzl, "Algorithms for Model Reduction of Large
%         Dynamical Systems", Technische Universitaet Chemnitz, SFB 393,
%         Preprint SFB393/99-40, 1999
%------------------------------------------------------------------
if(strmatch('NN18',ex,'exact'))
  nx=1006; nu=1; ny=1;
  d=[1;1;1;1;1;1;[1:1000]'];
  sd=[  -100;0;-200;0;-400;zeros(1001,1)];
  ud=[0; 100;0; 200;0; 400;zeros(1000,1)];
  A=spdiags([sd,-d,ud],-1:1,1006,1006);
  C =[10*ones(1,6),ones(1,1000)]; B=C';
  C1=[C; zeros(1,1006)];
  B1=0.01*B;
  D11=zeros(2,1);
  D12=[0;1];
  D21=zeros(1,1);
end

%------------------------------------------------------------------
% (HF2D1): 2D heat flow example, HF2D_JetControlExamples
%          thermal properties of copper, sd:=0.3825
%------------------------------------------------------------------
if(strmatch('HF2D1',ex,'exact'))
 nx=3796; nu=2; ny=3;
 load hf2d1;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D1_M316): medium scale version of (HF2D1), sd:=0.615
%------------------------------------------------------------------
if(strmatch('HF2D1_M316',ex,'exact'))
 nx=316; nu=2; ny=3;
 load hf2d1_m316;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D1_M541): medium scale version of (HF2D1), sd:=0.615
%------------------------------------------------------------------
if(strmatch('HF2D1_M541',ex,'exact'))
 nx=541; nu=2; ny=3;
 load hf2d1_m541;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D2): 2D heat flow example, HF2D_JetControlExamples
%          thermal properties of platinum, sd:=0.5325
%------------------------------------------------------------------
if(strmatch('HF2D2',ex,'exact'))
 nx=3796; nu=2; ny=3;
 load hf2d2;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D2_M316): medium scale version of (HF2D2), sd:=0.72
%------------------------------------------------------------------
if(strmatch('HF2D2_M316',ex,'exact'))
 nx=316; nu=2; ny=3;
 load hf2d2_m316;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D2_M541): medium scale version of (HF2D2), sd:=0.72
%------------------------------------------------------------------
if(strmatch('HF2D2_M541',ex,'exact'))
 nx=541; nu=2; ny=3;
 load hf2d2_m541;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (HF2D3): 2D heat flow example, rectangular domain
%          thermal properties of copper
%------------------------------------------------------------------
if(strmatch('HF2D3',ex,'exact'))
 nx=4489; nu=2; ny=4;
 load hf2d3;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D4): 2D heat flow example, rectangular domain 
%          thermal properties of platinum
%------------------------------------------------------------------
if(strmatch('HF2D4',ex,'exact'))
 nx=2025; nu=2; ny=4;
 load hf2d4;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);  
end

%------------------------------------------------------------------
% (HF2D5): 2D heat flow example, rectangular domain with perturbation
%          operator added
%          thermal properties of copper, sd=0.3825
%------------------------------------------------------------------
if(strmatch('HF2D5',ex,'exact'))
 nx=4489; nu=2; ny=4;
 load hf2d5;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D5_M289): medium scale version of (HF2D5), sd:=?
%------------------------------------------------------------------
if(strmatch('HF2D5_M289',ex,'exact'))
 nx=289; nu=2; ny=4;
 load hf2d5_m289;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D5_M529): medium scale version of (HF2D5), sd:=?
%------------------------------------------------------------------
if(strmatch('HF2D5_M529',ex,'exact'))
 nx=529; nu=2; ny=4;
 load hf2d5_m529;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (HF2D6): 2D heat flow example, rectangular domain with perturbation
%          operator added
%          thermal properties of of platinum, sd:=1.725
%------------------------------------------------------------------
if(strmatch('HF2D6',ex,'exact'))
 nx=2025; nu=2; ny=4;
 load hf2d6;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D6_M289): medium scale version of (HF2D6), sd:=0.9825
%------------------------------------------------------------------
if(strmatch('HF2D6_M289',ex,'exact'))
 nx=289; nu=2; ny=4;
 load hf2d6_m289;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D6_M529): medium scale version of (HF2D6), sd:=1.8825
%------------------------------------------------------------------
if(strmatch('HF2D6_M529',ex,'exact'))
 nx=529; nu=2; ny=4;
 load hf2d6_m529;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (HF2D7): 2D heat flow example, rectangular domain with perturbation
%          operator added and extended to a nonlinear model
%          thermal properties of copper
%------------------------------------------------------------------
if(strmatch('HF2D7',ex,'exact'))
 nx=4489; nu=2; ny=4;
 load hf2d7;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D8): 2D heat flow example, rectangular domain with perturbation
%          operator added and extended to a nonlinear model
%          thermal properties of of platinum
%------------------------------------------------------------------
if(strmatch('HF2D8',ex,'exact'))
 nx=2025; nu=2; ny=4;
 load hf2d8;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);  
end

%------------------------------------------------------------------
% (HF2D9): 2D heat flow example, Distributed control of a perturbed
%          linear heat flow model 
%------------------------------------------------------------------
if(strmatch('HF2D9',ex,'exact'))
 nx=3481; nu=2; ny=2;
 load hf2d9;
 [nx,nw]=size(B1); [nz,nx]=size(C1);  
 D11=sparse(nz,nw); D21=sparse(ny,nw); 
end

%------------------------------------------------------------------
% (HF2D9_M256): medium scale version of (HF2D9), sd:=0.4829 
%------------------------------------------------------------------
if(strmatch('HF2D9_M256',ex,'exact'))
 nx=256; nu=2; ny=2;
 load hf2d9_m256;
 [nx,nw]=size(B1); [nz,nx]=size(C1);  
 D11=sparse(nz,nw); D21=sparse(ny,nw); 
end

%------------------------------------------------------------------
% (HF2D9_M484): medium scale version of (HF2D9), sd:=? 
%------------------------------------------------------------------
if(strmatch('HF2D9_M484',ex,'exact'))
 nx=484; nu=2; ny=2;
 load hf2d9_m484;
 [nx,nw]=size(B1); [nz,nx]=size(C1);  
 D11=sparse(nz,nw); D21=sparse(ny,nw); 
end



%------------------------------------------------------------------
% (HF2D10): reduced version of (HF2D1)
%------------------------------------------------------------------
if(strmatch('HF2D10',ex,'exact'))
 nx=5; nu=2; ny=3;
 load hf2d10;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw); 
end

%------------------------------------------------------------------
% (HF2D11): reduced version of (HF2D2)
%------------------------------------------------------------------
if(strmatch('HF2D11',ex,'exact'))
 nx=5; nu=2; ny=3;
 load hf2d11;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);  
end

%------------------------------------------------------------------
% (HF2D12): reduced version of (HF2D3)
%------------------------------------------------------------------
if(strmatch('HF2D12',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d12;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw); 
end

%------------------------------------------------------------------
% (HF2D13): reduced version of (HF2D4)
%------------------------------------------------------------------
if(strmatch('HF2D13',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d13;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);   
end

%------------------------------------------------------------------
% (HF2D14): reduced version of (HF2D5)
%------------------------------------------------------------------
if(strmatch('HF2D14',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d14;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);  
end

%------------------------------------------------------------------
% (HF2D15): reduced version of (HF2D6)
%------------------------------------------------------------------
if(strmatch('HF2D15',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d15;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);  
end

%------------------------------------------------------------------
% (HF2D16): reduced version of (HF2D7)
%------------------------------------------------------------------
if(strmatch('HF2D16',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d16;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);   
end

%------------------------------------------------------------------
% (HF2D17): reduced version of (HF2D8)
%------------------------------------------------------------------
if(strmatch('HF2D17',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d17; 
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);   
end

%------------------------------------------------------------------
% (HF2D18): reduced version of (HF2D9)
%------------------------------------------------------------------
if(strmatch('HF2D18',ex,'exact'))
 nx=5; nu=2; ny=2;
 load hf2d18;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=zeros(nz,nw); D21=zeros(ny,nw);   
end

%------------------------------------------------------------------
% (HF2D_CD1): 2D convection diffusion equation (linear) with 
%             distributed control input on a rectangular domain
%             thermal properties of copper, 
%             eps1=0.29555, eps2=2.9555,
%             v_t=kappa(v_{xx}+v_{yy})-eps1*(v_x*v_y)+eps2*v
%                 + Sum_i( u_i b_i)
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_CD1',ex,'exact'))
 nx=3600; nu=2; ny=2;
 load hf2d_cd1;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD1_M256): medium scale version of (HF2D_CD1), 
%                  eps1=0.29555, eps2=2.9555
%------------------------------------------------------------------
if(strmatch('HF2D_CD1_M256',ex,'exact'))
 nx=256; nu=2; ny=2;
 load hf2d_cd1_m256;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD1_M484): medium scale version of (HF2D_CD1), 
%                  eps1=0.29555, eps2=2.9555
%------------------------------------------------------------------
if(strmatch('HF2D_CD1_M484',ex,'exact'))
 nx=484; nu=2; ny=2;
 load hf2d_cd1_m484;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD2): 2D convection diffusion equation (linear) with 
%             distributed control input on a rectangular domain
%             thermal properties of platin, 
%             eps1=0.29555, eps2=2.1755,
%             v_t=kappa(v_{xx}+v_{yy})-eps1*(v_x*v_y)+eps2*v
%                 + Sum_i( u_i b_i)
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_CD2',ex,'exact'))
 nx=3600; nu=2; ny=2;
 load hf2d_cd2;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD2_M256): medium scale version of (HF2D_CD2), 
%                  eps1=0.29555, eps2=2.1755
%------------------------------------------------------------------
if(strmatch('HF2D_CD2_M256',ex,'exact'))
 nx=256; nu=2; ny=2;
 load hf2d_cd2_m256;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD2_M484): medium scale version of (HF2D_CD2), 
%                  eps1=0.29555, eps2=2.1755
%------------------------------------------------------------------
if(strmatch('HF2D_CD2_M484',ex,'exact'))
 nx=484; nu=2; ny=2;
 load hf2d_cd2_m484;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD3): 2D convection diffusion equation (linear) with 
%             distributed control input on a rectangular domain
%             thermal properties of yttrium, 
%             eps1=0.29555, eps2=4.35,
%             v_t=kappa(v_{xx}+v_{yy})-eps1*(v_x*v_y)+eps2*v
%                 + Sum_i( u_i b_i)
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_CD3',ex,'exact'))
 nx=4096; nu=2; ny=2;
 load hf2d_cd3;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD3_M324): medium scale version of (HF2D_CD3), 
%                  eps1=0.29555, eps2=4.35
%------------------------------------------------------------------
if(strmatch('HF2D_CD3_M324',ex,'exact'))
 nx=324; nu=2; ny=2;
 load hf2d_cd3_m324;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD3_M576): medium scale version of (HF2D_CD3), 
%                  eps1=0.29555, eps2=4.35
%------------------------------------------------------------------
if(strmatch('HF2D_CD3_M576',ex,'exact'))
 nx=576; nu=2; ny=2;
 load hf2d_cd3_m576;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (HF2D_CD4): POD version of (HF2D_CD1), dense 
%             eps1=0.29555, eps2=2.9555
%------------------------------------------------------------------
if(strmatch('HF2D_CD4',ex,'exact'))
 nx=7; nu=2; ny=2;
 load hf2d_cd4;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD5): POD version of (HF2D_CD2), dense 
%             eps1=0.29555, eps2=2.1755
%------------------------------------------------------------------
if(strmatch('HF2D_CD5',ex,'exact'))
 nx=7; nu=2; ny=2;
 load hf2d_cd5;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_CD6): POD version of (HF2D_CD3), dense 
%             eps1=0.29555, eps2=4.35
%------------------------------------------------------------------
if(strmatch('HF2D_CD6',ex,'exact'))
 nx=7; nu=2; ny=2;
 load hf2d_cd6;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS1): 2D instable heat equation (nonline./linear) with 
%             boundary control input on a rectangular domain
%             thermal properties of copper, 
%             eps2=?, eps1=0 for linear model, =1 for nonlin.
%             v_t=kappa(v_{xx}+v_{yy})+eps2*v
%                 -eps1*v^3 + boundary control input
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_IS1',ex,'exact'))
 nx=4489; nu=2; ny=4;
 load hf2d_is1;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS1_M361): medium scale version of (HF2D_IS1), 
%                  eps2=0.735
%------------------------------------------------------------------
if(strmatch('HF2D_IS1_M361',ex,'exact'))
 nx=361; nu=2; ny=4;
 load hf2d_is1_m361;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS1_M529): medium scale version of (HF2D_IS1), 
%                  eps2=0.635
%------------------------------------------------------------------
if(strmatch('HF2D_IS1_M529',ex,'exact'))
 nx=529; nu=2; ny=4;
 load hf2d_is1_m529;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS2): 2D instable heat equation (nonline./linear) with 
%             boundary control input on a rectangular domain
%             thermal properties of platinum, 
%             eps2=1.9, eps1=0 for linear model, =1 for nonlin.
%             v_t=kappa(v_{xx}+v_{yy})+eps2*v
%                 -eps1*v^3 + boundary control input
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_IS2',ex,'exact'))
 nx=4489; nu=2; ny=4;
 load hf2d_is2;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS2_M361): medium scale version of (HF2D_IS2), 
%                  eps2=1.9
%------------------------------------------------------------------
if(strmatch('HF2D_IS2_M361',ex,'exact'))
 nx=361; nu=2; ny=4;
 load hf2d_is2_m361;
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS2_M529): medium scale version of (HF2D_IS2), 
%                  eps2=1.9
%------------------------------------------------------------------
if(strmatch('HF2D_IS2_M529',ex,'exact'))
 nx=529; nu=2; ny=4;
 load hf2d_is2_m529; 
 invE=inv(E);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS3): 2D instable heat equation (nonline./linear) with 
%             distributed control input on a rectangular domain
%             thermal properties of copper, 
%             eps2=2.5, eps1=0 for linear model, =1 for nonlin.
%             v_t=kappa(v_{xx}+v_{yy})+eps2*v
%                 -eps1*v^3 + distributed control input
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_IS3',ex,'exact'))
 nx=3600; nu=2; ny=2;
 load hf2d_is3;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS3_M256): medium scale version of (HF2D_IS3), 
%                  eps2=2.5
%------------------------------------------------------------------
if(strmatch('HF2D_IS3_M256',ex,'exact'))
 nx=256; nu=2; ny=2;
 load hf2d_is3_m256;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS3_M484): medium scale version of (HF2D_IS3), 
%                  eps2=2.5
%------------------------------------------------------------------
if(strmatch('HF2D_IS3_M484',ex,'exact'))
 nx=484; nu=2; ny=2;
 load hf2d_is3_m484;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS4): 2D instable heat equation (nonline./linear) with 
%             distributed control input on a rectangular domain
%             thermal properties of platinum, 
%             eps2=1.0, eps1=0 for linear model, =1 for nonlin.
%             v_t=kappa(v_{xx}+v_{yy})+eps2*v
%                 -eps1*v^3 + distributed control input
%             (see Leibfritz, Volkwein: 
%              "Numerical feedback controller design for PDE systems 
%               using model reduction: techniques and case studies"
%               TR 2004, submitted)
%------------------------------------------------------------------
if(strmatch('HF2D_IS4',ex,'exact'))
 nx=3600; nu=2; ny=2;
 load hf2d_is4;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS4_M256): medium scale version of (HF2D_IS4), 
%                  eps2=1.0
%------------------------------------------------------------------
if(strmatch('HF2D_IS4_M256',ex,'exact'))
 nx=256; nu=2; ny=2;
 load hf2d_is4_m256;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS4_M484): medium scale version of (HF2D_IS4), 
%                  eps2=1.0
%------------------------------------------------------------------
if(strmatch('HF2D_IS4_M484',ex,'exact'))
 nx=484; nu=2; ny=2;
 load hf2d_is4_m484;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (HF2D_IS5): POD version of (HF2D_IS1), dense 
%------------------------------------------------------------------
if(strmatch('HF2D_IS5',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d_is5; invE=inv(E_pod);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS6): POD version of (HF2D_IS2), dense 
%------------------------------------------------------------------
if(strmatch('HF2D_IS6',ex,'exact'))
 nx=5; nu=2; ny=4;
 load hf2d_is6; invE=inv(E_pod);
 A=invE*A; B1=invE*B1; B=invE*B;
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS7): POD version of (HF2D_IS3), dense 
%------------------------------------------------------------------
if(strmatch('HF2D_IS7',ex,'exact'))
 nx=5; nu=2; ny=2;
 load hf2d_is7; 
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end

%------------------------------------------------------------------
% (HF2D_IS8): POD version of (HF2D_IS4), dense 
%------------------------------------------------------------------
if(strmatch('HF2D_IS8',ex,'exact'))
 nx=5; nu=2; ny=2;
 load hf2d_is8; 
 [nx,nw]=size(B1); [nz,nx]=size(C1); 
 D11=sparse(nz,nw); D21=sparse(ny,nw);
end


%------------------------------------------------------------------
% (CM1): Linear cable mass problem of order 20
%        J. A. Burns and B.B. King,
%        "A reduced bases approach to the design of low order
%        feedback controllers for nonlinear continuous systems",
%        ICAM Virginia Polytechnic Institute and State University,
%        Blacksburg
%        Note: System matrix A is Hurwitz, but
%              max. real part of eig(A) is close to zero
%------------------------------------------------------------------
if(strmatch('CM1',ex,'exact'))
  nx=40; nu=1; ny=2;
 load cm1;    
end

%------------------------------------------------------------------
% (CM2): To (CM1) corresponding example of order 60.
%------------------------------------------------------------------
if(strmatch('CM2',ex,'exact'))
 nx=60; nu=1; ny=2;
 load cm2;    
end

%------------------------------------------------------------------
% (CM3): To (CM1) corresponding example of order 120.
%------------------------------------------------------------------
if(strmatch('CM3',ex,'exact'))
 nx=120; nu=1; ny=2;
 load cm3;    
end

%------------------------------------------------------------------
% (CM4): To (CM1) corresponding example of order 240.
%------------------------------------------------------------------
if(strmatch('CM4',ex,'exact'))
 nx=240; nu=1; ny=2;
 load cm4;    
end

%------------------------------------------------------------------
% (CM5): To (CM1) corresponding example of order 480.
%------------------------------------------------------------------
if(strmatch('CM5',ex,'exact'))
 nx=480; nu=1; ny=2;
 load cm5;    
end

%------------------------------------------------------------------
% (CM6): To (CM1) corresponding example of order 960.
%------------------------------------------------------------------
if(strmatch('CM6',ex,'exact'))
 nx=960; nu=1; ny=2;
 load cm6;    
end

%------------------------------------------------------------------
% (CM1_IS): Perturbetd linear cable mass problem of order 20
%        J. A. Burns and B.B. King,
%        "A reduced bases approach to the design of low order
%        feedback controllers for nonlinear continuous systems",
%        ICAM Virginia Polytechnic Institute and State University,
%        Blacksburg
%        Note: (CM1_IS) instable version of (CM1), i.e. 
%              A of (CM1) is redefined by A:=(A+dA),
%              where dA:=[0 | 0; 0 | ds*eye(0.5*nx)] such that
%              the system matrix A:=(A+dA) is not Hurwitz!
%              Max. real part of eig(A) > 0,
%              sd:=9.5825e-003
%------------------------------------------------------------------
if(strmatch('CM1_IS',ex,'exact'))
  nx=40; nu=1; ny=2;
 load cm1_is;    
end

%------------------------------------------------------------------
% (CM2_IS): To (CM1_IS) corresponding example of order 60,
%           sd:=9.5825e-003
%------------------------------------------------------------------
if(strmatch('CM2_IS',ex,'exact'))
 nx=60; nu=1; ny=2;
 load cm2_is;    
end

%------------------------------------------------------------------
% (CM3_IS): To (CM1_IS) corresponding example of order 120,
%           sd:=9.5825e-003
%------------------------------------------------------------------
if(strmatch('CM3_IS',ex,'exact'))
 nx=120; nu=1; ny=2;
 load cm3_is;    
end

%------------------------------------------------------------------
% (CM4_IS): To (CM1_IS) corresponding example of order 240,
%           sd:=4.79125e-002
%------------------------------------------------------------------
if(strmatch('CM4_IS',ex,'exact'))
 nx=240; nu=1; ny=2;
 load cm4_is;    
end

%------------------------------------------------------------------
% (CM5_IS): To (CM1_IS) corresponding example of order 480,
%           sd:=9.5825e-002
%------------------------------------------------------------------
if(strmatch('CM5_IS',ex,'exact'))
 nx=480; nu=1; ny=2;
 load cm5_is;    
end

%------------------------------------------------------------------
% (CM6_IS): To (CM1_IS) corresponding example of order 960,
%           sd:=1.72485e-001
%------------------------------------------------------------------
if(strmatch('CM6_IS',ex,'exact'))
 nx=960; nu=1; ny=2;
 load cm6_is;    
end

%------------------------------------------------------------------
% (TMD): 2-degree-of-freedom (2DOF) tuned-mass damper (TMD);
%        L. Zuo and S. A. Nayfeh, "The Multi-Degree-of-Freedom Tuned-
%        Mass-Damper for Suppression of Single-Mode Vibration Under
%        Random and Harmonic Excitation", Dept. of Mech. Engineering,
%        Massachusetts Institute of Technology, Cambridge, Massachusetts 02139,
%        Draft for 2003 ASME Design Engineering Technical Conferences",
%------------------------------------------------------------------
if(strmatch('TMD',ex,'exact'))
md=2; ms=5; cs=0.05; ks=3.0; d=0.6; rho=0.4;
 Mp=[ 0.5*[md md; -md md] zeros(2,1); zeros(1,2) ms ];
 Cp=zeros(3,3); Cp(3,3)=cs;
 Kp=zeros(3,3); Kp(3,3)=ks;
 Bv=[0;0;cs]; Bp=[0;0;ks]; Bd=[0;0;1];
 Bu=[-1 -1; d*d/(rho*rho) -d*d/(rho*rho); 1 1];
 
 A = [zeros(3,3) eye(3); -inv(Mp)*Kp -inv(Mp)*Cp];
 B1= [inv(Mp)*Bv; inv(Mp)*(Bp-Cp*inv(Mp)*Bv)];
 B = [zeros(3,2); inv(Mp)*Bu];
 C1= [0 0 1 0 0 0; zeros(2,6)];
 C = [1 0 -1 0 0 0; 0 0 0 1 0 -1; 0 1 -1 0 0 0; 0 0 0 0 1 -1];
 D11=zeros(3,1);
 D12=[0 0; eye(2)];
 D21=C*[zeros(3,1); inv(Mp)*Bv];
end

%------------------------------------------------------------------
% (FS): Flexible satellite
%       Buschek, Calise: "mu-controllers: mixed and fixed"
%       Proc. AIAA Guidance, Nav. and Control Conf. Baltimore, 1995
%------------------------------------------------------------------
if(strmatch('FS',ex,'exact'))
 M=[77511 248.1; 248.1 1]; D=[0 0; 0 0.002288]; K=[0 0; 0 0.098696];
 Bh=[1;0];

 A=[zeros(2,2)   eye(2)  zeros(2,1);
    -inv(M)*K  -inv(M)*D zeros(2,1);
    1 0 0 0 0];
 B=[0;0;inv(M)*Bh;0];
 C=[1 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1]; 
end

%------------------------------------------------------------------
% (DLR1): "Plate Experiment" for the active vibration damping of large
%         flexible space structures, example of order 10
%         J. Bals, "Aktive Schwingungsdaempfung flexibler Strukturen",
%         Universitaet Karlsruhe, Fakultaet fuer Elektrotechnik,
%         Germany, 1989
%         reduced system
%------------------------------------------------------------------
if(strmatch('DLR1',ex,'exact'))
 nx=10;nu=2;ny=2;
 A=[0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1;
    -8.4268,0,9.6557,0,-5.083,-0.0253,0,0.0155,0,-0.0112;
    0,-20.2022,0,20.0736,0,0,-0.0244,0,0.0151,0;
    -23.9425,0,-10.7354,0,-147.0685,-0.0049,0,-0.0359,0,-0.0849;
    0,126.1547,0,-132.8028,0,0,0.0947,0,-0.1089,0;
    -39.905,0,6.607,0,-188.4411,-0.0247,0,0.0016,0,-0.1368];
 B=[0,0;0,0;0,0;0,0;0,0;-0.001,-0.001;-0.0133,0.0133;0.048,0.048;
    -0.0516,0.0516;0.0213,0.0213];
 C=[0,0,0,0,0,0.0115,-0.0536,0.9713,-0.2009,-0.5746;
     0,0,0,0,0,0.0115,0.0536,0.9713,0.2009,-0.5746];
 B1=[0,0;0,0;0,0;0,0;0,0;-0.0076,-0.0076;-0.0351,0.0351;0.0972,0.0972;-0.1824,0.1824;0.0748,0.0748];
 C1=[0,1,0,0,0,0,0,0,0,0;-0.8084,0,0.7509,0,-0.9501,0,0,0,0,0];
 D12=[1,0;0,1];
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D11=zeros(nz,nw);
 D21=[0 0;0.0972 0.7509];
end 

%------------------------------------------------------------------
% (DLR2): like (DLR1), original system
%         .mat file generated by
%         /export/home/leibfr/tankred/work/DLR_40_data.m  on msun10
%------------------------------------------------------------------
if(strmatch('DLR2',ex,'exact'))
 n=40; nu=2; ny=2;
 load dlr2_3;  %% ehemals "DLR_Order40.mat"  
end

%------------------------------------------------------------------
% (DLR3): like (DLR2), change in matrix C
%------------------------------------------------------------------
if(strmatch('DLR3',ex,'exact'))
 nx=40; nu=2; ny=2;
 load dlr2_3;  %% ehemals "DLR_Order40.mat"
 C=C_2;  
end

%------------------------------------------------------------------
% (ISS1): International Space Station
%        SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.11
%        W. Draijer, M. Steinbuch, O.H. Bosgra
%        and
%        "Approximation of the International Space Station 1R and 12A
%        flex models",
%        S. Gugercin, A. C. Antoulas and N. Bedrossian , 2001
%------------------------------------------------------------------
if(strmatch('ISS1',ex,'exact'))
  nx=270; nu=3; ny=3;
  load iss1_2;
  A=full(A); B=full(B); C=full(C);
  B1=0.01*B(:,1); C1=sqrt(1e5)*[eye(nx); zeros(3,nx)];
  D12=[zeros(nx,3);eye(nu)]; D21=[0;0;0.05]; D11=zeros(nx+3,1); 
end

%------------------------------------------------------------------
% (ISS2): like (ISS1) with a change in the sensor matrix C in the
%         first row.
%------------------------------------------------------------------
if(strmatch('ISS2',ex,'exact'))
  nx=270; nu=3; ny=3;
  load iss1_2;
  A=full(A); B=full(B); C=full(C);
  B1=0.01*B(:,1); C1=sqrt(1e5)*[eye(nx); zeros(3,nx)];
  D12=[zeros(nx,3);eye(nu)]; D21=[0;0;0.05]; D11=zeros(nx+3,1);
  Ch=C; Ch(1,1:135)=C(1,136:270); Ch(1,136:270)=zeros(1,135);
  C=Ch;
end

%------------------------------------------------------------------
% (CBM): Clamped beam model
%        Los Angeles University Hospital
%        SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.13
%        W. Draijer, M. Steinbuch, O.H. Bosgra
%        and
%        "A survey of model reduction methods for large--scale systems"
%        A. C. Antoulas, D. C. Sorensen and S. Gugercin, 2000   --> Ex. 4.5
%------------------------------------------------------------------
if(strmatch('CBM',ex,'exact'))
 nx=348; nu=1; ny=1; 
 load cbm;
 A=full(A);
 B1=0.09*B; C1=[C; zeros(1,nx)]; D11=zeros(2,1); D12=[0;sqrt(2)];
 D21=0.05; 
end

%------------------------------------------------------------------
% (LAH): Los Angeles University Hospital
%        SLICOT Working note 2002-2 Y. Chahlaoui, P. Van Dooren --> Ex. 2.11
%        W. Draijer, M. Steinbuch, O.H. Bosgra
%        and
%        "A survey of model reduction methods for large--scale systems"
%        A. C. Antoulas, D. C. Sorensen and S. Gugercin, 2000
%------------------------------------------------------------------
if(strmatch('LAH',ex,'exact'))
  nx=48; nu=1; ny=1;
  load lah;
  A=full(A);
  B1=0.01*B; C1=zeros(3,nx); C1(1,1)=1; C1(1,25)=1;
  D12=[0;0;1]; D21=[0.05]; D11=zeros(3,1);
end

%------------------------------------------------------------------
% (ROC1): Four-disk control system 
%         K. Zhou, J. C. Doyle, K. Glover, "Robust and optimal control",
%         Prentice Hall, 1996; p. 517, nc=1
%------------------------------------------------------------------
if(strmatch('ROC1',ex,'exact'))
 A=[-0.161,-6.004,-0.58215,-9.9835,-0.40727,-3.982, 0, 0;
        1,     0,       0,      0,       0,     0, 0, 0;
        0,     1,       0,      0,       0,     0, 0, 0;
        0,     0,       1,      0,       0,     0, 0, 0;
        0,     0,       0,      1,       0,     0, 0, 0;
        0,     0,       0,      0,       1,     0, 0, 0;
        0,     0,       0,      0,       0,     1, 0, 0;
        0,     0,       0,      0,       0,     0, 1, 0];
 B=[1;0;0;0;0;0;0;0];
 C=[0,0,6.4432e-3,2.3196e-3,7.1252e-2,1.0002,0.10455,0.99551];
 q1=1e-6; q2=1;
 B1=sqrt(q2)*[1,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0];
 C1=sqrt(q1)*[0,0,0,0,0.55,11,1.32,18;0,0,0,0,0,0,0,0];
 D12=[0;1]; D21=[0,1];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);
 [nx,nw]=size(B1);[nz,nx]=size(C1);
 D11=zeros(nz,nw);
  
 nc=1; 
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end

%------------------------------------------------------------------
% (ROC2): Transport Aircraft model (Boing flight condition VFC/MFC)
%         D. Gangsaas, K. R. Bruce, J. D. Blight and U.-L. Ly,
%         "Application of Modern Synthesis to Aircraft Control:
%         Three Case Studies", TOAC, Vol.31, Nr.11, pp.995-1014, 1986        
%         Case study III 1), nc=1
%------------------------------------------------------------------
if(strmatch('ROC2',ex,'exact'))
 A=[-0.00702, 0.06339, 0.00518,-0.55566,-0.06112,  0, 0.00712,-0.00566, 0;
    -0.01654,-0.38892, 1.00570, 0.00591,-0.04632,  0, 0.01654, 0.04018, 0;
     0.00061, 0.35210,-0.47381,       0, 1.78620,  0,-0.00061,-0.03638, 0;
           0,       0,       1,       0,       0,  0,       0,       0, 0;
           0,       0,       0,       0,-20.0000, 20,       0,       0, 0;
           0,       0,       0,       0,       0,-30,       0,       0, 0;
           0,       0,       0,       0,       0,  0,-0.55454,       0, 0;
           0,       0,       0,       0,       0,  0,       0,-0.55454, 0.00555;
           0,       0,       0,       0,       0,  0,       0,-0.00555,-0.55454];
 B=[0;0;0;0;0;30;0;0;0];    
 C=[0.005,0.11679,-0.00172, 0,-0.01413, 0,-0.005,-0.01207, 0;
        0,      0,       1, 0,       0, 0,     0,       0, 0];
 C1=(1/sqrt(2))*C(1:1,:);
 B1=[0,0, 0, 0;0,0, 0, 0;0,0, 0, 0;0,0, 0, 0;
     0,0, 0, 0;0,0, 1, 0;1.0531,0, 0, 0;
     0,1.28981, 0, 0;0,-54.51400, 0, 0];
 D12=(1/sqrt(2))*[1];
 D21=[0 0 0 0; 0 0 0 1];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);
 [nx,nw]=size(B1);[nz,nx]=size(C1);   
 D11=zeros(nz,nw);
 
 nc=1; 
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end

%--------------------------------------------------------------------
% (ROC3): Output feedback problem: Wang and Rosenthal   %%ehemals (ROC8)
%         "Output feedback pole placemant with dynamic compensatores" 
%         TOAC, vol.41, Nr. 6, pp. 830-843, 1996;  Example 3.21, nc=2
%--------------------------------------------------------------------
if(strmatch('ROC3',ex,'exact'))
  A=[0,0,0,0,0,0,0,0,-1;1,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,-1;
     0,0,1,0,0,0,0,0,0;0,0,0,1,0,0,0,0,0;0,0,0,0,-1,0,0,0,0;
     0,0,0,0,0,1,0,0,0;0,0,0,0,0,0,1,0,0;0,0,0,0,0,0,0,1,0];
  B=-[0,1;1,0;0,0;0,0;0,1;0,0;0,1;0,0;1,0];
  C=[0,0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0,1];
  
  nc=2;
  [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);    
  A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)]; 
  B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)]; 
  C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
end 

%------------------------------------------------------------------
% (ROC4): Four disk control system
%         D. S. Bernstein and W. M. Haddad, "LQG control with an H_inf
%         performance bound: A Riccati equation approach",
%         TOAC, Vol. 34, Nr. 3, pp. 293-305, 1989;   nc=1 
%------------------------------------------------------------------
if(strmatch('ROC4',ex,'exact'))
  A=[-0.1610 1 0 0 0 0 0 0;
     -6.0040 0 1 0 0 0 0 0;
     -0.5822 0 0 1 0 0 0 0;
     -9.9835 0 0 0 1 0 0 0;
     -0.4073 0 0 0 0 1 0 0;
     -3.9820 0 0 0 0 0 1 0;
      0      0 0 0 0 0 0 1;
      0      0 0 0 0 0 0 0];
 B=[0;0;0.0064;0.00235;0.0713;1.002;0.1045;0.995];
 C=[1 0 0 0 0 0 0 0];
 B1=[B zeros(8,1)];
 C1=[0 0 0 0 0.55 11 1.32 18;
     0 0 0 0  0    0   0   0];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C); 
 [nx,nw]=size(B1);[nz,nx]=size(C1);     
 D11=zeros(nz,nw);
 D12=zeros(nz,nu);
 D21=zeros(ny,nw);
       
 nc=1;
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end

%------------------------------------------------------------------
% (ROC5): Free Gyro-stabilized mirror system 
%         B. M. Chen, "H_inf Control and Its Applications",
%         Springer-Verlag, "Lecture Notes in Control and Information Sciences",
%         Vol.235, 1998; p. 310,  nc=1
%------------------------------------------------------------------
if(strmatch('ROC5',ex,'exact'))
 ab=0.004; bb=0.00128; cb=0.00098; db=0.02;
 eb=0.0049; fb=0.0025; gb=0.00125; lb=0.0032;
 kb=0.0025;
 th3=1; % spin velocity of th flywheel
 psi1=0.011; psi2=0.0011; % some constants
 
 n1=ab+bb+0.5*(eb+gb)+lb; n2=cb+0.25*fb+lb;
 epsA=5e-5;
 A=[0,1,0,0,0,0; 0,0,0,-kb*th3/n1,0,0;0,0,0,1,0,0;
    0,kb*th3/n2,0,0,0,0;0,0,0,0,-epsA,0;0,0,0,0,0,-epsA];
 B=[0,0;1/n1,0;0,0;0,1/n2;0,0;0,0];
 C=[1,0,0,0,0,0;0,0,1,0,0,0;0,0,0,0,1,0;0,0,0,0,0,1];
 B1=[0, 0,0;-1, 0,0;0, 0,0;0,-1,0;0, 0,psi1;0, 0,psi2];
 C1=[1 0 0 0 -1  0;0 0 1 0  0 -1];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C); 
 [nx,nw]=size(B1);[nz,nx]=size(C1);
 D11=zeros(nz,nw);
 D12=zeros(nz,nu);
 D21=zeros(ny,nw); 

 nc=1; 
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end

%------------------------------------------------------------------
% (ROC6): P. Gahinet, "Reliable computation of H_inf central
%         controllers near the optimum", Institut National de Recherche
%         en Informatique et en Automatique, INRIA, Rocquencourt,
%         1992; Example 7.3,   nc=2 
%------------------------------------------------------------------
if(strmatch('ROC6',ex,'exact'))
 A=[1 -1 0;1 1 -1; 0 1 -2];
 B=[1; 0; 1];
 C=[0 -1 1];
 B1=[1 2 0; 0 -1 0; 1 1 0];
 C1=[0 0 0; 1 1 0; -1 0 1];
 D12=[1; 0; 0];
 D21=[0 0 1];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);    
 [nx,nw]=size(B1);[nz,nx]=size(C1);
 D11=zeros(nz,nw); 
 
 nc=2;
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end

%------------------------------------------------------------------
% (ROC7): Flexible actuator 
%         B. Fares, P. Apkarian and D. Noll,
%         "An Augmented Lagrangian Method for a Class of LMI-Constrained
%         Problems in Robust Control Theory",
%         IJOC, Vol. 74, Nr. 4, pp. 348-360;  nc=1
%------------------------------------------------------------------
if(strmatch('ROC7',ex,'exact'))
 A=[0 1 0 0; -1 0 0 0; 0 0 0 1.02; 0.2 0 0 0];
 B=[0;-0.2;0;1];
 C=[1 0 0 0;0 0 1 0];
 C1=[0.1 0 0 0;0 0 0.1 0; 0 0 0 0]; [nz,nx]=size(C1);
 B1=[0;1;0;-0.2]; 
 D12=[0;0;0.2];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);    
 [nx,nw]=size(B1);[nz,nx]=size(C1);
 D11=zeros(nz,nw); 
 D21=zeros(ny,nw);
 
 nc=1;
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)]; 
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)]; 
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)]; 
 C1=[C1 zeros(nz,nc)]; B1=[B1; zeros(nc,nw)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];
end 

%------------------------------------------------------------------
% (ROC8): Augmented three mass spring system      %%ehemals (ROC3)
%         L. El Ghaoui, F. Oustry and M. AitRami,
%         "A cone complementarity linearization algorithm for static
%         output feedback and related problems",
%         TOAC, Vol. 42, Nr. 8, pp. 1171-1176, 1997;   nc=3
%------------------------------------------------------------------
if(strmatch('ROC8',ex,'exact'))
 k=1;
 A=[0 1 0 0 0 0;-k 0 k 0 0 0;0 0 0 1 0 0;
    k 0 -2*k 0 k 0;0 0 0 0 0 1;0 0 k 0 -k 0];
 B=[0;1;0;0;0;0];
 C=[0 0 0 0 1 0];
 B1=[0;0;0;0;0;1];
 C1=eye(6); CC=[0 0 0 0 0 0]; C1=[C1; CC];
 [nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);
 [nx,nw]=size(B1);[nz,nx]=size(C1);
 D11=zeros(nz,nw); 
 D12=[0;0;0;0;0;0;1];
 D21=[1];
 
 nc=3;  
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)];
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)];
 B1=[B1; zeros(nc,nw)]; C1=[C1 zeros(nz,nc)];
 D12=[zeros(nz,nc) D12]; D21=[zeros(nc,nw); D21];  
end

%--------------------------------------------------------------------
% (ROC9): Augmented two mass spring system
%         M. Chilali and P. Gahinet, "H_inf design with pole placement
%         constraints: An LMI approach", TOAC, Vol.41, Nr.3, pp.358-367
%--------------------------------------------------------------------
if(strmatch('ROC9',ex,'exact'))
 k=1; % spring constant 
 A=[0 1 0 0;-k 0 k 0;0 0 0 1;k 0 -k 0];
 B1=[0;0;0;1];
 B=[0;1;0;0];
 C1=eye(4); CC=[0 0 0 0];
 C1=[C1; CC];
 C=[0 0 1 0];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw);
 D12=[0;0;0;0;1];
 D21=[1];
 
 nc=2;
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)]; 
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)]; 
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)]; 
 C1=[C1 zeros(nz,nc)]; B1=[B1; zeros(nc,nw)];
 D12=[zeros(nz,nc) D12];  D21=[zeros(nc,nw); D21]; 
end

%----------------------------------------------------------------------
% (ROC10): Inverted pendulum
%          P. Apkarian and H. D. Tuan, "Robust Conrol via Concave
%          Minimization, Local and Global Algorithms", TOAC, Vol. 45, 
%          Nr. 2, pp. 299-305, 2000; Ex. 1
%----------------------------------------------------------------------
if(strmatch('ROC10',ex,'exact'))
 A=[0, 1, 0, 0, 0; 48.9844, 0, -48.9844, 0, 0; 0, 0, 0, 0.18494, 0; 0, 0, 0, -50.0, 0; 0, 0, -0.5, 0, 0];
 B1=[0,0;0,0;0,0;0,0;0,0.5];
 B=[0;0;50;0;0];
 C1=[0,0,0,0.0036988,0;0,0,0,0,1];
 C=[0,0,1,0,0;1,0,-1,0,0;0,0,0,0,1];
 [nx,nx]=size(A); [nx,nw]=size(B1); [nx,nu]=size(B);
 [nz,nx]=size(C1); [ny,nx]=size(C);
 D11=zeros(nz,nw);
 D12=zeros(nz,nu);
 D21=zeros(ny,nw);
 
 nc=1;
 A=[A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)]; 
 B=[zeros(nx,nc) B; eye(nc) zeros(nc,nu)]; 
 C=[zeros(nc,nx) eye(nc); C zeros(ny,nc)]; 
 C1=[C1 zeros(nz,nc)]; B1=[B1; zeros(nc,nw)];
 D12=[zeros(nz,nc) D12];  D21=[zeros(nc,nw); D21]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,nx]=size(A); [nx,nu]=size(B); [ny,nx]=size(C);

if(isempty(B1))
 B1=eye(nx); C1=eye(nx);
 [nx,nw]=size(B1); [nz,nx]=size(C1);
 D12=[zeros(nz-nu,nu);eye(nu)];  
 D11=zeros(nz,nw); D21=zeros(ny,nw);
else
 [nx,nw]=size(B1); [nz,nx]=size(C1);
end
