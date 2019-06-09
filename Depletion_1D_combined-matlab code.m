close all
clear all

%-----------------------------Parameters----------------------------------
%     k1 reaction rate: H2O2+Amplexred---Resofurin;  
%     k2 reaction rate: H2O2+Resofurin--Resazurin;
%   A(H2O2); B(Amplexred); C(Resofurin); D(Resazurin);

L=0.25;     %Length of the box, 
%including excitation area and hydrogel flow [m]
stepx=250;      %Number of step on the x-axis

dx=L/(stepx-1);  %Width of space step(x) [m]
center=round(stepx/2);   %Location of centre
x=0:dx:L;           %Range of x and specifying the grid points

dt=10;               %Width of time step [s] 
stept=1800;        %Number of time steps;  stept*dt= total time;

timer=0;
delay_tstep=1;

%--- A (H2O2)conditions---
input_fractA=20;  %percentage area of input region of A (H2O2)
vis_A=5e-8;       %A (H2O2) Diffusion coefficient/viscocity [m^2/s]

init_conc_A=1e-4;        %A (H2O2),initial concentration [mmol]; 
%init_conc_A /dt=constant rate of H2O2 generation

u_A=zeros(1,stepx);                  %Preallocating u_A
un_A=zeros(1,stepx);                 %Preallocating un_A
stockA = zeros(stept,stepx);

%--- B (Amplexred) conditions---
input_fractB= stepx-input_fractA;    %percentage area of input region of B (Amplexred);
init_conc_B=2.5e-3;             %B(Amplexred)  initial concentration [mmol]

u_B=zeros(1,stepx);             %Preallocating u_B
un_B=zeros(1,stepx);            %Preallocating un_B
stockB = zeros(stept,stepx);

%--- C (Resorufin) conditions---
input_fractC=0;       %percentage area of input region of C (Resorufin) 
init_conc_C=0;            %C (Resorufin) initial concentration [mmol]

u_C=zeros(1,stepx);                  %Preallocating u_C
un_C=zeros(1,stepx);                 %Preallocating un_C
stockC = zeros(stept,stepx);

%--- D (Resazurin) conditions---
input_fractD=0;        %percentage area of input region of D (Resazurin) 
init_conc_D=0;         %D (Resazurin)  initial concentration [mmol]
u_D=zeros(1,stepx);                  %Preallocating u_D
un_D=zeros(1,stepx);                 %Preallocating un_D
stockD = zeros(stept,stepx);

%---Boundary conditions---
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)

% ----------------------------Initial conditions---------------------------
u_A(1:input_fractA) = init_conc_A;
u_initA = u_A;

u_B(input_fractA+1:end) = init_conc_B;
u_initB = u_B;

% ----------------------------1D-Diffusion---------------------------------
i=2:stepx-1;

for k=0:stept
    
    u_tot_A(k+1)=sum(sum(u_A));
    un_A = u_A;
    
    stockA(k+1,:) = u_A;
    stockB(k+1,:) = u_B;
    stockC(k+1,:) = u_C;
    stockD(k+1,:) = u_D;
    
    %---Fick's 2nd law in 1D---
    
    u_A(i)=un_A(i)+(vis_A*dt*(un_A(i+1)-2*un_A(i)+un_A(i-1))/(dx*dx));

    %---Boundary conditions---
    
    u_A(1)=u_A(2)-UnW*dx;
    u_A(stepx)=u_A(stepx-1)+UnE*dx;
    
    %---reaction process
    
    for m=1:stepx-1;
        [u_A(m),u_B(m),u_C(m),u_D(m)]=reaction3(u_A(m),u_B(m),u_C(m),u_D(m),dt);
    end
 
  timer = timer + 1;
    
    if (timer==delay_tstep)
        
        u_A(1:input_fractA) = init_conc_A + u_A(1:input_fractA);
        
        timer = 0;

    end
        
end

% ----------------------------Plots---------------------------------------

