close all
clear all
%-------------------%
%% Add MESH2D tool to path
%The MESH2D folder should be added to the path 
%It can be downloaded from https://matlab.mathworks.com/open/fileexchange/v1?id=25555
prompt1 = 'Indicate the full path to the MESH2D folder';
folder = input(prompt1);
addpath(genpath(folder))
run custom_cmap.m
%% Import velocity data for Re=90 or Re=180 VIV case
disp('Importing data')
prompt2 = 'Indicate the full path to the simulation data folder';
dd1 = input(prompt2);
dd=strcat(dd1,'\Results\');
d=dir(fullfile(dd, '*.csv'));
k_dir=natsortfiles({d.name}); %Sort data over time
A = importdata(strcat(dd,'\',char(k_dir(1))));
xUVt=A.data;
n=ceil(size(xUVt,1)); %Spatial Dimension
dt=0.01; %Timestep
start=1; %Starting time (avoid effect of solution initialization)
endt=length(k_dir); %End time
M=ceil((endt-start)/1.5); %Time used for building the ROM
M1=ceil(endt-start); %Overall time

%Allocate the matrices for (x,y) mesh positions, (ux,uy) velocity components
x=zeros(n,M1);
y=zeros(n,M1);
velx=zeros(n,M1);
vely=zeros(n,M1);
for i=1:M1
    %Import solution file for each timestep
    A = importdata(strcat(dd,'\',char(k_dir(i+start))));
    %Convention: ux="u001", uy="u002", x="Points:0",y="Points:1"
    xUVt=A.data;
    x(:,i)=xUVt(:,7);
    y(:,i)=xUVt(:,8);
    velx(:,i)=xUVt(:,2);
    vely(:,i)=xUVt(:,3);
end

%% Load solid displacement data
dOsc=strcat(dd1,'\functionals.txt');
Osc_data=table2array(readtable(dOsc));
dx1_osc=Osc_data(start+1:endt,4); %x displacement
dy1_osc=Osc_data(start+1:endt,5); %y displacement
vx_osc=Osc_data(start+1:endt,6); %ux velocity
vy_osc=Osc_data(start+1:endt,7); %uy velocity
vx_avg=mean(vx_osc(1:M));
vy_avg=mean(vy_osc(1:M));
vosc_avg=[vx_avg; vy_avg];
%% Force Data
mass=0.01319; %mass
k_osc=10; %spring constant
Stil=k_osc/mass;
gtil=(1.2-1.0)/1.2;

%Force calculation (equilibrium position set as (x_e,y_e)=(0,-0.01)
Fx=dx1_osc*k_osc+ddt(vx_osc,dt,41)*mass;
Fy=(dy1_osc-0.01)*k_osc+ddt(vy_osc,dt,41)*mass+gtil*mass*(-0.1);

%Compute time average and center the data
Fx_avg=mean(Fx(1:M)); Fy_avg=mean(Fy(1:M));
Fx=Fx-Fx_avg; Fy=Fy-Fy_avg;
F_avg=[Fx_avg; Fy_avg];
% mesh displacement data
dx_grid=dx1_osc-dx1_osc(1); %x displacement
dy_grid=dy1_osc-dy1_osc(1); %y displacement

%% Construct the new grid
disp('Mesh and ALE map (Laplace) computation')
%We manually locate the position of the solid (t=0)
Rx=0.070;
Ry=0.050;

%We focus on a truncated domain
x_pl=0.5;
x_ph=2.6;
y_pl=0.05;
y_ph=0.95;

%Truncated domain coordinates
xc=1.5+dx1_osc(1)-x_pl;
yc=0.5003+dy1_osc(1)-y_pl;
xmax=x_ph-x_pl;
ymax=y_ph-y_pl;

%Mesh computation (Delauney triangulation) and Laplace solver solution.
[mdisp,points,ind_BC1,ind_BC2,ind_cl,k_cl,mesh]=FE_solver_v3b(xc,yc,Rx,Ry,xmax,ymax);

%% Construct a Euclidean domain (for plotting purposes)
Nx=1500;
Ny=250;
X=linspace(0,xmax,Nx); %Cut a segment in the x direction
Y=linspace(0,ymax,Ny); %Cut a segment in the y direction
[X1,Y1]=meshgrid(X,Y);
X2=linspace(0,max(x(:,1)),Nx); %Cut a segment in the x direction
Y2=linspace(0,max(y(:,1)),Ny); %Cut a segment in the y direction
[X3,Y3]=meshgrid(X2,Y2);
%% Interpolate velocities to the constructed mesh
disp('Flowfield interpolation to MESH2D grid')
n=length(points);
ux_grid=zeros(n,M1);
uy_grid=zeros(n,M1);

for i=1:M1
    %mesh displacement
    ux=dx_grid(i)*mdisp;
    uy=dy_grid(i)*mdisp;
    
    %Moving grid solution
    xi=points(:,1)+ux;
    yi=points(:,2)+uy;
    
    ux_grid(:,i) = griddata(x(:,i)-x_pl,y(:,i)-y_pl,velx(:,i),xi,yi,'linear');
    uy_grid(:,i) = griddata(x(:,i)-x_pl,y(:,i)-y_pl,vely(:,i),xi,yi,'linear');
end

%Ensure no interpolation error was made on the FSI boundary velocities.
ux_grid(ind_BC1,:)=ones(length(ind_BC1),1)*vx_osc.';
uy_grid(ind_BC1,:)=ones(length(ind_BC1),1)*vy_osc.';

%Data centering
ux_mean=mean(ux_grid(:,1:M),2);
uy_mean=mean(uy_grid(:,1:M),2);
ux_grid=ux_grid-ux_mean;
uy_grid=uy_grid-uy_mean;

%% Centering effect (SVD)
[Uxx,Sx1,Vxx]=svd(velx(:,1:M),'econ');
[Uyy,Sy1,Vyy]=svd(vely(:,1:M),'econ');
ux_data_avg=mean(velx(:,1:M),2);
uy_data_avg=mean(vely(:,1:M),2);
[Uxx2,Sx2,Vxx2]=svd(velx(:,1:M)-ux_data_avg,'econ');
[Uyy2,Sy2,Vyy2]=svd(vely(:,1:M)-uy_data_avg,'econ');

%Check the error between the average field and the rank-1 reconstruction by
%the first SVD mode of ux
reproj1 = Uxx(:,1)*Sx1(1,1)*Vxx(:,1)';
cb=find(sqrt(sum((reproj1-repmat(ux_data_avg,1,M)).^2,1))==max(sqrt(sum((reproj1-repmat(ux_data_avg,1,M)).^2,1))));
basis1 = griddata(x(:,1),y(:,1),ux_data_avg-reproj1(:,cb),X3,Y3,'linear');
figure(1)
set(gcf,'position',[10 10 800 280])
[C,h]=contourf(X3,Y3,basis1/max(ux_data_avg)*100,50);
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)
a=colorbar;
xlim([x_pl,x_ph])
ylabel(a,'$\max u_x$ error $( 100 \%)$','Interpreter','latex','FontSize', 20)
colormap(map2)
hColourbar.Label.Position(1) = 3;

set(h,'LineColor','none')
set(h,'LineColor','none')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

%Check the error between the average field and the rank-1 reconstruction by
%the third SVD mode of uy
reproj2 = Uyy(:,3)*Sy1(3,3)*Vyy(:,3)';
cb=find(sqrt(sum((reproj2-repmat(uy_data_avg,1,M)).^2,1))==max(sqrt(sum((reproj2-repmat(uy_data_avg,1,M)).^2,1))));
basis2 = griddata(x(:,1),y(:,1),uy_data_avg-reproj2(:,cb),X3,Y3,'linear');
figure(2)
set(gcf,'position',[10 10 800 280])
[C,h]=contourf(X3,Y3,basis2/max(uy_data_avg)*100,50);
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)
a=colorbar;
xlim([x_pl,x_ph])
ylabel(a,'$\max u_y$ error $( 100 \%)$','Interpreter','latex','FontSize', 20)
colormap(map)
hColourbar.Label.Position(1) = 3;
set(h,'LineColor','none')
set(h,'LineColor','none')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


figure(3)
set(groot,'defaultAxesTickLabelInterpreter','latex');
semilogy(diag(Sx1),'or','MarkerSize',9,'LineWidth',1.5);
hold on
semilogy(diag(Sx2),'*r','MarkerSize',9,'LineWidth',1.5);
semilogy(diag(Sy1),'ob','MarkerSize',9,'LineWidth',1.5);
semilogy(diag(Sy2),'*b','MarkerSize',9,'LineWidth',1.5);

xlim([0 41])
xlabel('$\#$','Interpreter','latex','FontSize', 30)
leg1=legend('${u_x}$','${\tilde{u}_x}$','${u_y}$','${\tilde{u}_y}$');
leg1.FontSize = 22;
set(leg1,'Interpreter','latex');
ylabel('$\sigma$','Interpreter','latex','FontSize', 30);
hold off

%% Construct connectivity matrices
disp('Connectivity matrices computation')
R1=zeros(n,n); %Connectivity matrix used for linear term and H_a quadratic term
Rp=zeros(n,n); %Connectivity matrix used for H_b quadratic term
for j=1:n
    [k1,dum]=find(mesh==j);
    ele1=unique(mesh(k1,:));
    ele1(ele1==j)=[];
    for i=1:length(ele1(:)) %second degree connectivity
        [k2add,dum]=find(mesh==ele1(i));
        k1=[k1; k2add];
    end
    k1=unique(k1);
    ele=unique(mesh(k1,:));
    R1(j,ele(:))=1; %Second-order connectivity
    Rp(j,ele1(:))=1; %First-order connectivity
end
R1=sparse(R1);
Rp=sparse(Rp);

%% Dirichlet BC enforcement on sFOM
%Dirichlet BC on the Fluid-Structure interface
kEL=ind_BC1; %Indices of DOFs on Fluid-Structure interface
R3f=zeros(2*n,2);
R3f(kEL,1)=1; %u_x component
R3f(n+kEL,2)=1; %u_y component
R3f=sparse(R3f);

%Dirichlet BC on the domain inlet
[kin,dum]=find(points(:,1)==0); %Indices of inlet points
in_1=length(kin);
R4f=zeros(2*n,2*in_1);
R4f(kin,1:in_1)=eye(in_1,in_1); %u_x component
R4f(n+kin,in_1+1:2*in_1)=eye(in_1,in_1); %u_y component
R4f=sparse(R4f);

%Eliminate the connectivity of these DOFs from the connectivity matrices
R1(kEL,:)=0;
R1(kin,:)=0;
Rp(kEL,:)=0;
Rp(kin,:)=0;
R1f=[R1 R1; R1 R1];
Rpf=[Rp Rp; Rp Rp];

%% sFOM inference (optimal regularization DOFs)
disp('Flowfield sFOM inference')
Uold=[ux_grid(1:n,1:M-1); uy_grid(1:n,1:M-1)]; %Velocity data matrix (timestep k)
Unew = [ux_grid(1:n,2:M); uy_grid(1:n,2:M)];  %%Velocity data matrix (timestep k+1)
A2=Unew; %LS right-hand side (Ax=b)

%Allocation of sparse inferred matrices
Aexp=zeros(2*n,2*n); %Linear term
Cexp=zeros(2*n,1); %Constant term
Hexp=zeros(2*n,4*n); %Quadratic term ("self" node)
nonz=zeros(n,1);
for i=1:n
    nonz(i)=length(find(Rpf(i,:)));
end
nz_pres=max(nonz);
Hpres=zeros(2*n,(nz_pres+1)*(nz_pres)/2); %Quadratic term ("neighbouring" nodes)
Kexp=zeros(2*n,4*n); %Bilinear term


%Optimal regularization values
R_spA=zeros(2*n,1);
R_spH=zeros(2*n,1);

%L-curve computations
n_reg=20;
b_count1=zeros(2*n,n_reg); %Global L-curve error
x_count1=zeros(2*n,n_reg); %Global L-curve norm
b_count2=zeros(2*n,n_reg); %Local L-curve error
x_count2=zeros(2*n,n_reg); %Local L-curve norm
b_ch=zeros(2*n,1);
x_ch=zeros(2*n,1);
th=0:pi/100:pi/2;
%Regularization values
c=0.01;
c_n=logspace(log10(0.0001),log10(0.5),n_reg);

%Percentage of DOFs for optimal lambda computation
indR=ceil(n/10);
m1=randperm(n,indR);
m=[m1 m1+n];

%sFOM inference with optimal regularization computation
for j=1:2*indR
    i=m(j);
    jp=n*floor(i/(n+1));
    %Connectivity computations
    k=find(R1f(i,:)==1);
    kp=find(Rpf(i,:)==1);
    q=length(k);
    qp=length(kp);
    z=floor(i/(n+1));
    
    if isempty(k)==0
        
        Uself1= repmat(Unew(i-n*z,:),size(k,2),1);
        Uself2= repmat(Unew(i+n-n*z,:),size(k,2),1);
        Umeshx= repmat(vx_osc(2:M).'-vx_avg,size(k,2),1);
        Umeshy= repmat(vy_osc(2:M).'-vy_avg,size(k,2),1);
        %Kronecker product computation
        vkron=[];
        for kk=1:M-1
            longa=[];
            for j=1:qp %Unique kronecker product terms computation
                long1=Unew(kp(j),kk)*Unew(kp(j:end),kk);
                longa=[longa;long1];
            end
            vkron=[vkron, longa];
        end
        qv=size(vkron,1);
        %LS left-hand side (Ax=b)
        A1=[ones(1,M-1); Uold(k,:); Umeshx.*Unew(k,:); Umeshy.*Unew(k,:); vkron; Uself1.*Unew(k,:); Uself2.*Unew(k,:)];
        %Solution via normal equations
        AA=A1*A1.';
        if ismember(i,m)==1
            %L-curve computations
            b_res=zeros(length(c_n),1);
            x_norm=zeros(length(c_n),1);
            for jj=1:length(c_n)
                Reg=eye(5*q+qv+1,5*q+qv+1);
                Reg(1:3*q+1,1:3*q+1)=c_n(jj)*Reg(1:3*q+1,1:3*q+1);
                Reg(3*q+2:end,3*q+2:end)=c_n(jj)*Reg(3*q+2:end,3*q+2:end);
                RR=Reg.'*Reg;
                G1=(AA+RR)\A1*A2(i,:).'; %LS solution
                b_res(jj)=norm(A2(i,:)-G1.'*A1,2); %Solution error
                x_norm(jj)=norm(G1,2); %Solution norm
            end
            %Normalization for solution error and solution norm
            b_resN=(b_res-min(b_res))/(max(b_res)-min(b_res));
            x_normN=(x_norm-min(x_norm))/(max(x_norm)-min(x_norm));
            %Global L-curve
            b_count1(i,:)=b_res(:,1);
            x_count1(i,:)=x_norm(:,1);
            %Optimal regularization computation (L-curve)
            rmin=sqrt(min(x_normN.^2+b_resN.^2));
            i_b=find(x_normN.^2+b_resN.^2==min(x_normN.^2+b_resN.^2));
            b_ch(i)=b_res(i_b);
            x_ch(i)=x_norm(i_b);
            c=c_n(i_b);
            %L-curve plots
            figure(4)
            plot(b_resN(:,1),x_normN(:,1),'DisplayName','L2 regularization')
            xlabel({'$\left\| \hat{\vec{b}} \right\|_{2}$'},'Interpreter','latex')
            ylabel({'$\left\| \hat{\vec{x}} \right\|_{2}$'},'Interpreter','latex')
            hold on
            plot(cos(th)*rmin,sin(th)*rmin,'DisplayName','$ min \left( {\left\| \hat{\vec{b}} \right\|}_{2}^2+{\left\| \hat{\vec{x}}_{2}^2 \right\|}\right)$')
            hold off
            leg1=legend('show');
            set(leg1,'Interpreter','Latex')
            set(gca, 'FontName', 'Times New Roman')
            set(gca, 'FontSize', 16)
            drawnow
        end
        %Solution for optimal regularization parameter
        R_spH(i)=c;
        Reg=eye(5*q+qv+1,5*q+qv+1);
        Reg(1:3*q+1,1:3*q+1)=c*Reg(1:3*q+1,1:3*q+1);
        Reg(3*q+2:end,3*q+2:end)=c*Reg(3*q+2:end,3*q+2:end);
        RR=Reg.'*Reg;
        AA=A1*A1.';
        G1=(AA+RR)\A1*A2(i,:).';
        %Sparse matrices completion
        Cexp(i)=G1(1);
        Aexp(i,k)=G1(2:q+1);
        Kexp(i,[k 2*n+k])=G1(q+2:3*q+1);
        Hpres(i,1:qv)=G1(3*q+2:3*q+qv+1);
        Hexp(i,[k 2*n+k])=G1(3*q+qv+2:end);
    end
    disp(i)
end

%% Interpolation of optimal regularization values
REGu_h = scatteredInterpolant(points(m1,1),points(m1,2),R_spH(m1),'linear','linear');
REGv_h = scatteredInterpolant(points(m1,1),points(m1,2),R_spH(m1+n),'linear','linear');
R_intH=zeros(2*n,1);
R_intH(1:n)= REGu_h(points(:,1),points(:,2));
R_intH(n+1:2*n)= REGv_h(points(:,1),points(:,2));
R_intH(R_intH<min(c_n))=min(c_n);

%% sFOM inference (interpolated regularization DOFs)
nr=1:n;
nr(m1)=[];

for j=1:2*(n-indR)
    z=floor(j/(n-indR+1));
    i=nr(j-z*(n-indR))+z*n;
    k=find(R1f(i,:)==1);
    kp=find(Rpf(i,:)==1);
    q=length(k);
    qp=length(kp);
    z=floor(i/(n+1));
    
    if isempty(k)==0
        Uself1= repmat(Unew(i-n*z,:),size(k,2),1);
        Uself2= repmat(Unew(i+n-n*z,:),size(k,2),1);
        Umeshx= repmat(vx_osc(2:M).'-vx_avg,size(k,2),1);
        Umeshy= repmat(vy_osc(2:M).'-vy_avg,size(k,2),1);
        %Kronecker product computation
        vkron=[];
        for kk=1:M-1
            longa=[];
            for jj=1:qp %Unique kronecker product terms computation
                long1=Unew(kp(jj),kk)*Unew(kp(jj:end),kk);
                longa=[longa;long1];
            end
            vkron=[vkron, longa];
        end
        qv=size(vkron,1);
        %LS left-hand side (Ax=b)
        A1=[ones(1,M-1); Uold(k,:); Umeshx.*Unew(k,:); Umeshy.*Unew(k,:); vkron; Uself1.*Unew(k,:); Uself2.*Unew(k,:)];
        regH=R_intH(i);
        Reg=eye(5*q+qv+1,5*q+qv+1);
        Reg(1:3*q+1,1:3*q+1)=regH*Reg(1:3*q+1,1:3*q+1);
        Reg(3*q+2:end,3*q+2:end)=regH*Reg(3*q+2:end,3*q+2:end);
        RR=Reg.'*Reg;
        %Solution via normal equations
        AA=A1*A1.';
        %LS solution
        G1=(AA+RR)\A1*A2(i,:).';
        b_ch(i)=norm(A2(i,:)-G1.'*A1,2); %Solution error
        x_ch(i)=norm(G1,2); %Solution norm
        %Sparse matrices completion
        Cexp(i)=G1(1);
        Aexp(i,k)=G1(2:q+1);
        Kexp(i,[k 2*n+k])=G1(q+2:3*q+1);
        Hpres(i,1:qv)=G1(3*q+2:3*q+qv+1);
        Hexp(i,[k 2*n+k])=G1(3*q+qv+2:end);
    end
    disp(i)
end

figure(5)
scatter(points(:,1),points(:,2),[],R_intH(1:n),'filled')
hcb = colorbar;
hcb.Title
hcb.Title.String = '\lambda';
colormap jet
xlabel('x (m)')
ylabel('y (m)')
title('u_x')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
drawnow

figure(6)
scatter(points(:,1),points(:,2),[],R_intH(n+1:2*n),'filled')
hcb = colorbar;
hcb.Title
hcb.Title.String = '\lambda';
colormap jet
xlabel('x (m)')
ylabel('y (m)')
title('u_x')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)
drawnow

figure(7)
x_tot1=sum(x_count1.^2,1);
b_tot1=sum(b_count1.^2,1);
plot(x_tot1,b_tot1,'bo')
hold on
x_tot2=sum(x_count2.^2,1);
b_tot2=sum(b_count2.^2,1);
plot(x_tot2,b_tot2,'ro')
plot(sum(x_ch.^2),sum(b_ch.^2),'*','Linewidth',4)
hold off

%% Construct solid-proximity connectivity matrix (for body force sFOM)
disp('Forces sFOM inference')
n_bc=length(ind_cl);
R2=zeros(n_bc,n_bc);
for j=1:n_bc
    i=ind_cl(j);
    [k1,dum]=find(mesh==i);
    k1=unique(k1);
    ele=unique(mesh(k1,:));
    for q=1:length(ele)
        qq=find(ind_cl(:)==ele(q));
        R2(j,qq)=1;
    end
end
R2f=[R2 R2; R2 R2];

%% Force sFOM inference
A2f=[Fx(2:M).' ; Fy(2:M).']; %LS right-hand side (Ax=b)
%Quadratic term computation
quad_f=[];
for i=1:2*n_bc
    k=find(R2f(i,:)==1);
    jk=floor(k/(n_bc+1));
    kr=unique(ind_cl(k-jk*n_bc));
    kr=[kr; kr+n];
    q=length(kr);
    j2=n_bc*floor(i/(n_bc+1));
    Uself1= repmat(Unew(ind_cl(i-j2),:),size(kr,2),1);
    Uself2= repmat(Unew(ind_cl(i-j2)+n,:),size(kr,2),1);
    quad_f=[quad_f; Uself1.*Unew(kr,:); Uself2.*Unew(kr,:)];
end

A1f=[Unew(ind_cl,:); Unew(ind_cl+n,:); quad_f]; %LS left-hand side (Ax=b)
dx_F=2*n_bc+size(quad_f,1);
%L-curve regularization values
cF_n=[0.1 0.15 0.2 0.3];
kF_n=1;
bF_res=zeros(4,1); %Solution error
xF_norm=zeros(4,1); %Solution norm
for j=1:length(cF_n)
    Reg=cF_n(j)*eye(dx_F,dx_F);
    RR=Reg.'*Reg;
    %Solution via normal equations
    AA=A1f*A1f.';
    %LS solution
    G=(AA+RR)\A1f*A2f.';
    bF_res(j)=norm(A2f-G.'*A1f,2);
    xF_norm(j)=norm(G,2);
end
%Optimal regularization computation (L-curve)
bF_resN=(bF_res-min(bF_res))/(max(bF_res)-min(bF_res));
xF_normN=(xF_norm-min(xF_norm))/(max(xF_norm)-min(xF_norm));
i_bF=find(xF_normN.^2+bF_resN.^2==min(xF_normN.^2+bF_resN.^2));
cF=cF_n(i_bF);

%Optimal regularization solution
Reg=cF*eye(dx_F,dx_F);
RR=Reg.'*Reg;
AA=A1f*A1f.';
G=(AA+RR)\A1f*A2f.';
G=G.';
Ds=G(:,1:2*n_bc); %Linear term
Qs=G(:,2*n_bc+1:end); %Quadratic term

%% Flowfield sFOM-POD
disp('sFOM projection (sFOM-POD)')
[U1,S1,V] = svd(Unew,'econ');

semilogy(diag(S1),'*')
r1 = 30;
U1 = U1(:,1:r1);
S1 = S1(1:r1,1:r1);

Atilde=U1'*Aexp*U1; %Linear term
Ctilde=U1'*Cexp; %Constant term
Btilde=U1'*R3f; %Channel inlet BC term
Itilde=U1'*R4f; %Fluid-Structure interface BC term
Ktilde=U1'*Kexp*kron(eye(2,2),U1); %Bilinear term
H_int=zeros(2*n,r1^2); %Quadratic term
for i=1:2*n
    k=find(R1f(i,:)==1);
    kp=find(Rpf(i,:)==1);
    qp=length(kp);
    if isempty(k)==0
        z=floor(i/(n+1));
        U_kron=kron(U1([i-z*n i+n-z*n],:),U1(k,:));
        H_int1=Hexp(i,[k k+2*n])*U_kron; %H_a contribution ("self node")
        
        U_neig=[];
        for j=1:qp
            long1=kron(U1(kp(j),:),U1(kp(j:end),:));
            U_neig=[U_neig;long1];
        end
        H_int2=Hpres(i,1:(1+qp)*qp/2)*U_neig; %H_b contribution ("neighbouring nodes")
        H_int(i,:)=H_int1+H_int2;
    end
end
Htilde=U1'*H_int;

%% Body forces sFOM-POD
Ut_kron=[]; %Quadratic term
for i=1:2*n_bc
    k=find(R2f(i,:)==1);
    jk=floor(k/(n_bc+1));
    kr=unique(ind_cl(k-jk*n_bc));
    kr=[kr; kr+n];
    j2=n_bc*floor(i/(n_bc+1));
    U_kron=kron(U1([ind_cl(i-j2) ind_cl(i-j2)+n],:),U1(kr,:));
    Ut_kron=[Ut_kron; U_kron];
end
Qtilde=Qs*Ut_kron;
Dtilde=Ds*U1([ind_cl;n+ind_cl],:); %Linear term

%% sFOM-POD Simulation
disp('sFOM-POD simulation')

%Initial conditions
u_inlet=[ux_grid(kin,1:M1); uy_grid(kin,1:M1)]; %Inlet velocity
osc_old=[vx_osc(1); vy_osc(1)];
disp_old=[dx1_osc(1); dy1_osc(1)];
u_old=U1'*[ux_grid(1:n,1); uy_grid(1:n,1)];
Fold=[Fx(1);Fy(1)];

endtime=M1; %Simulation end-time
oscillation=zeros(4,endtime);
oscillation(:,1)=[osc_old; disp_old];
errU_rel=zeros(endtime,1);
errV_rel=zeros(endtime,1);

for t=2:endtime
    u_in=u_inlet(:,t);
    %Implicit formulation initialization
    eps=1;
    eps2=1;
    eps3=1;
    disp_test=disp_old;
    osc_test=osc_old;
    u_test=u_old;
    %For each timestep, successive under-relaxation until convergence
    while eps2>10^(-12) || eps>10^(-12) || eps3>10^(-12)
        kk=1;
        %Fluid velocity model
        u_tilde=Atilde*u_old+Btilde*(osc_test-vosc_avg)+Itilde*u_in+Ktilde*kron(osc_test-vosc_avg,u_test)+Htilde*kron(u_test,u_test)+Ctilde;
        %Force model
        F=Dtilde*u_tilde+Qtilde*kron(u_tilde,u_tilde);
        %Solid oscillation
        disp_new=(osc_old+osc_test)*dt/2+disp_old;
        osc_new=osc_old-dt*Stil/2*(disp_new+disp_old-2*[0; 0.01])+dt*gtil*[0; 0.1]+dt/2*(F+Fold+2*F_avg)/mass;
        
        %Successive under-relaxation
        disp_test=0.5*disp_new+0.5*disp_test;
        osc_test=0.5*osc_new+0.5*osc_test;
        u_test=0.5*u_tilde+0.5*u_test;
        eps=norm(osc_test-osc_new,2);
        eps2=norm(u_test-u_tilde,2);
        eps3=norm(disp_test-disp_new,2);
    end
    %Convergence
    disp_old=disp_new;
    osc_old=osc_new;
    u_old=u_tilde;
    Fold=F;
    
    oscillation(:,t)=[osc_test; disp_test]; %solid dynamics
    ux=(disp_old(1)-dx1_osc(1))*mdisp(1:n); %mesh displacement
    uy=(disp_old(2)-dy1_osc(1))*mdisp(1:n);
    u_new=U1*u_tilde; %reconstructed flowfield
    
    if mod(t,floor(endtime/2))==0
        ux1 = griddata(points(1:n,1)+ux,points(1:n,2)+uy,ux_grid(:,t)+ux_mean,X1,Y1,'linear');
        Ux_cDMD = griddata(points(1:n,1)+ux,points(1:n,2)+uy,u_new(1:n)+ux_mean,X1,Y1,'linear');
        uy1 = griddata(points(1:n,1)+ux,points(1:n,2)+uy,uy_grid(:,t)+uy_mean,X1,Y1,'linear');
        Uy_cDMD = griddata(points(1:n,1)+ux,points(1:n,2)+uy,u_new(n+1:2*n)+uy_mean,X1,Y1,'linear');
        
        h1=figure(8);
        h1.Position = [10 10 800 280];
        [C,h]=contourf(X1,Y1,Ux_cDMD,50);
        set(gca, 'FontName', 'Times New Roman')
        set(gcf,'renderer','opengl');
        set(gca, 'FontSize', 18)
        a=colorbar;
        ylabel(a,'u_x (m/s)','FontSize',18,'Rotation',90);
        hColourbar.Label.Position(1) = 3;
        set(h,'LineColor','none')
        set(h,'LineColor','none')
        %c.LineWidth = 1;
        colormap(map2)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(h1,'papersize',[24 10]);
        
        h2=figure(9);
        h2.Position = [10 10 800 280];
        [C,h]=contourf(X1,Y1,ux1,50);
        a=colorbar;
        set(gca, 'FontName', 'Times New Roman')
        set(gca,'FontSize', 18)
        ylabel(a,'u_x (m/s)','FontSize',18,'Rotation',90);
        hColourbar.Label.Position(1) = 3;
        set(h,'LineColor','none')
        set(h,'LineColor','none')
        colormap(map2)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(h2,'papersize',[24 10]);
        drawnow
        
        h3=figure(10);
        h3.Position = [10 10 800 280];
        [C,h]=contourf(X1,Y1,Uy_cDMD,50);
        a=colorbar;
        set(gca, 'FontName', 'Times New Roman')
        set(gca,'FontSize', 18)
        ylabel(a,'u_y (m/s)','FontSize',18,'Rotation',90);
        hColourbar.Label.Position(1) = 3;
        set(h,'LineColor','none')
        set(h,'LineColor','none')
        colormap(map)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(h3,'papersize',[24 10]);
        drawnow
        
        h4=figure(11);
        h4.Position = [10 10 800 280];
        [C,h]=contourf(X1,Y1,uy1,50);
        a=colorbar;
        set(gca, 'FontName', 'Times New Roman')
        set(gca,'FontSize', 18)
        ylabel(a,'u_y (m/s)','FontSize',18,'Rotation',90);
        hColourbar.Label.Position(1) = 3;
        set(h,'LineColor','none')
        set(h,'LineColor','none')
        colormap(map)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(h4,'papersize',[24 10]);
        drawnow
    end
    errU_rel(t)=mean(abs(ux_grid(:,t)-u_new(1:n))./max(ux_grid(:,t)));
    errV_rel(t)=mean(abs(uy_grid(:,t)-u_new(n+1:end))./max(uy_grid(:,t)));
    
end

hf=figure(12);
plot((1:M1)*dt,oscillation(1,:),'r','LineWidth',2)
hold on
plot((1:M1)*dt,vx_osc,'b','LineWidth',2)
plot([M,M]*dt,[5*min(vx_osc) 5*max(vx_osc)],'--')
set(gcf,'position',[10 10 800 280])
legend('ROM r=30','CFD')
xlabel('Time (s)')
ylabel('u_{{osc}_x} (m/s)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)
set(hf,'papersize',[22 8])

hr=figure(13);
plot((1:M1)*dt,oscillation(2,:),'r','LineWidth',2)
hold on
plot((1:M1)*dt,vy_osc,'b','LineWidth',2)
plot([M,M]*dt,[5*min(vy_osc) 5*max(vy_osc)],'--')
set(gcf,'position',[10 10 800 280])
legend('ROM r=30','CFD')
xlabel('Time (s)')
ylabel('u_{{osc}_y} (m/s)')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 18)
set(hr,'papersize',[22 8]);

figure(14)
plot((1:M1)*dt,100*errU_rel,'g','LineWidth',2)
hold on
plot((1:M1)*dt,100*errV_rel,'b','LineWidth',2)
plot([M,M]*dt,100*[0 max(max(errU_rel),max(errV_rel))],'--')
xlabel('Time (s)')
ylabel('Average relative error (100%)')
legend('u','v')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 16)