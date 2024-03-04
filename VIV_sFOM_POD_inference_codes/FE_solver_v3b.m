function[u,points,ind_BC1,ind_BC2,ind_cl,k_cl,tria]=FE_solver_v3b(xc,yc,Rx,Ry,xmax,ymax)

%Includes geometry definition
[points,etri,tria,dum]=Ellipse_Mesh(xc,yc,Rx,Ry,xmax,ymax);

%Find domain boundary elements
ind_b1=find(points(:,1)==min(points(:,1)));
ind_b2=find(points(:,2)==max(points(:,2)));
ind_b3=find(points(:,1)==max(points(:,1)));
ind_b4=find(points(:,2)==min(points(:,2)));

ind_bp=unique([ind_b1; ind_b2; ind_b3; ind_b4]); %boundary points (complete domain)

%Distinguish between ellipse and domain bounds
i_del=ismember(etri(:,1),ind_bp);
ind_b=find(i_del==1); %Indices of boundary edges in etri

%The remain indices correspond to ellipse points
ind_e1=find(i_del==0); %Indices of ellipse edges in etri

edge_e=etri(ind_e1,:);
edge_b=etri(ind_b,:);

%boundary values for displacement
dbx1=ones(size(ind_e1,1),1);
dbx2=zeros(size(ind_b,1),1);

%First ellipse BCs then domain
g1=[dbx1.' dbx2.'];
%Correct indexing of BCs to points

ind_BC1=unique(edge_e);
ind_BC2=unique(edge_b);

ind_BC=[ind_BC1; ind_BC2]; %First ellipse then domain

%We want to find the points "close" to the body
%There we will have both linear and quadratic terms

%<4
ind_cl1=find((points(:,1)-xc).^2/Rx^2+(points(:,2)-yc).^2/Ry^2<2);

k_cl = boundary(points(ind_cl1,1),points(ind_cl1,2));
k_cl=ind_cl1(k_cl);
% figure(2)
% plot(points(k,1),points(k,2),'og','LineWidth',4)
% hold on
% ind_cl2=find(points(:,1)>xc & abs(points(:,2)-yc)<3*Ry);

ind_cl=ind_cl1; %unique([ind_cl1; ind_cl2]);

%--------------------------%
%Since the problem is linear we just need to solve it once!
nodes=tria.';
p=points.';
e=[edge_e.' edge_b.'];


[ijunk,nelem] = size(nodes); %elements
[ijunk,nnode] = size(p); %nodes

gk=zeros(nnode,nnode); %matrix A
gf = zeros(nnode,1); %vector b


for nel = 1:nelem % Begin to assemble by element.
    
    for j=1:3 % The coordinates of the nodes in the
        jj = nodes(j,nel); % element.
        xx(j) = p(1,jj); %x position
        yy(j) = p(2,jj); %y position
    end
    
    for i=1:3
        j = i+1 - fix((i+1)/3)*3; %nearest int->0
        if j == 0
            j = 3;
        end
        m = i+2 - fix((i+2)/3)*3;
        if m == 0
            m = 3;
        end
        %we build each test function
        a(i) = xx(j)*yy(m) - xx(m)*yy(j);
        b(i) = yy(j) - yy(m);
        c(i) = xx(m) - xx(j);
    end
    
    delta = ( c(3)*b(2) - c(2)*b(3) )/2.0; % Area of element
    
    for ir = 1:3
        ii = nodes(ir,nel);
        for ic=1:3
            %matrix elements
            ak = (b(ir)*b(ic) + c(ir)*c(ic))/(4*delta);
            jj = nodes(ic,nel);
            %on i,i elements, there are many terms
            gk(ii,jj) = gk(ii,jj) + ak;
        end
        %For our case f=0
    end
end % End assembling by element.
%------------------------------------------------------
% Now deal with the Dirichlet BC
[ijunk,npres] = size(e); % e(1,:), e(2,:) boundary edge elements
for i=1:npres
    nod=ind_BC(i);
    for k=1:nnode
        gf(k) = gf(k) - gk(k,nod)*g1(i);
        gk(nod,k) = 0;
        gk(k,nod) = 0;
    end
    gk(nod,nod) = 1;
    gf(nod) = g1(i);
end

u=gk\gf; % Solve the linear system.

figure(1)
triplot(tria,points(:,1)+u*0.1,points(:,2)+u*0.1,'b');
hold on
triplot(tria,points(:,1),points(:,2),'k');
hold off
end
