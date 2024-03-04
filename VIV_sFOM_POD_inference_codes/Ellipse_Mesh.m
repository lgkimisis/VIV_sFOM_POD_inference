function[vert,etri,tria,tnum]=Ellipse_Mesh(xc,yc,Rx,Ry,xmax,ymax)
%DEMO0 a very simple example to start with -- mesh a square
%domain with a square hold cut from its centre.

%First we define the edge points
t=-pi:0.07:pi; %0.1
x=xc+Rx*cos(t);
y=yc+Ry*sin(t);
edgep=[x.',y.'];

node_out = [0, 0 ; xmax, 0; xmax, ymax; 0, ymax ];

node=[node_out; edgep];
%------------------------------------------- setup geometry
m=size(node, 1);

edge_out=[1,2; 2,3; 3,4; 4,1];

edge_cyl=[[5:m].' [6:m,5].'];

edge = [edge_out; edge_cyl] ;

%------------------------------------------- call mesh-gen.
   [vert,etri, ...
    tria,tnum] = refine2(node,edge);
%------------------------------------------- call mesh-gen.
    hfun = +.025; %.03           % uniform "target" edge-lengths

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun) ;

end