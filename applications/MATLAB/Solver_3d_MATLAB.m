function [] = Solver_3d_MATLAB(ProcID)
% ProcID=1
%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%
% ProcID=1
% str_folder="/home/philippe/Desktop/3d_Matlab/";
% str_interm=str_folder+"Interm/3d/";

str_interm="/scratch/leuven/330/vsc33032/Interm/3d/";
str_B="/scratch/leuven/330/vsc33032/Matrices/3d/";

fileLevel=strcat(str_interm,strcat("Index_",num2str(ProcID)));
fileRES=strcat(str_interm,strcat("Res_",num2str(ProcID)));

level=dlmread(fileLevel);
level_num=level;
HigherOrderHandle=strcat(str_interm,strcat("HighOrder_",num2str(ProcID)));
HigherOrder=dlmread(HigherOrderHandle);




if(length(level)==1)
    
    str=strcat('_',num2str(level(1)));
    level=str;
    
    
    if(HigherOrder==0) %Pure h refinement
        Order=1;
        
        str_elements=str_folder+"Mesh/Slope/h_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
    elseif(HigherOrder==1) %Pure p refinement
        Order=level_num+1;
         str_elements=str_folder+"Mesh/3d_beam/p_refinement/";
%          str_elements=str_folder+"GenerateMeshes_1_p_dev_beam/p_ref/Matlab/p_refinement/";

        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        files_Integration_points=strcat(str_elements,"GaussPoints_L",level,".txt");
        InterpPoints=dlmread(files_Integration_points);
        InterpPoints=InterpPoints(:,2:4);
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        %nelem TODO
    end
else
    
    if(HigherOrder==3) %MultiIndex hp refinement
        Order=level(2)+1;
        str_elements=str_folder+"Mesh/hp_refinement/Slope/";
        str='_';
        for len=1:length(level)
            str=strcat(str,num2str(level(len)));
            
        end
        
        files_elements=strcat(str_elements,"Elements_L",str,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        level=str;
        
    end
end













fileC=strcat(str_interm,"Cohesion",level,"_",num2str(ProcID));
file_nodes=strcat(str_elements,"Nodes_L",level,".txt");



c=dlmread(fileC);

if(HigherOrder==0)
    Nodes=dlmread(file_nodes);
    file_nodes_centers=strcat(str_elements,"ElementsCenter_L",level,".txt");
    nodes_centers=dlmread(file_nodes_centers);
    c=griddata(Nodes(:,2),Nodes(:,3),c,nodes_centers(:,1),nodes_centers(:,2),'cubic');  % Cohesion in N/m^2


else
    x=0:0.025:5;
    y=0:0.005:1;
    z=0:0.0025:0.5;
    [X,Y,Z]=meshgrid(x,y,z);

 x=reshape(X,size(X,1)*size(X,2)*size(X,3),1);
 y=reshape(Y,size(Y,1)*size(Y,2)*size(Y,3),1);
 z=reshape(Z,size(Z,1)*size(Z,2)*size(Z,3),1);

Values=reshape(c,size(X,1),size(X,2),size(X,3));



vq=interp3(X,Y,Z,Values,InterpPoints(:,1),InterpPoints(:,2),InterpPoints(:,3));

reshapeValue=size(vq,1)/nElem;
    
uncertainVal=reshape(vq,reshapeValue,nElem);
uncertainVal=transpose(uncertainVal)*1000;
    
end




file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
files_elements=strcat(str_elements,"Elements_L",level,".txt");
Elements=dlmread(files_elements);
Nodes=dlmread(file_nodes);











switch level_num
    
    case 0
    increment=1;
    elementPlane='solid8';
    EltOpts.nXi=2;
    
    case 1     
    increment=2;
    elementPlane='solid27_h';
     EltOpts.nXi=3;
     
    case 2
    increment=3;
    elementPlane='solid64_h';
    EltOpts.nXi=4;
    
    case 3
    increment=4;
    elementPlane='solid125_h';
    EltOpts.nXi=5;
    
    case 4
    increment=5;
    elementPlane='solid216_h';
    EltOpts.nXi=6;
    
end









%%%%%%%%%%

%% Parameters
Ek=ones(nElem,1).*30e6;         % Young's modulus in N/m^2
nu=0.15;        % Poisson ratio
el1=Elements(1,5:end);
nodesPos=Nodes(el1,:);
nod=nodesPos(:,2:end);
%     figure
%     for i=1:length(nod)
%         hold on
%         
%         plot3(nod(i,1),nod(i,2),nod(i,3),'*b')
%         text(nod(i,1),nod(i,2),nod(i,3),num2str((i)))
%         hold on        
%         
%         i=i+1;
%         
%     end



StepSize=0.1/10;


%% Create model
EltOpts.problem='3d';
EltOpts.nl='linear';
EltOpts.Coeffs=dlmread('coeffs_se125');
EltOpts.bendingmodes=false;
StrintPt='CC_SparseGrid_3d';
% StrintPt='gausslegendre_3d';

intPt=str2func( StrintPt );


EltOpts.IntegrationFunc=intPt;
EltOpts.detJac=dlmread(str_B+"detJac");
EltOpts.B=dlmread(str_B+"B_"+num2str(level_num+1));
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);







 for k = 1:nElem, Materials(k,:) = {k 'linear' {'isotropic' {uncertainVal(k,:) nu} }}; end

%  Ek=ones(nElem,1).*30e6;         % Young's modulus in N/m^2
%  for k = 1:nElem, Materials(k,:) = {k 'linear' {'isotropic' [Ek(k,:) nu] }}; end

Sections=[1];
%




LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
BottomNodes=selectnode(Nodes,-inf,-1e-6,-inf,inf,1e-6,inf);



% BottomNodes2=model.getNodes('lines',1);
% TopNodes=model.getNodes('lines',3);

LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
DOF = getdof(Elements,Types);
sdof = [0.04;0.05;0.06;LeftNodes(:,1)+0.01;RightNodes(:,1)+0.01;LeftNodes(:,1)+0.02;RightNodes(:,1)+0.02;LeftNodes(:,1)+0.03;RightNodes(:,1)+0.03]; 
DOF = removedof(DOF,sdof);
sdog_algo=unique(floor(sdof));





%%%%%%%%%%%%%%%%%%%%%
% forceNodes=selectnode(Nodes,min(Nodes(:,2))-1e-6,min(Nodes(:,3))-1e-6,max(Nodes(:,4))-1e-6,max(Nodes(:,2))+1e-6,max(Nodes(:,3))+1e-6,max(Nodes(:,4))+1e-6);
% forceN=forceNodes(:,2:end);
% [srted,posid]=sortrows(forceN,[1 2]);
% n_nodes_x=length(uniquetol(Nodes(:,2),10^-7));
% n_nodes_y=length(uniquetol(Nodes(:,3),10^-7));
% n_nodes_z=length(uniquetol(Nodes(:,4),10^-7));
% forceN=forceNodes(posid,:);
% forceN=forceN(:,1);
% forceNodes=forceN+0.02;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
forceNodes=selectnode(Nodes,min(Nodes(:,2))-1e-6,max(Nodes(:,3))-1e-6,max(Nodes(:,4))./2-1e-6,max(Nodes(:,2))+1e-6,max(Nodes(:,3))+1e-6,max(Nodes(:,4))./2+1e-6);
forceN=forceNodes(:,2:end);
[srted,posid]=sortrows(forceN,[1 2]);
forceN=forceNodes(posid,:);
forceN=forceN(1:increment:end,:);
forceN=forceN(:,1);
forceNodes=forceN+0.02;
%%%%%%%%%%%%%%%%%

% forceNodes=selectnode(Nodes,max(Nodes(:,2))/2-1e-6,max(Nodes(:,3))./2-1e-6,max(Nodes(:,4))-1e-6,max(Nodes(:,2))/2+1e-6,max(Nodes(:,3))./2+1e-6,max(Nodes(:,4))+1e-6);
% forceN=forceNodes(:,1:end);
% forceN=forceN(:,1);
% forceNodes=forceN+0.02;



%% Load




%%%%%%%%%
%     Force=10000000*(0:1:1);
%     PLoad=Force;
%     
%  ar_x=[0:2*pi/(n_nodes_x-1):2*pi];
%  ar_y=[0:2*pi/(n_nodes_y-1):2*pi];
% vec_x=(-cos(ar_x)+1)/sum((-cos(ar_x)+1));
% vec_y=(-cos(ar_y)+1)/sum((-cos(ar_y)+1));
% mat=kron(vec_x,vec_y');
% vec=reshape(mat,[1,size(mat,2)*size(mat,1)]);
% P = nodalvalues(DOF,forceNodes,-vec);
%%%%%%%%%%%%%
    Force=10000*(0:1:1);
    PLoad=Force;
ar=[0:2*pi/(length(forceNodes)-1):2*pi];
vec=(-cos(ar)+1)/sum((-cos(ar)+1));
P = -nodalvalues(DOF,forceNodes',vec);
%%%%%%%%%%
%     Force=10000000*(0:1:1);
%     PLoad=Force;
% P = -nodalvalues(DOF,forceNodes',-1);
%  figure;
%  plotnodes(Nodes);
% 
%  figure;
%  plotelem(Nodes,Elements,Types);

%% Compute equilibrium path
Options=struct;
Options.verbose=false;


% filename=string(elementPlane)+"_"+string(StrintPt)+"_"+string(EltOpts.nXi)+".mat";
% filename=char(filename);
% if exist(filename, 'file') == 2
% ElemCache=load(filename,'ElemCache');
% ElemCache=ElemCache.ElemCache;
% else
% ElemCache={};
% end



[U,HistPar,critPnt,ElemCache,Divergence]=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,PLoad,Options,[],[],[],[],[]);

% if exist(filename, 'file') == 2
% 
% else
% save(filename,'ElemCache')
% end
% 
% [K,L]=size(U);




% QoINodes=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6,inf,min(Nodes(:,3))+1e-6,1e-6);
% QoINodes_2=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6+0.125,inf,min(Nodes(:,3))+1e-6,1e-6+0.125);
QoINodes_3=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6+0.25,inf,min(Nodes(:,3))+1e-6,1e-6+0.25);
% QoINodes_4=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6+0.375,inf,min(Nodes(:,3))+1e-6,1e-6+0.375);
% QoINodes_5=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-1e-6+0.5,inf,min(Nodes(:,3))+1e-6,1e-6+0.5);

u=U(:,end);
% u_out=selectdof(DOF,QoINodes(:,1)+0.02)*U(:,end);
% u_out_2=selectdof(DOF,QoINodes_2(:,1)+0.02)*U(:,end);
u_out_3=selectdof(DOF,QoINodes_3(:,1)+0.02)*U(:,end);
% u_out_4=selectdof(DOF,QoINodes_4(:,1)+0.02)*U(:,end);
% u_out_5=selectdof(DOF,QoINodes_5(:,1)+0.02)*U(:,end);


[s,id]=sort(QoINodes_3(:,2));
QoINodes_3=QoINodes_3(id,2);
u_out_3=u_out_3(id);
% figure
% plot(QoINodes_3,u_out_3,'-*')
 
%  figure
%  plot(QoINodes_3(1:increment:end),u_out_3(1:increment:end),'-*')
dlmwrite(fileRES,min(u_out_3),'delimiter',' ','precision',15);


% xlabel('x axis')
% ylabel('y axis')
% zlabel('z axis')
end

