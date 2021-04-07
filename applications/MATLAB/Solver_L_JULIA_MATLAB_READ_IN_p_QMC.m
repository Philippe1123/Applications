
function [] = Solver_L_JULIA_MATLAB_READ_IN(ProcID)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%ProcID=2
%NumberOfGaussPoints=[9 25 49 81]; % paper 
NumberOfGaussPoints=[25    81   121   225   289   441   625   729];
str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data/";


str_interm=str_folder+"Interm2/Beam/";

fileQoI=strcat(str_interm,strcat("NQoI"));
fileMaxLevelI=strcat(str_interm,strcat("MaxLVL"));

MaxLVL=dlmread(fileMaxLevelI);
NQoI=dlmread(fileQoI);

fileRES=strcat(str_interm,strcat("Res_",num2str(ProcID)));
fileLevel=strcat(str_interm,strcat("Index_",num2str(ProcID)));
level=dlmread(fileLevel);
level_num=level;
HigherOrderHandle=strcat(str_interm,strcat("HighOrder_",num2str(ProcID)));
HigherOrder=dlmread(HigherOrderHandle);
%filedoRestriction=strcat(str_interm,strcat("doRestriction_",num2str(ProcID)));

t=250;
nu=0.15;





%dores=dlmread(filedoRestriction);




if(length(level)==1)
    
    str=strcat('_',num2str(level(1)));
    level=str;
    
    
    if(HigherOrder==0) %Pure h refinement
        Order=1;
        
        str_elements=str_folder+"Mesh/Beam/h_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
    elseif(HigherOrder==1) %Pure p refinement
        EltOpts.GP='1';

        Order=level_num+1;
        str_elements=str_folder+"Mesh/Beam/p_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        %nelem TODO
    end
else
    
    if(HigherOrder==3) %MultiIndex hp refinement
        Order=level(2)+1;
        str_elements=str_folder+"Mesh/Beam/hp_refinement/";
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

fileE=strcat(str_interm,"E",level,"_",num2str(ProcID));
E=dlmread(fileE);
file_nodes=strcat(str_elements,"Nodes_L",level,".txt");



switch Order
    case 1
        elementPlane='plane4';
    case 2
        elementPlane='plane9';
    case 3
        elementPlane='plane16';
    case 4
        elementPlane='plane25';
    case 5
        elementPlane='plane36';
    case 6
        elementPlane='plane49';
    case 7
        elementPlane='plane64';
    case 8
        elementPlane='plane81';
    case 9
        elementPlane='plane100';
end

EltOpts.nl='linear';
EltOpts.bendingmodes=false;
Types = {1 elementPlane EltOpts};   % {EltTypID EltName EltOpts} 
Sections = [1 t];               % [SecID t] 

    Options=struct;
    Options.verbose=false;

%% Mesh
% Line1 = [0 0 0; Lx 0 0];
% Line2 = [Lx 0 0; Lx Ly 0];
% Line3 = [Lx Ly 0; 0 Ly 0];
% Line4 = [0 Ly 0; 0 0 0];
% 
% [Nodes,Elements,Edge1,Edge2,Edge3,Edge4,Normals] = makemesh(Line1,Line2,Line3,Line4,m,n,Types(1,:),1,1);
% Nodes(:,4) = 0;
if(HigherOrder==0)
    Nodes=dlmread(file_nodes);
    file_nodes_centers=strcat(str_elements,"ElementsCenter_L",level,".txt");
    nodes_centers=dlmread(file_nodes_centers);
   % Ek=griddata(Nodes(:,2),Nodes(:,3),E,nodes_centers(:,1),nodes_centers(:,2),'cubic');  % Cohesion in N/m^2
   F=scatteredInterpolant(Nodes(:,2),Nodes(:,3),E);
   Ek=F(nodes_centers(:,1),nodes_centers(:,2));
%   Ek=E;
else
    %  c=reshape(c,nElem,NumberOfGaussPoints(Order));
    Ek=transpose(reshape(E,NumberOfGaussPoints(Order),nElem));
   % Ek=ones(nElem,NumberOfGaussPoints(Order))*3e4;
end

file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
files_elements=strcat(str_elements,"Elements_L",level,".txt");
Elements=dlmread(files_elements);
Nodes=dlmread(file_nodes);


filesNodesMax=strcat(str_elements,"Nodes_L_",num2str(MaxLVL),".txt");
NodesMax=dlmread(filesNodesMax);
% 
% nElem=size(Elements,1);
% Elements(:,4) = 1:nElem;
% 
% 
% if(size(E,2)==1)
% Ek = E*ones(nElem,1);
% else
% %[j,k]=size(E);
% %size(E)
% %nElem
% Ek = reshape(E,[nElem,1]);
% end


Materials = cell(nElem,3);
% for k = 1:nElem, Materials(k,:) = {k 'elastoplastic' {'isotropic' [Ek(k) nu] 'vm' [fy] 'multilinear' [1 fy+E_in/20]}}; end
for k = 1:nElem, Materials(k,:) = {k 'linear' {'isotropic' [Ek(k,:) nu] }}; end

%  figure;
%  plotnodes(Nodes);
% 
%  figure;
%  plotelem(Nodes,Elements,Types);
% pause(10)
%% DOFs
LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
DOF = getdof(Elements,Types);
sdof = [0.03;0.04;0.05;0.06;LeftNodes(:,1)+0.01;RightNodes(:,1)+0.01;LeftNodes(:,1)+0.02;RightNodes(:,1)+0.02]; 
DOF = removedof(DOF,sdof);
sdog_algo=unique(floor(sdof));



%forceNodes=selectnode(Nodes,max(Nodes(:,2))/2-1e-6,max(Nodes(:,3))-1e-6,-inf,max(Nodes(:,2))/2+1e-6,max(Nodes(:,3))+1e-6,inf);
forceNodes=selectnode(Nodes,max(Nodes(:,2))/2-1e-6,-inf,-inf,max(Nodes(:,2))/2+1e-6,inf,inf);

forceNodes=forceNodes(:,1);

[posy,srt]=sort(Nodes(forceNodes,3),'descend');
forceNodes=forceNodes(srt);
% forceNodes=forceNodes(1:(end-1));
forceNodes=forceNodes(1);

forceNodes=forceNodes+0.02;
seldof=forceNodes';
%% Load

    Force=10000000*(0:1:1);
    PLoad=Force;
% P = nodalvalues(DOF,seldof,[-1*0.5 -1*ones(1,size(seldof,2)-2) -1*0.5]./(size(seldof,2)-1));
P = nodalvalues(DOF,seldof,[-1*ones(1,size(seldof,2))]);

%% Compute displacements
U=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,PLoad,Options);
[K,L]=size(U);

%fprintf('before loop')

% % % % for ik=1:L
% % % % 
% % % % u=U(:,ik);
% % % % u=u';
% % % %         matxDir=zeros(nely+1,nelx+1);
% % % %         matyDir=zeros(nely+1,nelx+1);
% % % %         x=1;  y=1;  z=2;  t=1;v=1;
% % % %         while(x<=(nelx+1))
% % % %             y=1;
% % % %             while(y<=(nely+1))
% % % %                 if(sum(ismember(v,sdog_algo))==0)
% % % %                     matxDir(y,x)=u(t);
% % % %                     matyDir(y,x)=u(z);
% % % %                     z=z+2;  t=t+2;
% % % %                 else
% % % %                     matxDir(y,x)=0;
% % % %                     matyDir(y,x)=0;
% % % %                 end
% % % %                 v=v+1;
% % % %                 y=y+1;
% % % %             end
% % % %             x=x+1;
% % % %         end
% % % %         Matx=matxDir;
% % % %         Maty=matyDir;
% % % %         
% % % %         [n,m]=size(matxDir);
% % % %     maxNdof=2*(nelx+1)*(nely+1);    
% % % %             Pc=zeros(maxNdof,1);
% % % %             x=1;
% % % %             i=1;
% % % %             y=1;
% % % %             while(y<=m)
% % % %                 x=1;
% % % %                 while(x<=n)
% % % %                     
% % % %                     
% % % %                     Pc(i)=matxDir(x,y);
% % % %                     Pc(i+1)=matyDir(x,y);
% % % %                     i=i+2;
% % % %                     x=x+1;
% % % %                 end
% % % %                 y=y+1;
% % % %             end 
% % % %             
% % % %             if(ik==1)
% % % %             u_new=Pc;
% % % %             else
% % % %             u_new=[u_new Pc];    
% % % %             
% % % %             end
% % % % end

%u=U(:,end);
%u_out=min(u);
U=abs(U);
if(NQoI==1)

res_node=selectnode(Nodes,max(Nodes(:,2))/2-1e-6,min(Nodes(:,3))-1e-6,-inf,max(Nodes(:,2))/2+1e-6,min(Nodes(:,3))+1e-6,inf);
u_out=selectdof(DOF,res_node(:,1)+0.02)*U(:,end);

dlmwrite(fileRES,u_out,'delimiter',' ','precision',15);


else
res_node=selectnode(Nodes,-inf,min(Nodes(:,3))-1e-6,-inf,inf,min(Nodes(:,3))+1e-6,inf);

[Nodes_x,I]=sort(res_node(:,2));

u_out=selectdof(DOF,res_node(I,1)+0.02)*U(:,end); 



res_node_max=selectnode(NodesMax,-inf,min(NodesMax(:,3))-1e-6,-inf,inf,min(NodesMax(:,3))+1e-6,inf);

Nodes_x_max=sort(res_node_max(:,2));

Nodes_y_max=res_node_max(:,3);

u_out=interp1(Nodes_x,u_out,Nodes_x_max);

    
dlmwrite(fileRES,u_out,'delimiter',' ','precision',15);
    
end


% [Matx,Maty] = Extract_matrix_Displacements_From_Vector(U,nelx,nely);
%% Plot displacements
% figure
% plotdisp(Nodes,Elements,Types,DOF,U(:,end))
% 
% figure;
% animdisp(Nodes,Elements,Types,DOF,U)
% 
% figure;
% seldof = find(Nodes(:,2)==Lx/2 & Nodes(:,3)==0 & round(Nodes(:,4)*1000)/1000==0)+0.02;
% plot(selectdof(DOF,seldof)*U,PLoad)
% xlabel('Displacement [mm]')
% ylabel('Force [N]')
% hold on


% seldof = find(Nodes(:,2)==Lx/2 & Nodes(:,3)==Ly & round(Nodes(:,4)*1000)/1000==0)+0.02;
% plot(selectdof(DOF,seldof)*U,PLoad)
% xlabel('Displacement [mm]')
% ylabel('Force [N]')
% legend('Lower level','Upper level')

end


