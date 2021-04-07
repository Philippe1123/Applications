function [ ] = Slope_Slover_JULIA_MATLAB_SF(ProcID)
% ProcID=5;
NumberOfGaussPoints=[7 16 28 37 61];
%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/packages/MultilevelEstimators/KoRDo/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%55


str_interm=str_folder+"Interm/Slope/";


fileC=strcat(str_interm,strcat("Cohesion_",num2str(ProcID)));
fileRES=strcat(str_interm,strcat("Res_",num2str(ProcID)));
fileLevel=strcat(str_interm,strcat("Index_",num2str(ProcID)));
filedoRestriction=strcat(str_interm,strcat("doRestriction_",num2str(ProcID)));
dores=dlmread(filedoRestriction);



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
        str_elements=str_folder+"Mesh/Slope/p_refinement/";
        
        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
        nElem=size(Elements,1);
        EltOpts.GP_k_dp=1;
        %nelem TODO
    end
else
    
    if(HigherOrder==3) %MultiIndex hp refinement
        Order=level(2)+1;
        str_elements=str_folder+"Mesh/Slope/hp_refinement/";
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


% dores=0;
% if(dores==0)
fileC=strcat(str_interm,"Cohesion",level,"_",num2str(ProcID));
file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
% else
%     fileC=strcat(str_interm,strcat(strcat("Cohesion_",strcat(num2str(level+1),"_")),num2str(ProcID)));
%     file_nodes=strcat(str_elements,strcat(strcat("Nodes_L_",num2str(level+1)),".txt"));
%     files_elements=strcat(str_elements,strcat(strcat("Elements_L_",num2str(level+1)),".txt"));

% end


c=dlmread(fileC);

if(HigherOrder==0)
    Nodes=dlmread(file_nodes);
    file_nodes_centers=strcat(str_elements,"ElementsCenter_L",level,".txt");
    nodes_centers=dlmread(file_nodes_centers);
    c=griddata(Nodes(:,2),Nodes(:,3),c,nodes_centers(:,1),nodes_centers(:,2),'cubic');  % Cohesion in N/m^2
else
    %  c=reshape(c,nElem,NumberOfGaussPoints(Order));
    c=transpose(reshape(c,NumberOfGaussPoints(Order),nElem));
end




file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
files_elements=strcat(str_elements,"Elements_L",level,".txt");
Elements=dlmread(files_elements);
Nodes=dlmread(file_nodes);






%
%  figure
%  trisurf(Elements(:,5:7),Nodes(:,2),Nodes(:,3), c_nodes);
% figure
% stem3(centers_coord(:,1),centers_coord(:,2),c_centers,'r')
%%%%%%%%%%

%% Parameters
E=30e6;         % Young's modulus in N/m^2
nu=0.25;        % Poisson ratio
rho=1330;       % Density [kg?m^3]
phi=20;         % Friction angle in degrees

% Match Mohr-Coulomb at shear meridian
alpha_dp=sind(phi)/3;
k_dp=c.*cosd(phi);
% 
% switch Order
%     case 1
%         elementPlane='plane6';
%     case 2
%         elementPlane='plane15';
%     case 3
%         elementPlane='plane28';
%     case 4
%         elementPlane='plane45';
%            case 5
%         elementPlane='plane66';
% end




switch Order
    case 1
        elementPlane='plane3';
    case 2
        elementPlane='plane10';
    case 3
        elementPlane='plane21';
    case 4
        elementPlane='plane36t';
           case 5
        elementPlane='plane55';
end

StepSize=0.1/10;


%% Create model
EltOpts.problem='2dstrain';
EltOpts.nl='elastoplastic';
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);



nElem=size(Elements,1);



% for k =1:nElem, Materials(k,:)={k 'elastoplastic' {'isotropic' [E nu rho] 'dp' [alpha_dp k_dp(k,:)] 'multilinear' [1 k_dp(k,:)+E/300]}};end
for k =1:nElem, Materials(k,:)={k 'elastoplastic' {'isotropic' [E nu rho] 'dp' [alpha_dp k_dp(k,:)] 'multilinear' [1 k_dp(k,:)]}};end

Sections=[1];
%




LeftNodes=selectnode(Nodes,-1e-6,-inf,-inf,1e-6,inf,inf);
RightNodes=selectnode(Nodes,max(Nodes(:,2))-1e-6,-inf,-inf,max(Nodes(:,2))+1e-6,inf,inf);
BottomNodes=selectnode(Nodes,-inf,-1e-6,-inf,inf,1e-6,inf);



% BottomNodes2=model.getNodes('lines',1);
% TopNodes=model.getNodes('lines',3);

% Create list of DOFs
DOF=getdof(Elements,Types);
seldof=[LeftNodes(:,1)+0.01;RightNodes(:,1)+0.01;BottomNodes(:,1)+0.01;BottomNodes(:,1)+0.02];
InitalDof=DOF;
DOF=removedof(DOF,seldof);


% Load vector
DLoads=accel([0 9.81 0],Elements,Types,Sections,{1 E nu rho});
P=elemloads(DLoads,Nodes,Elements,Types,DOF);

%% Compute equilibrium path
% SolverOpts.initalpha=0.2;
% SolverOpts.maxstep=25;
% SolverOpts.maxiter=25;
SolverOpts.initalpha=0.01;
SolverOpts.maxstep=500;
SolverOpts.maxiter=25;
SolverOpts.verbose=false;
SolverOpts.StopConstr=-selectdof(DOF,9+0.02);
%SolverOpts.tolerance=1e-8;

[U,alpha]=solver_arclen(Nodes,Elements,Types,Sections,Materials,DOF,P,SolverOpts);

CornerNode=selectnode(Nodes,15,14,0); % QOI



Options=struct;

%Options.verbose=false;
Div=1;
% while (Div==1)
%     Pload=(0:StepSize:1.0); %load
%     [U,~,~,~,Div]=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,Pload,Options);
%     
%     res=-selectdof(DOF,8+0.02)*[zeros(size(DOF)) U];
%     res=res(end);
%     count=1;
%     if(Div==1||res>1||res<0)
%         StepSize=StepSize/2;
%         Div=1;
%         count=count+1;
%         if(count>10)
%             
%             Div=0;
%             U=zeros(length(Nodes),1);
%         end
%     end
%     
%     
% end
res1=-selectdof(DOF,CornerNode(1)+0.02)*[zeros(size(DOF)) U];
res2=-selectdof(DOF,9+0.02)*[zeros(size(DOF)) U];
res=[res1 res2];
res=res1(end);
%figure
% plot(abs(res1),(0:SolverOpts.maxstep),'--b')

res=alpha(end);

% figure
% plotdisp(Nodes,Elements,Types,DOF,U(:,end))
% 
% U_New=zeros(length(Nodes),1);
% F(160) = struct('cdata',[],'colormap',[]);
% 
% 
% for indx=1:160
%     U_New=zeros(length(Nodes),1);
% 
% for nodeNum=1:length(Nodes(:,1))
% 
%     [n,m]=selectdof(DOF,nodeNum) ;
% 
%     if m==0
% 
%         U_New(nodeNum)=0;
% 
%     elseif(length(m)==1)
% 
%         U_New(nodeNum)=abs(U(m,indx));
% 
% 
%     else
% 
%         for id=1:length(m)
%             U_New(nodeNum)=U(m(id),indx)^2+U_New(nodeNum);
% 
% 
%         end
%         U_New(nodeNum)=sqrt(U_New(nodeNum));
%     end
% 
% end
% 
% trisurf(Elements(:,5:7),Nodes(:,2),Nodes(:,3), U_New);
% view(2)
% F(indx) = getframe(gcf);
% drawnow
% end
dlmwrite(fileRES,res,'delimiter',' ','precision',15);


end

