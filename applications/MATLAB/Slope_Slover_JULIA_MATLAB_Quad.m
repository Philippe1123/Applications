function [res1 ] = Slope_Slover_JULIA_MATLAB_Quad(ProcID)

% close all
NumberOfGaussPoints=[3 5 9 17 17]; %Paper
%NumberOfGaussPoints=[1 7 16 28];
%%%%%%%%%SELECT CORRECT SYSTEM
%str_folder="/home/philippe/Desktop/Slope_Matlab/Mesh_Quad/p_ref/Matlab/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
%%%%%%%%%%%%%55


str_interm=str_folder+"Interm10/Slope/";


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
         str_elements=str_folder+"Mesh/Slope/Quad/p_refinement/";

        
        files_elements=strcat(str_elements,"Elements_L",level,".txt");
        Elements=dlmread(files_elements);
        
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


% dores=0;
% if(dores==0)
fileC=strcat(str_interm,"Cohesion",level,"_",num2str(ProcID));
file_nodes=strcat(str_elements,"Nodes_L",level,".txt");
% else
%     fileC=strcat(str_interm,strcat(strcat("Cohesion_",strcat(num2str(level+1),"_")),num2str(ProcID)));
%     file_nodes=strcat(str_elements,strcat(strcat("Nodes_L_",num2str(level+1)),".txt"));
%     files_elements=strcat(str_elements,strcat(strcat("Elements_L_",num2str(level+1)),".txt"));

% end

% fileC="/home/philippe/Desktop/Slope_Matlab/Mesh_Quad/Field"+level;
c=dlmread(fileC);
% c=ones(nElem,NumberOfGaussPoints(1)^2)*8000;
if(HigherOrder==0)
    Nodes=dlmread(file_nodes);
    file_nodes_centers=strcat(str_elements,"ElementsCenter_L",level,".txt");
    nodes_centers=dlmread(file_nodes_centers);
  %  c=c(1:nElem);
    c=griddata(Nodes(:,2),Nodes(:,3),c,nodes_centers(:,1),nodes_centers(:,2),'cubic');  % Cohesion in N/m^2
% F=scatteredInterpolant(Nodes(:,2),Nodes(:,3),c);
% c=F(nodes_centers(:,1),nodes_centers(:,2));

else
    
    %c=c(1:NumberOfGaussPoints(Order)*nElem);
    %  c=reshape(c,nElem,NumberOfGaussPoints(Order));
    c=transpose(reshape(c,NumberOfGaussPoints(level_num+1)^2,nElem));
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

switch Order
    case 1
        elementPlane='plane9_h';
        EltOpts.nXi=NumberOfGaussPoints(1);


    case 2
        elementPlane='plane16_h';
              EltOpts.nXi=NumberOfGaussPoints(2);


    case 3
        elementPlane='plane25_h';
        EltOpts.nXi=NumberOfGaussPoints(3);

        
    case 4
        elementPlane='plane36_h';
        EltOpts.nXi=NumberOfGaussPoints(4);

        
    case 5
        elementPlane='plane49_h';
        EltOpts.nXi=NumberOfGaussPoints(5);

        
    case 7
        elementPlane='plane64';
        EltOpts.nXi=NumberOfGaussPoints(1);

        
    case 8
        elementPlane='plane81';
        EltOpts.nXi=NumberOfGaussPoints(1);

        
    case 9
        elementPlane='plane55';
                EltOpts.nXi=NumberOfGaussPoints(1);

end



StepSize=0.1/10;


%% Create model
EltOpts.problem='2dstrain';
EltOpts.maxiter=50;
EltOpts.nl='elastoplastic';
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);


% figure
% plotelem(Nodes,Elements,Types,'GCS','off','LineWidth',3);
nElem=size(Elements,1);



 for k =1:nElem, Materials(k,:)={k 'elastoplastic' {'isotropic' [E nu rho] 'dp' [alpha_dp k_dp(k,:)] 'multilinear' [1 k_dp(k,:)+E/300]}};end

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




CornerNode=selectnode(Nodes,15,14,0); % QOI



Options=struct;

Options.verbose=false;
Options.maxiter=200;
StepSize=0.1/10;

Div=0;
% while (Div==1)
    Pload=(0:StepSize:1); %load
     [U,~,~,~,Div]=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,Pload,Options,[],[],[],[],[]);
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
%res2=-selectdof(DOF,9+0.02)*[zeros(size(DOF)) U];
%res=[res1 res2];
res=res1(end);

% figure
% plotdisp(Nodes,Elements,Types,DOF,U(:,end),'DispScal',20)


% DisplacedNodes=Nodes;
% 
% for i=1:length(Nodes)
%     
% xdis=selectdof(DOF,Nodes(i,1)+0.01)*[zeros(size(DOF)) U(:,end)];    
% xdis=xdis(end);
% 
% ydis=selectdof(DOF,Nodes(i,1)+0.02)*[zeros(size(DOF)) U(:,end)];    
% ydis=ydis(end);
% 
% DisplacedNodes(i,2)=DisplacedNodes(i,2)+xdis*20;
% DisplacedNodes(i,3)=DisplacedNodes(i,3)+ydis*20;
% 
% 
% end
% hold on
% plotnodes(DisplacedNodes,'numbering','off','*r','GCS','off','LineWidth',3)
% 
dlmwrite(fileRES,res,'delimiter',' ','precision',15);


end

