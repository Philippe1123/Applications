function [ ] = Slope_Slover_JULIA_MATLAB_MDPI_2(ProcID)
%ProcID=2

NumberOfGaussPoints=[7 13  19 25 28 33 37 61 73]; %Paper
%NumberOfGaussPoints=[1 7 16 28];
%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%55


str_interm=str_folder+"Interm7/Slope/";


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
        
        str_elements=str_folder+"Mesh/Slope/MLMC_3/h_refinement/";
        
        
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


c=dlmread(fileC);

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

switch Order
    case 1
        elementPlane='plane3';
        EltOpts.nXi=NumberOfGaussPoints(1);
        EltOpts.Sort_nXi=[5,1,3,4,2,7,6];
    case 2
        elementPlane='plane6';
         EltOpts.nXi=NumberOfGaussPoints(2);
         EltOpts.Sort_nXi=[6,3,13,1,4,7,10,8,11,9,2,12,5];

    case 3
        elementPlane='plane10';
        EltOpts.nXi=NumberOfGaussPoints(3);
        EltOpts.Sort_nXi=[10,1,16,6,5,9,14,17,19,18,8,15,11,13,2,3,12,7,4];

    case 4
        elementPlane='plane15';
        EltOpts.nXi=NumberOfGaussPoints(4);
        EltOpts.Sort_nXi=[7,10,4,19,11,9,22,3,13,12,20,21,25,2,5,18,1,16,17,8,14,6,15,24,23];
        
    case 5
        elementPlane='plane21';
        EltOpts.nXi=NumberOfGaussPoints(5);
        EltOpts.Sort_nXi=[26,13,19,15,12,27,17,20,22,21,11,18,2,4,8,9,3,16,10,28,1,25,6,23,5,14,7,24];
        
    case 6
        elementPlane='plane28';
        EltOpts.nXi=NumberOfGaussPoints(6);
        EltOpts.Sort_nXi=[28,8,27,6,9,29,24,22,25,23,7,26,13,14,1,3,15,5,2,31,19,33,12,30,32,4,11,10,16,21,17,18,20];
        
    case 7
        elementPlane='plane36t';
        EltOpts.nXi=NumberOfGaussPoints(7);
        EltOpts.Sort_nXi=[21,10,26,18,9,25,28,33,36,31,8,29,2,4,14,15,3,19,16,24,30,20,12,22,23,17,13,11,27,32,37,34,35,5,6,1,7];
        
    case 8
        elementPlane='plane45';
        EltOpts.nXi=NumberOfGaussPoints(8);
        EltOpts.Sort_nXi=[4,29,14,40,30,26,46,13,25,24,44,1,61,9,10,43,8,41,42,27,31,5,33,57,54,19,34,56,20,21,32,55,51,28,38,35,39,7,18,16,2,12,22,3,11,37,6,23,15,17,52,50,47,36,48,53,49,45,58,59,60];
        
    case 9
        elementPlane='plane55';
        EltOpts.nXi=NumberOfGaussPoints(9);
        EltOpts.Sort_nXi=[1,37,21,50,36,29,64,20,35,34,51,63,73,8,11,52,7,49,53,30,41,2,43,69,66,22,44,68,26,25,40,62,61,33,47,42,48,6,19,17,3,14,31,4,15,46,5,32,16,18,59,57,54,45,55,60,56,58,70,71,72,12,9,28,23,39,10,13,24,27,38,65,67];
end



StepSize=0.1/10;


%% Create model
EltOpts.problem='2dstrain';
EltOpts.nl='elastoplastic';
Types={1 elementPlane EltOpts};
Materials = cell(nElem,3);



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
Div=1;
% while (Div==1)
    Pload=(0:StepSize:1.0); %load
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

dlmwrite(fileRES,res,'delimiter',' ','precision',15);


end

