function [ ] = Slope_Slover_JULIA_MATLAB_GlobalNested(ProcID)
%ProcID=2

NumberOfGaussPoints=[16 19 28 37 61 73 126]; %Paper
%NumberOfGaussPoints=[1 7 16 28];
%%%%%%%%%SELECT CORRECT SYSTEM
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";
%str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
%str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
%%%%%%%%%%%%%55


str_interm=str_folder+"Interm3/Slope/";


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
 %   case 1
 %       elementPlane='plane3';
 %       EltOpts.nXi=NumberOfGaussPoints(1);
 %       EltOpts.Sort_nXi=[6,11,3,8,4,13,7,2,9,1,10,12,5];
    case 1
        elementPlane='plane6';
         EltOpts.nXi=NumberOfGaussPoints(1);
         EltOpts.Sort_nXi=[7,16,10,15,2,6,3,9,11,1,13,4,12,8,14,5];

    case 2
        elementPlane='plane10';
        EltOpts.nXi=NumberOfGaussPoints(2);
        EltOpts.Sort_nXi=[13,17,10,19,5,16,6,9,12,1,18,7,14,8,15,11,2,3,4];

    case 3
        elementPlane='plane15';
        EltOpts.nXi=NumberOfGaussPoints(3);
        EltOpts.Sort_nXi=[4,26,7,22,14,25,15,6,3,1,21,16,17,5,18,2,8,9,10,28,20,19,13,12,27,11,23,24];
        
    case 4
        elementPlane='plane21';
        EltOpts.nXi=NumberOfGaussPoints(4);
        EltOpts.Sort_nXi=[4,21,13,30,5,20,6,12,3,1,31,7,28,34,29,2,14,15,16,24,27,32,10,9,25,8,22,23,33,17,36,37,18,26,35,11,19];
        
    case 5
        elementPlane='plane28';
        EltOpts.nXi=NumberOfGaussPoints(5);
        EltOpts.Sort_nXi=[9,4,22,34,19,5,40,23,8,45,36,41,52,53,51,61,10,43,42,27,2,11,29,30,26,44,57,56,20,1,31,21,46,32,47,48,12,7,18,16,13,25,37,28,14,6,38,15,17,24,35,39,50,33,3,55,49,54,58,59,60];
        
    case 6
        elementPlane='plane36t';
        EltOpts.nXi=NumberOfGaussPoints(6);
        EltOpts.Sort_nXi=[8,1,31,39,22,2,50,32,7,42,38,49,64,60,63,73,11,52,53,30,14,15,37,36,29,51,69,68,26,4,41,25,46,40,54,55,58,6,19,17,20,35,44,33,21,5,47,16,18,34,59,48,57,43,45,62,56,61,70,71,72,12,9,23,28,3,10,13,24,27,67,65,66];
        
    case 7
        elementPlane='plane45';
        EltOpts.nXi=NumberOfGaussPoints(7);
        EltOpts.Sort_nXi=[29,61,17,70,4,66,6,18,30,10,68,5,69,16,71,28,1,3,2,64,73,78,14,15,62,13,63,65,109,102,112,114,40,110,111,113,41,49,52,82,103,106,76,7,108,54,9,80,50,104,39,8,75,74,38,105,77,107,81,51,53,85,67,20,88,97,72,90,21,86,87,19,89,55,58,26,79,31,23,43,91,37,34,115,42,118,94,46,48,121,36,96,124,126,100,11,120,84,60,12,24,27,56,116,99,122,123,92,32,125,44,98,45,93,117,101,119,33,95,47,22,35,83,25,57,59];
        
  %  case 9
  %      elementPlane='plane55';
  %      EltOpts.nXi=NumberOfGaussPoints(9);
  %      EltOpts.Sort_nXi=[113,155,21,8,11,157,115,10,153,9,100,102,112,140,22,20,2,4,3,59,68,134,145,57,69,156,72,60,38,137,15,16,14,136,135,154,138,50,71,53,158,43,139,163,161,41,55,73,13,12,51,40,160,159,39,162,96,42,52,70,54,89,98,24,152,101,103,25,87,99,23,88,58,36,56,110,74,27,62,143,92,146,77,164,65,97,149,167,18,67,170,79,151,95,173,175,169,61,19,28,37,111,141,165,171,172,94,147,75,174,63,17,166,148,93,64,168,76,150,142,66,144,26,78,114,35,44,33,86,47,30,104,128,107,80,116,131,122,127,5,119,85,83,121,133,109,49,125,91,1,7,31,34,124,45,105,129,117,82,81,118,6,123,126,130,120,84,106,132,108,29,46,90,32,48];
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

