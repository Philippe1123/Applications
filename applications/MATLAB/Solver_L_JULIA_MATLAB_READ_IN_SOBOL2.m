
function [] = Solver_L_JULIA_MATLAB_READ_IN_SOBOL2(ProcID)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%ProcID=2
%NumberOfGaussPoints=[9 25 49 81]; % paper
NumberOfGaussPoints=[25    81   121   225   289   441   625   729];
str_folder="/home/philippe/.julia/dev/Applications/applications/SPDE/data/";
str_folder="/home/philippeb/.julia/packages/MultilevelEstimators/l8j9n/applications/SPDE/data/";
str_folder="/vsc-hard-mounts/leuven-user/330/vsc33032/.julia/dev/MultilevelEstimators/applications/SPDE/data/";


str_interm=str_folder+"Interm8/Beam/";

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

% file_sort=strcat(str_elements,"SortVec_L",level,".txt");
% srt=dlmread(file_sort);
% E=E(srt);

switch Order
    case 1
        elementPlane='plane4';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(1));
        EltOpts.Sort_nXi=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25];
    case 2
        elementPlane='plane9';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(2));
        EltOpts.Sort_nXi=[11,17,10,18,12,65,71,64,72,66,2,8,1,9,3,74,80,73,81,75,20,26,19,27,21,31,29,35,33,28,34,36,30,32,13,15,16,14,67,69,70,68,49,47,53,51,46,52,54,48,50,4,6,7,5,58,56,62,60,55,61,63,57,59,76,78,79,77,22,24,25,23,40,38,44,42,37,43,45,39,41];
    case 3
        elementPlane='plane16';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(3));
        EltOpts.Sort_nXi=[85,81,78,82,86,41,37,34,38,42,8,4,1,5,9,52,48,45,49,53,96,92,89,93,97,109,107,103,101,100,102,104,108,110,87,79,80,88,43,35,36,44,21,19,15,13,12,14,16,20,22,10,2,3,11,32,30,26,24,23,25,27,31,33,54,46,47,55,98,90,91,99,120,118,114,112,111,113,115,119,121,105,106,83,84,65,63,61,59,57,56,58,60,62,64,66,39,40,17,18,6,7,28,29,50,51,76,74,72,70,68,67,69,71,73,75,77,94,95,116,117];
    case 4
        elementPlane='plane25';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(4));
        EltOpts.Sort_nXi=[145,141,136,142,146,85,81,76,82,86,10,6,1,7,11,100,96,91,97,101,160,156,151,157,161,209,205,201,199,196,200,202,206,210,149,139,140,150,89,79,80,90,59,55,51,49,46,50,52,56,60,14,4,5,15,74,70,66,64,61,65,67,71,75,104,94,95,105,164,154,155,165,224,220,216,214,211,215,217,221,225,203,204,143,144,119,115,113,111,109,106,110,112,116,148,120,83,84,53,54,8,9,68,69,98,99,134,130,158,126,124,121,125,127,131,159,135,190,191,218,219,179,207,175,173,171,169,197,166,198,170,172,176,208,180,178,177,167,168,118,147,137,138,174,117,107,108,88,87,77,78,58,57,47,48,30,29,27,25,23,21,19,17,16,18,20,22,26,28,13,24,12,2,3,45,44,42,40,38,36,34,32,31,33,35,37,41,43,39,73,72,62,63,103,102,92,93,133,162,152,153,163,194,182,183,195,222,212,213,223,192,188,132,186,184,122,181,123,185,187,189,193,129,128,114];
    case 5
        elementPlane='plane36';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(5));
        EltOpts.Sort_nXi=[199,193,188,194,200,97,91,86,92,98,12,6,1,7,13,114,108,103,109,115,216,210,205,211,217,271,267,261,259,256,260,262,268,272,203,191,192,204,101,89,90,102,67,63,57,55,52,56,58,64,68,16,4,5,17,84,80,74,72,69,73,75,81,85,118,106,107,119,220,208,209,221,288,284,278,276,273,277,279,285,289,265,266,197,198,169,165,163,159,157,154,158,160,164,166,170,95,96,61,62,10,11,78,79,112,113,186,182,180,176,174,171,175,177,181,183,187,214,215,282,283,237,269,233,231,227,225,257,222,258,226,228,232,234,270,238,235,223,224,236,201,189,190,202,167,155,156,168,99,87,88,100,65,53,54,66,33,31,29,27,23,21,19,18,20,22,24,28,30,34,32,14,2,3,15,50,48,46,44,40,38,36,35,37,39,41,45,47,51,49,82,70,71,83,116,104,105,117,184,172,173,185,218,206,207,219,254,274,275,255,286,252,250,248,244,242,240,239,241,243,245,249,251,287,253,263,264,229,230,195,196,161,162,135,133,131,129,127,125,123,121,120,122,124,126,128,130,132,136,134,93,94,59,60,25,26,8,9,42,43,76,77,110,111,152,150,148,146,144,142,140,138,137,139,141,143,145,147,149,153,151,178,179,212,213,280,281,246,247];
    case 6
        elementPlane='plane49';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(6));
        EltOpts.Sort_nXi=[287,281,274,282,288,161,155,148,156,162,14,8,1,9,15,182,176,169,177,183,308,302,295,303,309,375,371,365,361,358,362,366,372,376,291,277,278,292,165,151,152,166,81,77,71,67,64,68,72,78,82,18,4,5,19,102,98,92,88,85,89,93,99,103,186,172,173,187,312,298,299,313,396,392,386,382,379,383,387,393,397,369,370,285,286,249,245,243,239,235,232,236,240,244,246,250,159,160,75,76,12,13,96,97,180,181,270,266,264,260,256,253,257,261,265,267,271,306,307,390,391,419,417,413,411,407,403,401,400,402,404,408,412,414,418,420,377,359,360,378,293,275,276,294,251,233,234,252,167,149,150,168,83,65,66,84,41,39,35,33,29,25,23,22,24,26,30,34,36,40,42,20,2,3,21,62,60,56,54,50,46,44,43,45,47,51,55,57,61,63,104,86,87,105,188,170,171,189,272,254,255,273,314,296,297,315,398,380,381,399,440,438,434,432,428,424,422,421,423,425,429,433,435,439,441,409,410,367,368,283,284,241,242,209,207,203,201,199,197,193,191,190,192,194,198,200,202,204,208,210,157,158,73,74,31,32,10,11,52,53,94,95,178,179,230,228,224,222,220,218,214,212,211,213,215,219,221,223,225,229,231,262,263,304,305,388,389,430,431,415,405,406,416,373,363,364,374,335,333,331,329,327,325,323,321,319,317,316,318,320,322,324,326,328,330,332,334,336,289,279,280,290,247,237,238,248,205,195,196,206,163,153,154,164,125,123,121,119,117,115,113,111,109,107,106,108,110,112,114,116,118,120,122,124,126,79,69,70,80,37,27,28,38,16,6,7,17,58,48,49,59,100,90,91,101,146,144,142,140,138,136,134,132,130,128,127,129,131,133,135,137,139,141,143,145,147,184,174,175,185,226,216,217,227,268,258,259,269,310,300,301,311,356,354,352,350,348,346,344,342,340,338,337,339,341,343,345,347,349,351,353,355,357,394,384,385,395,436,426,427,437];
    case 7
        elementPlane='plane64';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(7));
        EltOpts.Sort_nXi=[391,385,376,386,392,241,235,226,236,242,16,10,1,11,17,266,260,251,261,267,416,410,401,411,417,547,541,535,531,526,532,536,542,548,397,381,382,398,247,231,232,248,147,141,135,131,126,132,136,142,148,22,6,7,23,172,166,160,156,151,157,161,167,173,272,256,257,273,422,406,407,423,572,566,560,556,551,557,561,567,573,539,540,389,390,347,341,339,335,331,326,332,336,340,342,348,239,240,139,140,14,15,164,165,264,265,372,366,364,360,356,351,357,361,365,367,373,414,415,564,565,599,597,591,589,585,581,577,576,578,582,586,590,592,598,600,549,527,528,550,399,377,378,400,349,327,328,350,249,227,228,250,149,127,128,150,49,47,41,39,35,31,27,26,28,32,36,40,42,48,50,24,2,3,25,74,72,66,64,60,56,52,51,53,57,61,65,67,73,75,174,152,153,175,274,252,253,275,374,352,353,375,424,402,403,425,574,552,553,575,624,622,616,614,610,606,602,601,603,607,611,615,617,623,625,587,588,537,538,387,388,337,338,299,297,291,289,287,285,281,277,276,278,282,286,288,290,292,298,300,237,238,137,138,37,38,12,13,62,63,162,163,262,263,324,322,316,314,312,310,306,302,301,303,307,311,313,315,317,323,325,362,363,412,413,562,563,612,613,593,583,584,594,543,533,534,544,449,447,443,441,439,437,435,433,431,427,426,428,432,434,436,438,440,442,444,448,450,393,383,384,394,343,333,334,344,293,283,284,294,243,233,234,244,199,197,193,191,189,187,185,183,181,177,176,178,182,184,186,188,190,192,194,198,200,143,133,134,144,43,33,34,44,18,8,9,19,68,58,59,69,168,158,159,169,224,222,218,216,214,212,210,208,206,202,201,203,207,209,211,213,215,217,219,223,225,268,258,259,269,318,308,309,319,368,358,359,369,418,408,409,419,474,472,468,466,464,462,460,458,456,452,451,453,457,459,461,463,465,467,469,473,475,568,558,559,569,618,608,609,619,595,579,580,596,545,529,530,546,499,497,495,493,491,489,487,485,483,481,479,477,476,478,480,482,484,486,488,490,492,494,496,498,500,445,429,430,446,395,379,380,396,345,329,330,346,295,279,280,296,245,229,230,246,195,179,180,196,145,129,130,146,99,97,95,93,91,89,87,85,83,81,79,77,76,78,80,82,84,86,88,90,92,94,96,98,100,45,29,30,46,20,4,5,21,70,54,55,71,124,122,120,118,116,114,112,110,108,106,104,102,101,103,105,107,109,111,113,115,117,119,121,123,125,170,154,155,171,220,204,205,221,270,254,255,271,320,304,305,321,370,354,355,371,420,404,405,421,470,454,455,471,524,522,520,518,516,514,512,510,508,506,504,502,501,503,505,507,509,511,513,515,517,519,521,523,525,570,554,555,571,620,604,605,621];
    case 8
        elementPlane='plane81';
        EltOpts.nXi=sqrt(NumberOfGaussPoints(8));
        EltOpts.Sort_nXi=[477,469,460,470,478,261,253,244,254,262,18,10,1,11,19,288,280,271,281,289,504,496,487,497,505,645,639,631,627,622,628,632,640,646,483,465,466,484,267,249,250,268,159,153,145,141,136,142,146,154,160,24,6,7,25,186,180,172,168,163,169,173,181,187,294,276,277,295,510,492,493,511,672,666,658,654,649,655,659,667,673,637,638,475,476,429,423,421,415,411,406,412,416,422,424,430,259,260,151,152,16,17,178,179,286,287,456,450,448,442,438,433,439,443,449,451,457,502,503,664,665,701,699,693,691,685,681,677,676,678,682,686,692,694,700,702,647,623,624,648,485,461,462,486,431,407,408,432,269,245,246,270,161,137,138,162,53,51,45,43,37,33,29,28,30,34,38,44,46,52,54,26,2,3,27,80,78,72,70,64,60,56,55,57,61,65,71,73,79,81,188,164,165,189,296,272,273,297,458,434,435,459,512,488,489,513,674,650,651,675,728,726,720,718,712,708,704,703,705,709,713,719,721,727,729,687,688,633,634,471,472,417,418,323,321,315,313,309,307,303,299,298,300,304,308,310,314,316,322,324,255,256,147,148,39,40,12,13,66,67,174,175,282,283,350,348,342,340,336,334,330,326,325,327,331,335,337,341,343,349,351,444,445,498,499,660,661,714,715,695,683,684,696,641,629,630,642,539,537,533,531,529,525,523,521,519,515,514,516,520,522,524,526,530,532,534,538,540,479,467,468,480,425,413,414,426,317,305,306,318,263,251,252,264,215,213,209,207,205,201,199,197,195,191,190,192,196,198,200,202,206,208,210,214,216,155,143,144,156,47,35,36,48,20,8,9,21,74,62,63,75,182,170,171,183,242,240,236,234,232,228,226,224,222,218,217,219,223,225,227,229,233,235,237,241,243,290,278,279,291,344,332,333,345,452,440,441,453,506,494,495,507,566,564,560,558,556,552,550,548,546,542,541,543,547,549,551,553,557,559,561,565,567,668,656,657,669,722,710,711,723,697,679,680,698,643,625,626,644,593,591,589,587,585,583,579,577,575,573,571,569,568,570,572,574,576,578,580,584,586,588,590,592,594,535,517,518,536,481,463,464,482,427,409,410,428,319,301,302,320,265,247,248,266,211,193,194,212,157,139,140,158,107,105,103,101,99,97,93,91,89,87,85,83,82,84,86,88,90,92,94,98,100,102,104,106,108,49,31,32,50,22,4,5,23,76,58,59,77,134,132,130,128,126,124,120,118,116,114,112,110,109,111,113,115,117,119,121,125,127,129,131,133,135,184,166,167,185,238,220,221,239,292,274,275,293,346,328,329,347,454,436,437,455,508,490,491,509,562,544,545,563,620,618,616,614,612,610,606,604,602,600,598,596,595,597,599,601,603,605,607,611,613,615,617,619,621,670,652,653,671,724,706,707,725,689,690,635,636,581,582,527,528,473,474,419,420,377,375,373,371,369,367,365,363,361,359,357,355,353,352,354,356,358,360,362,364,366,368,370,372,374,376,378,311,312,257,258,203,204,149,150,95,96,41,42,14,15,68,69,122,123,176,177,230,231,284,285,338,339,404,402,400,398,396,394,392,390,388,386,384,382,380,379,381,383,385,387,389,391,393,395,397,399,401,403,405,446,447,500,501,554,555,608,609,662,663,716,717];
    case 9
        elementPlane='plane100';
%not used at the moment
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



yposNodes=uniquetol(sort(Nodes(:,3)),10^-5);
% yposNodes=yposNodes([2,4,6]);%%%%%%%%%%%%debug
% ary=[0:2*pi/(length(yposNodes)-1):2*pi];
seldof=[];
vec=[];


for id=1:length(yposNodes)

forceNodes=selectnode(Nodes,min(Nodes(:,2))-1e-6,yposNodes(id)-1e-6,-inf,max(Nodes(:,2))+1e-6,yposNodes(id)+1e-6,inf);
forceNodes=forceNodes(:,1);
[posy,srt]=sort(Nodes(forceNodes,2),'ascend');
forceNodes=forceNodes(srt);
ar=[0:2*pi/(length(forceNodes)-1):2*pi];
% vec_int=(-cos(ar)+1)/sum((-cos(ar)+1)).*(-cos(ary(id))+1)/sum((-cos(ary)+1));
vec_int=(-cos(ar)+1)/sum((-cos(ar)+1)).*1/length(yposNodes);



forceNodes=forceNodes+0.02;
seldof_int=forceNodes';

seldof=[seldof seldof_int];
vec=[vec vec_int];

end


%% Load

Force=10000000*(0:1:1);


    PLoad=Force;
% % % % % % % % ar=[0:2*pi/(length(forceNodes)-1):2*pi];
% % % % % % % % vec=(-cos(ar)+1)/sum((-cos(ar)+1));
P = nodalvalues(DOF,seldof,vec);

%% Compute displacements
U=solver_nr(Nodes,Elements,Types,Sections,Materials,DOF,P,PLoad,Options,[],[],[],[],[]);
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
    

x=-1:0.2:1;
[X,Y]=meshgrid(x,x);
pt=[reshape(X,1,size(X,2)*size(X,1));reshape(Y,1,size(Y,2)*size(Y,1));ones(1,size(Y,2)*size(Y,1));2*ones(1,size(Y,2)*size(Y,1))];


TriangleCoordLocal=[[-1;-1;1;2],[1;-1;1;2],[1;1;1;1],[-1;1;1;2]];

pts=[];
Res=[];
for idx=1:nElem
    
  TriangleVertices=Nodes(Elements(idx,5:8),2:3);
  TriangleVertices=TriangleVertices';
  TriangleVertices=[TriangleVertices;[1 1 1 1];[2 2 2 2]];  
  TransformationMatrix=TriangleVertices*inv(TriangleCoordLocal);
  pt_global=TransformationMatrix*pt;

 pts=[pts; [pt_global(1:2,:)' ones(length( pt_global(1:2,:)'),1).*idx]]; 
end 
 
 %for jx=1:length(pt)
 
  idx=find(abs(pts(:,1)-2.5) < 0.001);
 pts_int=pts(idx,:);  
 midelem=pts_int(find(abs(pts_int(:,2)-0.5) < 0.001),end);

 Nod=Elements(midelem,5:end);
 
  for jx=1:length(pt)

 
 switch Order
    case 1
 Nwres(jx)=sh_qs4(pt(1,jx),pt(2,jx))'*selectdof(DOF,Nod+0.02)*U(:,end); 
     String="Order 1";
    case 2
 Nwres(jx)=sh_qs9(pt(1,jx),pt(2,jx))'*selectdof(DOF,Nod+0.02)*U(:,end); 
      String="Order 2";

    case 3
 Nwres(jx)=sh_qs16(pt(1,jx),pt(2,jx))'*selectdof(DOF,Nod+0.02)*U(:,end); 
       String="Order 3";
    case 4
 Nwres(jx)=sh_qs25(pt(1,jx),pt(2,jx))'*selectdof(DOF,Nod+0.02)*U(:,end);
       String="Order 4";
    case 5
 Nwres(jx)=sh_qs36(pt(1,jx),pt(2,jx))'*selectdof(DOF,Nod+0.02)*U(:,end); 
       String="Order 5";
        
end
 
  end
 
 
  Res=[Res;  Nwres'];
%end 
pts=pts(find(pts(:,3)==midelem),1:2);
  idx=find(abs(pts(:,1)-2.5) < 0.001);
out=[pts(idx,2) Res(idx)];
idx=find(abs(out(:,1)-0.5) < 0.001);
u_out=out(idx,2);
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




clear all;
clc;
close all;

end


