%program to calculate kinship structures
% this program:
% reads in qx and fx files
% constructs U and F matrices
% calls kinship function

%load qx and fx files from csv files
%qx = probability of death in age class x
%fx = age-specific fertility in age class x

qxyears=dlmread('qx_years.csv');
fxyears=dlmread('fx_years.csv');

%pick a year; here, choose year 1, which is 1947
qx=qxyears(:,1);
fx=fxyears(:,1);

%divide fertility by 2 to get female offspring
fx=fx/2;

%construct U matrix and F matrix
px=1-qx(:,1);

%remove last age class (qx=1, px=0)

U=diag(px(1:end-1),-1);
F=zeros(size(U));
F(1,:)=(px.*fx(:,1))'; %post-breeding birth pulse model

%call kinship program
kin_out=kinship_function_DR(U,F);

%names of kin
kinname={'daughter',...
    'granddaughter', ...
    'great-granddaughter', ...
    'mother', ...
    'grandmother', ...
    'great-grandmother', ...
    'older sisters', ...
    'younger sisters', ...
    'nieces (older)', ...
    'nieces (younger)', ...
    'aunts (older)', ...
    'aunts (younger)', ...
    'cousins (older)', ...
    'cousins (younger)'};


%names for combined kin (combining age; e.g., older and younger sisters
%become sisters)
kinname2={'daughters',...
    'granddaughters', ...
    'great-granddaughters', ...
    'mothers', ...
    'grandmothers', ...
    'great-grandmothers', ...
    'sisters', ...
    'nieces', ...
    'aunts', ...
    'cousins'};

%following code gives examples of calculation of various aspects of kin,
%from the output structure arrays

%for purposes of illustration, figures are created to plot the results for
%sisters (combining older and younger sisters), which are kin type 7.

om=kin_out.om;  %maximum age class

%age distributions of kin
% choose an age for Focal
ageFocal = 30;

for ikin=1:10 
    y_age(:,ikin) = (kin_out.allkin2(1:om,ageFocal,ikin))';
end %for ikin

figure
plot(y_age(:,7))
xlabel('Age of kin')
ylabel('Number of kin')
title('Sisters')


%numbers of each kin by age of Focal
for ikin=1:10
    y_numbers(:,ikin) = ( sum(kin_out.allkin2(1:om,:,ikin),1) )';
end % for ikin

figure
plot(y_numbers(:,7))
xlabel('Age of Focal')
ylabel('Number of kin')
title('Sisters')


% deaths of kin
% deaths will be experienced or cumulative, depending on setting in the
% kinship function program

for ikin=1:10
    %number of deaths
    y_deaths(:,ikin) = ( sum(kin_out.allkin2(om+1:end,:,ikin),1) )';
end %for ikin

figure
plot(y_deaths(:,7))
xlabel('Age of Focal')
ylabel('Number of deaths')
title('Sisters')


%dementia prevalence

% age-specific dementia prevalence
prev=dlmread('dementia2015.csv');

for ikin=1:10
    
    %number of relatives with dementia
    y_dementia(:,ikin) = ( prev'*kin_out.allkin2(1:om,:,ikin) )';
    
end %for ikin
 
figure
plot(y_dementia(:,7))
xlabel('Age of Focal')
ylabel('Number of kin with dementia')
title('Sisters')

