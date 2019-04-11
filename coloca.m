%% COLOCA.m
%%% Computing biological colocalization using distance to
%%% connected component iso-contour binning
%%% EV 2019


clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%USER PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

background_ratio=[]; %% A number between 0 and 1, leave blank for automatic
pixel_threshold=5; %% A positive integer to remove small connected components
max_range=40; %% Maximum number of pixel distance to compute in histogram
um_step=0.2; %% Histogram um binsize
max_um=2.2; %% Maximum localization distance to draw in histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[file,path]=uigetfile('Select images to analyze','MultiSelect','on','s','*');


if ~iscell(file)
    tmp=file;
    clear file
    file{1}=tmp;
end

for i=1:length(file)
    I=imread([path file{i}]);
    channels=input('Which channels do you want to analyze for colocalization? Colocalization will be done with respect to first channel selection,  [R=1,G=2,B=3,W=4] ex. [1 2], [2 1 3]...');
    scale=input('Enter scale (px/um): ');
    %%% ROI EXCLUSION %%%
    do_exclusion=input('Do you want to exclude areas?[y/n]','s');
    if strcmpi(do_exclusion,'y');
        exclusion_done='n';
        while strcmpi(exclusion_done,'n')
            figure(1)
            imshow(I)
            title('Please draw polygon of areas to EXCLUDE (double click final point to finish)')
            polyexclude=roipoly;
            I=I.*repmat(uint16(~polyexclude),[1 1 3]);
            exclusion_done=input('Are you finished excluding areas?[y/n]','s');
        end
    end
    %%% ROI INCLUSION %%%
    inclusion_done='n';
    while strcmpi(inclusion_done,'n')
        figure(2)
        imshow(I)
        title('Please draw polygon of areas to INCLUDE (double click final point to finish)')
        [x2,y2,polyinclude,xi2,yi2]=roipoly;
        Icandidate=I.*repmat(uint16(polyinclude),[1 1 3]);
        inclusion_done=input('Are you finished INCLUDING areas?[y/n]','s');
    end
    I=Icandidate(round(min(yi2)):round(max(yi2)),round(min(xi2)):round(max(xi2)),:);
    
    source_channel=channels(1);
    target_channels=channels(2:end);
    
    source=I;
    source(:,:,(1:size(I,3))~=channels(1))=0;
    
    for t=1:length(target_channels)
        target{t}=I;
        target{t}(:,:,(1:size(I,3))~=target_channels(t))=0;
    end
    
    figure
    subplot(length(channels),1,1)
    imshow(source)
    title([file{i} ':Source channel'])
    for t=1:length(target_channels)
        subplot(length(channels),1,1+t)
        imshow(target{t})
        title([file{i} ':Target channel'])
    end
    %%% DENOISING %%%
    if isempty(background_ratio)
        source_threshold=double(max(source(:)))*graythresh(double(source(source>0))/double(max(source(:))));
        for t=1:length(target_channels)
            target_threshold{t}=double(max(target{t}(:)))*graythresh(double(target{t}(target{t}>0))/double(max(target{t}(:))));
        end
    else
        source_threshold=double(max(source(:)))*background_ratio;
        for t=1:legth(target_channels)
            target_threshold{t}=double(max(target{t}(:)))*background_ratio;
        end
    end
    source_foreground=source.*uint16(source>source_threshold);
    for t=1:length(target_channels)
        target_foreground{t}=target{t}.*uint16(target{t}>target_threshold{t});
    end
    figure
    subplot(length(channels),1,1)
    imshow(source_foreground)
    title([file{i} ':Source denoised'])
    for t=1:length(target_channels)
        subplot(length(channels),1,1+t)
        imshow(target_foreground{t})
        title([file{i} ':Target denoised'])
    end
    
    %%% REMOVE SMALL COMPONENTS %%%
    
    source_connected_components = source_foreground.*uint16(bwareaopen(max(source_foreground,[],3),pixel_threshold));
    for t=1:length(target_channels)
        target_connected_components{t} = target_foreground{t}.*uint16(bwareaopen(max(target_foreground{t},[],3),pixel_threshold));
    end
    
    figure
    subplot(length(channels),1,1)
    imshow(source_connected_components)
    title([file{i} ': Source connected components'])
    for t=1:length(target_channels)
        subplot(length(channels),1,1+t)
        imshow(target_connected_components{t})
        title([file{i} ': Target connected components'])
    end
    
    figure
    merged=source_connected_components;
    for t=1:length(target_channels)
        merged=merged+target_connected_components{t};
    end
    imshow(merged)
    title('Merged image of connected components');
    
    %%% COMPUTE DISTANCE BASED ISO-CONTOURS RELATIVE TO TARGET %%%
    
    bw=uint16(max(source_connected_components,[],3)>0);
    zero_contour=bwperim(bw);
    positive_contour{1}=imdilate(bw,ones(3))& ~bw;
    for j=2:max_range
        positive_contour{j}=imdilate(bw,ones(j*2+1))-imdilate(bw,ones((j-1)*2+1));
    end
    negative_contour{1}=~imerode(bw,ones(3))& bw;
    for j=2:max_range
        negative_contour{j}=~imerode(bw,ones(j*2+1))& imerode(bw,ones((j-1)*2+1));
    end
    contour_image=zeros(size(bw));
    for j=1:max_range
        contour_image=contour_image+j*double(positive_contour{j});
        contour_image=contour_image-j*double(negative_contour{j});
    end
    %     contour_image(find(zero_contour))=0;
    
    figure
    imagesc(contour_image)
    colorbar
    title([file{i} ': Distance to nearest connected component'])
    
    %%% COMPUTE HISTOGRAM OF COLOCALIZATION %%%
    range=-max_range:max_range;
    for j=1:length(range)
        for t=1:length(target_channels)
            bins{t}(j,1)=range(j);
            bins{t}(j,2)=sum(sum(double(max(target_connected_components{t},[],3)>0).*(contour_image==range(j))));
        end
    end
    
    um_range=-max_um:um_step:max_um;
    
    for j=1:length(um_range)
        for t=1:length(target_channels)
            um_bins{t}(j,1)=um_range(j);
            try
                um_bins{t}(j,2)=sum(bins{t}(find(and(bins{t}(:,1)/scale<=um_bins{t}(j,1),bins{t}(:,1)/scale>um_bins{t}(j-1,1))),2));
            end
        end
    end
    colors={'r','g','b'};
    for t=1:length(target_channels)
        figure
        bar(um_bins{t}(:,1),um_bins{t}(:,2)./sum(um_bins{t}(:,2)),colors{target_channels(t)})
        title([file{i} '- Target channel ' num2str(target_channels(t)) ' : Histogram of distance of target image pixels to source image components ']);
        xlabel('Distance (um)');
        ylabel('Proportion of pixels [0-1]');
        grid on
    end
    
    %%% DISPLAY COLOCALIZATION PERCENTAGE
    for t=1:length(target_channels)
        display([file{i} ' has ' num2str(sum(bins{t}(bins{t}(:,1)<=0,2))/sum(bins{t}(:,2))*100) '% colocalization of channel ' num2str(target_channels(t)) ' in channel ' num2str(channels(1))]);
    end
    
    %%% WRITE COLOCALIZATION PERCENTAGE & HISTOGRAM BINS TO EXCEL FILE
    try
        existing_data=table2cell(readtable('colocalization.xls'));
    end
    
    headers=['Image','Source Channel','Target Channel','Localization %',num2cell(um_range)];
    data=[];
    for t=1:length(target_channels)
        tmp=[file{i},num2cell(source_channel),num2cell(target_channels(t)),num2cell(sum(bins{t}(bins{t}(:,1)<=0,2))/sum(bins{t}(:,2))*100),num2cell(um_bins{t}(:,2)')];
        data=[data;tmp];
    end
    try
        T=table([headers;existing_data;data]);
    catch
        T=table([headers;data]);
    end
    writetable(T,'colocalization.xls','WriteVariableNames',false);
end