% Genetic Data Analysis and Visualization
% =========================================================================
% Processes genetic data across species, analyzing gene presence/absence. 
% Supports customizable parameters for targeted analysis, including
% data summarization, visualization, and cross-correlation analysis. 
% Generates detailed plots and exports results to CSV.
%
% Usage:
% - Uncomment certain lines for standalone use in MATLAB.
% - Adjust parameters in the 'Parameters' section as needed in standalone 
%   use.
% - Ensure correct working directory.
%
% Outputs:
% - Visual plots in JPEG format (other formats supported, see Parameters.
% - CSV summary of gene analysis and cross-correlations.
%
% Requirements: MATLAB R2019b+, input data in specified folder structure.
%
% Author: Chris van der Ploeg & Manon de Visser | Date: 2024-03-26 | 
% Contact: vanderploeg.cj@gmail.com
% =========================================================================

% clear all;                    % UNCOMMENT for standalone running in Matlab
close all;  clc;
cd(current_cd);                 % COMMENT for standalone running in Matlab

% First, collect all species folders
total_results = {};
all_directories = dir;

% Parameters
filetype = 'jpg';
% save_file_name = 'all_data';                      % UNCOMMENT for standalone running in Matlab
% sum_for_each_subspecies = true;                   % UNCOMMENT for standalone running in Matlab
% plot_figures =            true;                   % UNCOMMENT for standalone running in Matlab
% contig_width =            100;                    % UNCOMMENT for standalone running in Matlab
% reads_threshold =         10;                     % UNCOMMENT for standalone running in Matlab
% categories =              {'ST','FT','hatch'};    % UNCOMMENT for standalone running in Matlab
% ctrl_category =           'hatch';                % UNCOMMENT for standalone running in Matlab
% common_identifier =       'DN';                   % UNCOMMENT for standalone running in Matlab

for i=1:length(categories)
    if contains(categories{i},ctrl_category)
        ctrl_category_index = i;
    end
end

if plot_figures
    mkdir(strcat('figures_',date))
    cd(strcat('figures_',date))
    figure_cd = cd;
    cd ..
end
count=1;
for i=1:length(all_directories)
    if ~contains(all_directories(i).name,'.')
       species_folder{count} = all_directories(i).name;
       count = count+1;
    end
end
co=1;
% Then, go into each species folder
ctrl_x_correlations = [];   
for species=1:length(species_folder)
    clear gene_file gene_data_unformatted
    cd(species_folder{species})
    % Import all sample names
    categories_tot_i = readcell('individuals.txt','Delimiter','|');
    categories_tot = categories_tot_i(contains(categories_tot_i,categories));
    if sum(contains(categories_tot_i,categories))<length(categories_tot_i)
        disp('I found more categories in individuals.txt than selected in the parameters')
    end
    for i = 1:length(categories)
        if sum(contains(categories_tot,categories(i)))==0
          error('Found less categories in the input files than selected in the parameters')  
        end
    end
    % Create an overview of all genes
    all_files = dir;
    count=1;
    for i=1:length(all_files)
        if contains(all_files(i).name,common_identifier)
           gene_file{count} = all_files(i).name;
           count = count+1;
        end
    end
    if count==1
        disp('Did not find any files here!')
        return;
    end
    % Then, loop for each gene and check for presence/absence
    for gene=1:length(gene_file) 
        clear coverage gene_data_unformatted gene_data_unformatted_i
        cov_count = 1;

        disp(strcat("working on species ",string(species)," of ",string(length(species_folder))))
        disp(strcat("working on gene ",string(gene)," of ",string(length(gene_file))))
        gene_data_unformatted_i = readcell(gene_file{gene},'FileType','text');
        % Check whether all samples are present in the gene data. Otherwise,
        % append with zeros
        sum_category = zeros(length(gene_data_unformatted_i),1);
        for u=1:length(categories_tot)
            if max(contains(string(char(gene_data_unformatted_i{:,4})),(categories_tot{u})))==1
                available(u)=1;
                sum_category = sum_category+contains(string(char(gene_data_unformatted_i{:,4})),(categories_tot{u}));
                continue;
            else
                available(u)=0;
            end
        end
        if sum(sum_category)<length(sum_category)
            disp('There are more categories available than selected in the parameters')
        end
        gene_data_unformatted = gene_data_unformatted_i(find(sum_category>0),:);
        uu_sto = 1;
        u=1;
        [tot_gen,~]=size(gene_data_unformatted);
        gene_length = tot_gen/sum(available);
        if sum(available)<length(available)
            while u<=length(available)-sum(available)
                for uu=1:length(categories_tot)
                    if available(uu)==0
                        for uuu = 1:gene_length
                            gene_data_unformatted{tot_gen+gene_length*(u-1)+uuu,1}=gene_data_unformatted{1,1};
                            gene_data_unformatted{tot_gen+gene_length*(u-1)+uuu,2}=uuu;
                            gene_data_unformatted{tot_gen+gene_length*(u-1)+uuu,3}=0;
                            gene_data_unformatted{tot_gen+gene_length*(u-1)+uuu,4}=categories_tot{uu};
                            tot_gen+gene_length*(u-1)+uuu;
                        end
                        u=u+1;
                    end
                end
            end
        end
        clear gene_data
        if plot_figures
            h1 = figure();
            hold on
        end
        title(strrep(gene_data_unformatted{1,1},'_','\_'))
        gene_data_tot = zeros(gene_length,length(categories));
        gene_data_ctrl = zeros(gene_length,length(categories_tot)/length(categories));
        ctrl_counter = 1;
        subplot_count = 1;
        for i=1:length(categories_tot)
            clear count_contig_v
            count_contig = 0;
            contig_width_counter = 0;
            for j  = 1:gene_length
                count_contig = count_contig+1;
                count_contig_v(j) = count_contig;
                if contig_width_counter==0 && gene_data_unformatted{j+(i-1)*gene_length,3}>0
                    contig_width_counter = contig_width_counter + 1;
                    contig_start_index = j;
                elseif contig_width_counter>0 && gene_data_unformatted{j+(i-1)*gene_length,3}==0
                    if contig_width_counter<contig_width
                        gene_data(contig_start_index:j,i) = 0;
                    else
                        gene_data(j,i)  = gene_data_unformatted{j+(i-1)*gene_length,3};
                    end
                    contig_width_counter = 0;
                else
                    gene_data(j,i)  = gene_data_unformatted{j+(i-1)*gene_length,3};
                    if gene_data(j,i)>0
                        contig_width_counter = contig_width_counter+1;
                    end
                    if j==gene_length && contig_width_counter>0
                        if contig_width_counter<contig_width
                            gene_data(contig_start_index:j,i) = 0;
                        else
                            gene_data(j,i)  = gene_data_unformatted{j+(i-1)*gene_length,3};
                        end
                        contig_width_counter = 0;
                    end
                end
            end
            coverage(cov_count) = max(gene_data(:,i));
            cov_count = cov_count+1;
            if norm(gene_data(:,i),inf)<reads_threshold
                gene_data(:,i) = zeros(length(gene_data(:,i)),1);    
            end
            % only if sequences should be merged
            for q = 1:length(categories)
                if contains(gene_data_unformatted{1+(i-1)*gene_length,4},categories{q}) && sum_for_each_subspecies==true
                    gene_data_tot(:,q) = gene_data_tot(:,q)+gene_data(:,i);
                end
            end
            %
            if contains(gene_data_unformatted{1+(i-1)*gene_length,4},categories{ctrl_category_index})
                gene_data_ctrl(:,ctrl_counter) = gene_data(:,i);
                ctrl_counter = ctrl_counter + 1;
            end
            sample_name{i} = gene_data_unformatted{1+(i-1)*gene_length,4};
            gene_name{i} = gene_data_unformatted{1+(i-1)*gene_length,1};
            if sum_for_each_subspecies==false
                for q = 1:length(categories)
                    contains_name(i,q)=contains(sample_name{i},categories{q});
                end
            else
                contains_name = eye(3);
            end
        end
        % Sort data and plot depth file histogram
        if sum_for_each_subspecies==true
            for i=1:length(categories)
                if max(gene_data_tot(:,i))~=min(gene_data_tot(:,i))
                    gene_data_tot(:,i) = gene_data_tot(:,i)/max(gene_data_tot(:,i));
                end
                if plot_figures 
                    subplot(length(categories),1,subplot_count);
                    hold on
                    plot(gene_data_tot(:,i),'LineWidth',1.5)
                    axis tight
                    ylabel(strcat({'Reads '},categories(i)))
                    if subplot_count==length(categories_tot)
                        xlabel('Target position [-]')
                    end
                    if subplot_count==1
                        title(strrep(gene_data_unformatted{1,1},'_','\_'))
                    end
                    subplot_count = subplot_count+1;
                end
            end
            clear gene_data
            gene_data = gene_data_tot;
        elseif sum_for_each_subspecies==false
            gene_data_tot = gene_data;
            for i=1:length(categories_tot)
                if max(gene_data_tot(:,i))~=min(gene_data_tot(:,i))
                    gene_data_tot(:,i) = gene_data_tot(:,i)/max(gene_data_tot(:,i));
                end
                if plot_figures 
                    subplot(length(categories_tot),1,subplot_count);
                    hold on
                    plot(gene_data_tot(:,i),'LineWidth',1.5)
                    axis tight
                    ylabel(strcat({'Reads '},categories_tot(i)))
                    if subplot_count==length(categories_tot)
                        xlabel('Target position [-]')
                    end
                    if subplot_count==1
                        title(strrep(gene_data_unformatted{1,1},'_','\_'))
                    end
                    subplot_count = subplot_count+1;
                end
            end
            clear gene_data
            gene_data = gene_data_tot;
        end
        [~, size_gen_dat] = size(gene_data);
        % First cross-correlate the control group
        [~, size_cross_ctrl] = size(gene_data_ctrl);
        x_corr_ctrl_counter = 1;
        clear x_corr_seq_ctrl_entry;
        for ctrl_count = 1:ctrl_counter-1
            if ctrl_count~=ctrl_counter
                for ctrl_count_2 = ctrl_count+1:ctrl_counter-1
                    if max(gene_data_ctrl(:,ctrl_count))==min(gene_data_ctrl(:,ctrl_count)) || max(gene_data_ctrl(:,ctrl_count_2))==min(gene_data_ctrl(:,ctrl_count_2))
                        x_corr_seq_ctrl = xcorr(gene_data_ctrl(:,ctrl_count),gene_data_ctrl(:,ctrl_count_2));
                    else
                        x_corr_seq_ctrl = xcorr(gene_data_ctrl(:,ctrl_count),gene_data_ctrl(:,ctrl_count_2),'normalized');
                    end
                    x_corr_seq_ctrl_entry(x_corr_ctrl_counter) = x_corr_seq_ctrl(floor(length(x_corr_seq_ctrl)/2));
                    x_corr_ctrl_counter = x_corr_ctrl_counter + 1;
                end
            end
        end
        
        ctrl_x_correlations = [ctrl_x_correlations;x_corr_seq_ctrl_entry'];
        % Cross-correlate single genes
        if plot_figures
            h2 = figure();
            hold on
            title(strrep(gene_data_unformatted{1,1},'_','\_'))
            ax=gca;
        end
        plotted=[];
        if sum_for_each_subspecies ==true
            categories_n = categories;
        else
            categories_n = categories_tot;
        end
        % Plot cross-correlation figures
        for i=1:length(categories_n)
            count=1;
            for j=1:length(categories_n)
                if max(gene_data(:,i))==min(gene_data(:,i)) || max(gene_data(:,j))==min(gene_data(:,j))
                    x_corr_seq = xcorr(gene_data(:,i),gene_data(:,j));
                else
                    x_corr_seq = xcorr(gene_data(:,i),gene_data(:,j),'normalized');
                end
                if i~=j
                    if sum_for_each_subspecies==true
                        val_corr(count,i) = x_corr_seq(gene_length);
                        combination{i,count} = [find(contains_name(i,:),1),find(contains_name(j,:),1)];
                        [size_plotted,~] = size(plotted);
                        not_plotted = 1;
                        if size_plotted == 0
                            plotted = num2str(combination{i,count});
                            gene_temp = combination{i,count};
                            cat = 1;
                            if plot_figures
                                plot(x_corr_seq,'Color',ax.ColorOrder(cat,:),'DisplayName',strcat(string(categories(gene_temp(1))),'->',string(categories(gene_temp(2))),' xcorr=',string(100*val_corr(count,i)),'%'),'LineWidth',2);
                            end
                            total_results{co,1} = string(species_folder(species));
                            total_results{co,2} = gene_data_unformatted{1,1};
                            total_results{co,3} = string(categories(gene_temp(1)));
                            total_results{co,4} = string(categories(gene_temp(2)));
                            total_results{co,5} = x_corr_seq(floor(length(x_corr_seq)/2));
                            total_results{co,6} = max(coverage);
                            co = co+1;
                        else
                            for q = 1:size_plotted
                                if contains(plotted(q,:),num2str(combination{i,count}))==0 && contains(plotted(q,:),flip(num2str(combination{i,count})))==0
                                    not_plotted = not_plotted+1;
                                else 
                                    cat = q;
                                end
                            end
                            if not_plotted>size_plotted
                                cat = size_plotted+1;
                                gene_temp = combination{i,count};
                                plotted = [plotted;num2str(combination{i,count})];
                                if plot_figures
                                    plot(x_corr_seq,'Color',ax.ColorOrder(cat,:),'DisplayName',strcat(string(categories(gene_temp(1))),'->',string(categories(gene_temp(2))),' xcorr=',string(100*val_corr(count,i)),'%'),'LineWidth',2)
                                end
                                total_results{co,1} = string(species_folder(species));
                                total_results{co,2} = gene_data_unformatted{1,1};
                                total_results{co,3} = string(categories(gene_temp(1)));
                                total_results{co,4} = string(categories(gene_temp(2)));
                                total_results{co,5} = x_corr_seq(floor(length(x_corr_seq)/2));
                                total_results{co,6} = max(coverage);
                                co = co+1;
                            end
                        end
                        count = count+1;
                else
                    val_corr(count,i) = x_corr_seq(gene_length);
                        combination{i,count} = [find(contains_name(i,:),1),find(contains_name(j,:),1)];
                        [size_plotted,~] = size(plotted);
                        not_plotted = 1;
                        if size_plotted == 0
                            if max(contains_name(i,:)+contains_name(j,:))==1
                                plotted = num2str(combination{i,count});
                                gene_temp = combination{i,count};
                                cat = 1;
                                if plot_figures
                                    plot(x_corr_seq,'DisplayName',strcat(string(categories_n(i)),'->',string(categories_n(j)),' xcorr=',string(100*val_corr(count,i)),'%'),'LineWidth',2);
                                end
                                total_results{co,1} = string(species_folder(species));
                                total_results{co,2} = gene_data_unformatted{1,1};
                                total_results{co,3} = string(categories_n(i));
                                total_results{co,4} = string(categories_n(j));
                                total_results{co,5} = x_corr_seq(floor(length(x_corr_seq)/2));
                                total_results{co,6} = max(coverage);
                                co = co+1;
                            end
                        else
                            for q = 1:size_plotted
                                if contains(plotted(q,:),num2str(combination{i,count}))==0 && contains(plotted(q,:),flip(num2str(combination{i,count})))==0
                                    not_plotted = not_plotted+1;
                                else 
                                    cat = q;
                                end
                            end
                            if max(contains_name(i,:)+contains_name(j,:))==1
                                cat = size_plotted+1;
                                gene_temp = combination{i,count};
                                plotted = [plotted;num2str(combination{i,count})];
                                if plot_figures
                                    plot(x_corr_seq,'DisplayName',strcat(string(categories_n(i)),'->',string(categories_n(j)),' xcorr=',string(100*val_corr(count,i)),'%'),'LineWidth',2)
                                end
                                total_results{co,1} = string(species_folder(species));
                                total_results{co,2} = gene_data_unformatted{1,1};
                                total_results{co,3} = string(categories_n(i));
                                total_results{co,4} = string(categories_n(j));
                                total_results{co,5} = x_corr_seq(floor(length(x_corr_seq)/2));
                                total_results{co,6} = max(coverage);
                                co = co+1;
                            end
                        end
                        count = count+1;
                    end
                    
                end
            end
        end
        % Save the depth file histograms and cross-correlation figures
        if plot_figures
            xlabel('Shift in target position [-]')
            ylabel('Normalized cross-correlation [-]')
            legend(gca,'show')
            curr_cd = cd;
            cd(figure_cd)
            file_name_fig = gene_data_unformatted{1,1};
            try
                file_name_fig = erase(file_name_fig,'.Tpyg');
            end
            plot([gene_length gene_length],[0 1],'--k','LineWidth',2)
            axis([0 gene_length*2 0 1.5])
            h1.OuterPosition = [681,   355,   300,   453];
            h1.InnerPosition = [689,   363,   284,   360];
            h2.InnerPosition = [984,   456,   292,   290];
            h2.OuterPosition = [976,   448,   308,   383];
            saveas(h1,strcat(file_name_fig,'_subplots'),filetype)
            saveas(h2,strcat(file_name_fig,'_xcorrs'),filetype)
            cd(curr_cd)
            close(h1) %
            close(h2) %
        end
    end
    cd ..
end
ctrl_plot_data = reshape((1-ctrl_x_correlations).*100,[],1);
ctrl_plot_data(find(ctrl_plot_data==100))=NaN;
% plot the histogram as a 'filetype' file (see Parameters)
if plot_figures
    h3=figure();
    hold on
    histogram((ctrl_plot_data),200)
    saveas(h3,'histogram_ctrl',filetype)
end
% save results to a .csv file
cell2csv(strcat(save_file_name,".csv"),total_results)
csvwrite(strcat(save_file_name,"_ctrl_xcorr",".csv"),ctrl_x_correlations)

function cell2csv(filename,cellArray,delimiter)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(filename,cellArray,delimiter)
%
% filename      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% delimiter = seperating sign, normally:',' (default)
%
% by Sylvain Fiedler, KA, 2004
% modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
if nargin<3
    delimiter = ',';
end

datei = fopen(filename,'w');
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)

        var = eval(['cellArray{z,s}']);

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end

        fprintf(datei,var);

        if s ~= size(cellArray,2)
            fprintf(datei,[delimiter]);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);
end
