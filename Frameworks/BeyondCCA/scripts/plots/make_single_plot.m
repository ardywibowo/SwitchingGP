function make_single_plot(xs, ys, ys_L, ys_U, algopts, plotopts)
%MAKE_SINGLE_PLOT Produces nice plots
%
%
% Copyright 2016, Anastasia Podosinnikova

  ver = version('-release');
  if str2double(ver(1:4)) < 2014
    warning('Use Matlab 2014 or later to reproduce the plots style correctly')
  end

  ff = figure('Visible','On'); hold on

    % the size and position of the figure on the screen
    position = plotopts.position;
    set(ff, 'Position', position)
    % font type for text and ticks
    set(gca, 'Fontname', plotopts.fontname)
    % font size for title, axis, ticks
    set(gca, 'FontSize', plotopts.fontsize)
    % whether latex interpreter for text is used
    if ~isfield(plotopts,'islatex'), plotopts.islatex = 0; end
    % setting title, if any
    if plotopts.istitle == 1
      if plotopts.islatex == 1
        title(plotopts.titlestr,'FontSize',plotopts.fontsize,'Interpreter','latex')
      else
        title(plotopts.titlestr,'FontSize',plotopts.fontsize)
      end
    end
    % x-axis name
    if plotopts.islatex == 1,
        xlabel(plotopts.xname,'Interpreter','latex')
    else
      xlabel(plotopts.xname)
    end
    % y-axis name
    if plotopts.islatex == 1,
      ylabel(plotopts.yname,'Interpreter','latex')
    else 
        ylabel(plotopts.yname)
    end
    % setting x- and y-limits
    if plotopts.isxlim==1, xlim(plotopts.xlims), end
    if plotopts.isylim==1, ylim(plotopts.ylims), end
    % setting x- and y-ticks values, if any
    if plotopts.isxticks==1, set(gca,'XTickLabel',plotopts.xticks); end
    if plotopts.isyticks==1, set(gca,'YTick',plotopts.yticks); end
    % set transparent background
    set(gcf, 'Color', 'None')
    
    % convert yys to the "errorbars format"
    len_x = length(xs);
    nalg = algopts.nalg;
    
    if (size(ys,1) ~= len_x) || (size(ys,2) ~= nalg) || ...
       (size(ys_L,1) ~= len_x) || (size(ys_L,2) ~= nalg) || ...
       (size(ys_U,1) ~= len_x) || (size(ys_U,2) ~= nalg)
      error('wrong size of the input')
    end
    
    % the order of plots
    alginds = 1:nalg;
    
    algspecs  = algopts.algspecs;
    algcolors = algopts.algcolors;
    markers   = algopts.markers;
    markcols  = algopts.markcols;
    marksizes = algopts.marksizes;
    
    linewidth = plotopts.linewidth;
    
    % actual plots
    for ii = 1:nalg
      ialg = alginds(ii);
      % plot the plots
      errorbar(xs,ys(:,ialg),ys_L(:,ialg),ys_U(:,ialg),...
         algspecs{ialg},'Color',algcolors{ialg},'LineWidth',linewidth);
      % compute the positions for the markers
      step_between_markers = 1/(nalg+1);
      mms = ys(:,ialg);
      mmlen = length(mms); 
      x_marker = zeros(mmlen - 1,1);
      y_marker = zeros(mmlen - 1,1); 
      for imm = 1:mmlen-1;
        y_marker(imm) = mms(imm) + (1-step_between_markers*ii)*(mms(imm+1)-mms(imm));
        x_marker(imm) = xs(imm) + (1-step_between_markers*ii)*(xs(imm+1)-xs(imm));
      end
      % add markers
      plot(x_marker,y_marker,'.','LineWidth',linewidth,...
        'Marker',markers{ialg},'MarkerSize',marksizes{ialg},...
        'MarkerFaceColor',markcols{ialg},'MarkerEdgeColor',algcolors{ialg});
    end

    % legend, if any
    if ~isfield(plotopts,'islegend'), plotopts.islegend = 0; end
    if plotopts.islegend == 1   
      ttemp = zeros(nalg,1);
      for ii = alginds
        ialg = alginds(ii);
        ttemp(ii) = plot(xs, Inf*ones(1,length(xs)), ...
           algspecs{ialg},'Color',algcolors{ialg},'LineWidth',linewidth,...
           'Marker',algopts.markers{ialg},'Markersize',marksizes{ialg},...
           'MarkerFaceColor',markcols{ialg},'MarkerEdgeColor',algcolors{ialg});
      end
      algnames = algopts.algnames;
      [leg, legobj] = legend(ttemp,algnames{alginds});
      %textobj = findobj(legobj, 'type', 'text' );
      set( leg, 'fontsize', plotopts.legendsize );
      %set(textobj, 'fontsize', plotopts.legendsize);
      %set(textobj, 'fontname', plotopts.fontname);
      if plotopts.islatex == 1, 
        set(leg,'Interpreter','latex'); 
      end
      %set(leg, 'position', plotopts.legposition);
      %set(leg, 'fontname', plotopts.fontname );
%     set(ll,'interpreter','latex');
%     set(ll,'fontname','timesnewroman');
      set(leg, 'location', 'NorthEast' );
      set(leg,'box','on')
    end
    
    box on
    
    %set(ff,'InvertHardcopy','off','Position', position,'PaperPositionMode','auto')
    %set(ff,'InvertHardcopy','off','Position',get(gcf,'Position'),'PaperPositionMode','auto')
    plotname2save = plotopts.plotname2save;
    saveas(ff,[plotname2save,'.eps'],'epsc')
  hold off
  
end 
