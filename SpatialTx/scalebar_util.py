def add_scalebar(ax, matchx=True, matchy=True, hidex=True, hidey=True, 
                 loc='lower right', fontsize=10, barcolor='black', barwidth=2, offset=(0.05, 0.05), **kwargs):
    """ Add scalebars to axes
    Adds a set of scale bars to *ax*, matching the size to the ticks of the plot
    and optionally hiding the x and y axes
    - ax : the axis to attach ticks to
    - matchx,matchy : if True, set size of scale bars to spacing between ticks
                    if False, size should be set using sizex and sizey params
    - hidex,hidey : if True, hide x-axis and y-axis of parent
    - loc : location ('lower right', 'lower left', 'upper right', 'upper left')
    - fontsize : font size for labels
    - barcolor : color of bars and text
    - barwidth : width of bars
    - offset : offset from corner as fraction of axes size
    - **kwargs : can include sizex, sizey, labelx, labely for manual specification
    Returns None
    """
    def f(axis):
        l = axis.get_majorticklocs()
        return len(l)>1 and (l[1] - l[0])
    
    sizex = 0
    sizey = 0
    labelx = None
    labely = None
    
    if matchx:
        sizex = f(ax.xaxis)
        labelx = str(sizex) if sizex else None
    else:
        sizex = kwargs.get('sizex', 0)
        labelx = kwargs.get('labelx', None)
        
    if matchy:
        sizey = f(ax.yaxis)
        labely = str(sizey) if sizey else None
    else:
        sizey = kwargs.get('sizey', 0)
        labely = kwargs.get('labely', None)
    
    # Use simple matplotlib elements
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Calculate base position based on location
    if 'right' in loc:
        x_base = xlim[1] - (xlim[1] - xlim[0]) * offset[0]
    else:  # left
        x_base = xlim[0] + (xlim[1] - xlim[0]) * offset[0]
    
    if 'lower' in loc:
        y_base = ylim[0] + (ylim[1] - ylim[0]) * offset[1]
    else:  # upper
        y_base = ylim[1] - (ylim[1] - ylim[0]) * offset[1]
    
    # Adjust positions based on what bars we're drawing
    if sizex and sizey:
        # Both bars - position them in an L shape
        if 'right' in loc:
            x_pos = x_base - sizex
        else:
            x_pos = x_base
            
        if 'lower' in loc:
            y_pos = y_base
            y_vert_pos = y_base
        else:
            y_pos = y_base
            y_vert_pos = y_base - sizey
            
    elif sizex:
        # Only horizontal bar
        if 'right' in loc:
            x_pos = x_base - sizex
        else:
            x_pos = x_base
        y_pos = y_base
        
    elif sizey:
        # Only vertical bar
        x_pos = x_base
        if 'lower' in loc:
            y_vert_pos = y_base
        else:
            y_vert_pos = y_base - sizey
    
    # Draw horizontal bar if requested
    if sizex:
        ax.plot([x_pos, x_pos + sizex], [y_pos, y_pos], 
                color=barcolor, linewidth=barwidth, solid_capstyle='butt')
        
        if labelx:
            # Position label below the bar
            label_y = y_pos - (ylim[1] - ylim[0]) * 0.015
            ax.text(x_pos + sizex/2, label_y, labelx, 
                   ha='center', va='top', fontsize=fontsize, color=barcolor)
    
    # Draw vertical bar if requested
    if sizey:
        ax.plot([x_pos, x_pos], [y_vert_pos, y_vert_pos + sizey], 
                color=barcolor, linewidth=barwidth, solid_capstyle='butt')
        
        if labely:
            # Position label to the left of the bar
            label_x = x_pos - (xlim[1] - xlim[0]) * 0.01
            ax.text(label_x, y_vert_pos + sizey/2, labely, 
                   ha='right', va='center', fontsize=fontsize, color=barcolor, rotation=90)

    if hidex : ax.xaxis.set_visible(False)
    if hidey : ax.yaxis.set_visible(False)
    if hidex and hidey: ax.set_frame_on(False)