from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.collections import *
import matplotlib.lines as lines
# from metrics import pos2BB

def transformXY(x, y, w, h):
    newX = x - w/2
    newY = y - h/2
    return newX, newY

def draw(pos, sizes, edges=[], ax=None, polygons=[], node_color=(0,0.8,1,0.5), edge_color=(0,0,0,0.2), edge_width=1, title="", nodes_label=[], label_font_size=11, figsize=(7,5), save=None, off_axis=False, nodes_first=True, overlapped_nodes=None, overlapped_color=(255/255, 51/255,51/255,0.5), overlapped_alpha=0.5, border_width=0):
    if(hasattr(pos, "numpy")):
            pos = pos.numpy()
    if(hasattr(sizes, "numpy")):
        sizes = sizes.numpy()
    pos = pos.squeeze()
    sizes = sizes.squeeze()
    
    assert pos.shape[0] == sizes.shape[0], f"pos: {pos.shape}, sizes: {sizes.shape}" 

    N = pos.shape[0]
    ax_was_None = ax == None
    if(ax is None):
        ax = plt.gca()
    if(figsize is not None and ax_was_None):
        plt.gcf().set_size_inches(figsize)

    xmin = np.inf
    xmax = -np.inf
    ymin = np.inf
    ymax = -np.inf

    do_label_nodes = False
    if(type(nodes_label) == list and len(nodes_label) == N):
        do_label_nodes = True

    node_color_as_list = False
    if(type(node_color) == list and len(node_color) == N):
        node_color_as_list = True

    if(nodes_first):
        nodes_zindex = 3
        edges_zindex = 2
    else:
        nodes_zindex = 2
        edges_zindex = 3

    if(overlapped_nodes is None):
        overlapped_nodes = [False for _ in range(N)]
        
    default_node_alpha = 0.9
    line_patches = []
    if edges is not None:
        for (src, tgt) in edges:
            (x1, y1) = pos[src]
            (x2, y2) = pos[tgt]
            (w1, h1) = sizes[src]
            (w2, h2) = sizes[tgt]
            # x1, y1 = transformXY(x1, y1, w1, h1)        
            # x2, y2 = transformXY(x2, y2, w2, h2)        
            # line = lines.Line2D((x1, x2), (y1,y2), linewidth=edge_width, color=edge_color, alpha=edge_color[-1], zorder=edges_zindex)
            # ax.add_line(line)
            line = ((x1, y1),(x2, y2))
            line_patches.append(line)
        ax.add_collection(LineCollection(line_patches, zorder=edges_zindex, linewidths=edge_width, colors=edge_color))

    _patches = []
    for nodeid in range(N):
        (x, y) = pos[nodeid]
        (w, h) = sizes[nodeid]
        x, y = transformXY(x, y, w, h)        
        col = node_color
        node_alpha = default_node_alpha

        if(node_color_as_list):
            col = node_color[nodeid]        
        if(overlapped_nodes[nodeid]):
            col = overlapped_color
            node_alpha = overlapped_alpha

        if(len(polygons) == 0):
            rect = patches.Rectangle((x, y), w, h, facecolor=col, linewidth=border_width, edgecolor='black', alpha=node_alpha, zorder=nodes_zindex)
        else:
            poly = np.expand_dims(pos[nodeid],0)+polygons[nodeid]
            rect = patches.Polygon(poly,facecolor=col, linewidth=border_width, edgecolor='black', alpha=node_alpha)
        # rect = patches.Ellipse((x, y), w, h, facecolor=col, linewidth=0, edgecolor='black', alpha=node_alpha, zorder=nodes_zindex)
        # ax.add_patch(rect)
        _patches.append(rect)
        if(do_label_nodes):
            ax.annotate(nodes_label[nodeid], (x+w/2, y+h/2), fontsize=label_font_size, color="black", va="center", ha="center")
        xmin = min(xmin, x-w)
        xmax= max(xmax, x+w)
        ymin = min(ymin, y-h)
        ymax= max(ymax, y+h)
        # xmin = ymin = 0.45
        # xmax = ymax = 0.55
    ax.add_collection(PatchCollection(_patches, match_original=True, zorder=nodes_zindex))    
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title(title)
    if(off_axis):
        ax.axis("off")
    if(save is not None):
        plt.tight_layout()
        plt.savefig(save, dpi=120)
    if(ax_was_None):
        plt.tight_layout()
        plt.show()
    return ax

# def draw_points(pos, sizes, ax=None, node_color=(0,0.8,1,0.5), title="", nodes_label=[], label_font_size=11, figsize=(7,5), save=None, off_axis=False, overlapped_nodes=None, overlapped_color="red", nodes_alpha=0.5, zoomX=None, zoomY=None, AR=1., title_size=20):
#     if(hasattr(pos, "numpy")):
#             pos = pos.numpy()
#     if(hasattr(sizes, "numpy")):
#         sizes = sizes.numpy()
#     pos = pos.squeeze()
#     sizes = sizes.squeeze()
    
#     assert pos.shape[0] == sizes.shape[0], f"pos: {pos.shape}, sizes: {sizes.shape}" 

#     N = pos.shape[0]
#     ax_was_None = ax == None
#     if(ax is None):
#         ax = plt.gca()
#     if(figsize is not None and ax_was_None):
#         plt.gcf().set_size_inches(figsize)

#     do_label_nodes = False
#     if(type(nodes_label) == list and len(nodes_label) == N):
#         do_label_nodes = True

#     node_color_as_list = False
#     if((type(node_color) == list or type(node_color) == np.ndarray) and len(node_color) == N):
#         node_color_as_list = True

#     if(overlapped_nodes is None):
#         overlapped_nodes = [False for _ in range(N)]
#     else:
#         if(not node_color_as_list):
#             node_color = np.full((N), node_color)
#         node_color[overlapped_nodes==True] = overlapped_color
        

#     col = EllipseCollection(widths=sizes[:,0], heights=sizes[:,1]*AR, angles=np.zeros(N), offsets=pos,transOffset=ax.transData, units="x",facecolors=node_color, alpha=nodes_alpha)
#     ax.add_collection(col)
    
    
#     # _patches = []
#     # for nodeid in range(N):
#     #     (x, y) = pos[nodeid]
#     #     (w, h) = sizes[nodeid]
#     #     x, y = transformXY(x, y, w, h)        
#     #     col = node_color
#     #     node_alpha = 0.9

#     #     if(node_color_as_list):
#     #         col = node_color[nodeid]        
#     #     if(overlapped_nodes[nodeid]):
#     #         col = overlapped_color
#     #         node_alpha = overlapped_alpha

#     #     rect = patches.Rectangle((x, y), w, h, facecolor=col, linewidth=0, edgecolor='black', alpha=node_alpha)
#     #     # ax.add_patch(rect)
#     #     _patches.append(rect)
#     #     if(do_label_nodes):
#     #         ax.annotate(nodes_label[nodeid], (x+w/2, y+h/2), fontsize=label_font_size, color="black", va="center", ha="center")
#     # ax.add_collection(PatchCollection(_patches, match_original=True))

#     # BB = pos2BB(pos, sizes, arnaud_approx=False)
#     np.array(pos)
#     if(zoomX is None):
#         ax.set_xlim(BB[0] - BB[2]/2,BB[0] + BB[2]/2)
#     else:
#         ax.set_xlim(*zoomX)
#     if(zoomY is None):
#         ax.set_ylim(BB[1] - BB[3]/2,BB[1] + BB[3]/2)
#     else:
#         ax.set_ylim(*zoomY)
#     ax.set_title(title, fontsize=title_size)
#     if(off_axis):
#         ax.axis("off")
#     if(save is not None):
#         plt.tight_layout()
#         plt.savefig(save, dpi=1000)
#     if(ax_was_None):
#         plt.tight_layout()
#         plt.show()
#     return ax

def matrix(A, cmap=None, values=False, axis_off=True, precision=1, figsize=None):
    if(figsize is not None):
        plt.gcf().set_size_inches(figsize)
    plt.imshow(A, cmap=cmap, origin="lower")
    if(axis_off):
        plt.axis("off")
    if(values):
        ax = plt.gca()
        for (i, j), z in np.ndenumerate(A):
            prec = '{:0.'+str(precision)+'f}'
            ax.text(j, i, prec.format(z), ha='center', va='center')
    plt.show()
    return

def monitor_pass(stress, overlap, scale, maxScaleRatio=1, figsize=(8,8), save=None, ax=None):
    ax_was_none = ax is None
    if(ax_was_none):
        ax = plt.gca()
    if(figsize is not None and ax_was_none):
        plt.gcf().set_size_inches(figsize)
    
    n_passes = len(stress)
    if(n_passes > 1):
        stress = 100*np.copy(stress) / np.max(stress)
        overlap = 100*np.copy(overlap) / np.max(overlap) 
        scale = 100*np.copy(scale) / maxScaleRatio
    X = np.arange(0, n_passes, 1)
    ax.plot(X, stress, label="stress")
    ax.plot(X, overlap, label="overlaps")
    ax.plot(X, scale, label="scale")
    ax.set_xlabel("number of passes")
    ticks = np.arange(0,110,10)
    ticks_labels = [str(t)+"%" for t in ticks]
    ax.set_yticks(ticks, ticks_labels)
    ax.legend()
    if(save is not None and ax_was_none):
        plt.savefig(save, dpi=1000)
    if(ax_was_none):
        plt.show()
    return ax

def monitor_iter(stress, overlap, scale, pass_length, maxScaleRatio=1, figsize=(6,6), save=None, ax=None, show=True, legend=False, yticks=False, xlabel=False):
    ax_was_none = ax is None
    if(ax_was_none):
        ax = plt.gca()
    if(figsize is not None and ax_was_none):
        plt.gcf().set_size_inches(figsize)
    
    n_passes = len(pass_length)
    n_iter = len(stress)
    if(n_passes > 1):
        stress = 100*np.copy(stress) / np.max(stress)
        overlap = 100*np.copy(overlap) / np.max(overlap) 
        scale = 100*np.copy(scale) / maxScaleRatio
    X = np.arange(0, n_iter, 1)
    ax.plot(X, stress, label="stress")
    ax.plot(X, overlap, label="overlaps")
    ax.step(X, scale, label="upscale", where="post")
    pl_label = []
    pl_tick = []
    old_scale = -1
    pass_num = 0
    lastOk = 0
    text_x_spacing = 5
    for i, s in enumerate(scale):
        if(s != old_scale):
            old_scale = s
            pass_num +=1
            pl_tick.append(i)
            if(len(pl_tick) > 1 and abs(pl_tick[-2] - pl_tick[-1]) <= text_x_spacing):
                pl_label.append("")
                pl_label[lastOk] = pl_label[lastOk] + ", " + str(pass_num)
            else:
                pl_label.append("pass "+str(pass_num))
                lastOk = pass_num-1
    for i, pl in enumerate(pl_label):
        ax.vlines(x=pl_tick[i], ymin=-5, ymax=105, color="black", ls=(0, (1, 5)))
        # ax.text(pl_tick[i]-text_x_spacing, 75, pl, rotation=90, verticalalignment='center')

    if(xlabel):
        ax.set_xlabel("Total number of iterations", fontsize=18)
    ticks = np.arange(0,110,20)
    ticks_labels = [str(t)+"%" for t in ticks]
    if(yticks):
        ax.set_yticks(ticks, ticks_labels)
    else:
        ax.set_yticks(ticks, ["" for t in ticks])
    # ax.set_xticks(np.arange(0, np.sum(pass_length), 10))
    if(legend):
        ax.legend(framealpha=1, prop={"size":18}) # framealpha=1,bbox_to_anchor=(0.07,0.95), handles=lns
    ax.set_ylim(-1, 105)
    ax.set_xlim(-8, np.sum(pass_length)+5)
    if(save is not None):
        plt.tight_layout()
        plt.savefig(save, dpi=1000)
    if(ax_was_none and show):
        plt.show()
    return ax


