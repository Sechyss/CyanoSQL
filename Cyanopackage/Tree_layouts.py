from ete3 import TreeStyle


def circular_tree():
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.mode = "c"
    ts.arc_start = -180
    ts.arc_span = 180

    return ts


def color_node(node, color):
    node.img_style["hz_line_color"] = node.img_style["vt_line_color"] = color
    for leaf in node.iter_descendants():
        leaf.img_style["hz_line_color"] = leaf.img_style["vt_line_color"] = color


def color_leaves(leaf, color):
    if leaf.is_leaf():
        leaf.img_style["hz_line_color"] = leaf.img_style["vt_line_color"] = color


