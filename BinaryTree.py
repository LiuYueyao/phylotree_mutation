class BinaryTree(object):

    def __init__(self, name_value):
        self.name = name_value  # name
        self.seq = []  # alignment sequence, list
        self.leftchild = None
        self.rightchild = None

    def set_name(self, name_value):
        self.name = name_value

    def set_seq(self, seq_value):
        self.seq = seq_value

    def append_seq(self, char):
        self.seq.append(char)

    def get_seq_i(self,i):
        return self.seq[i]

    def is_leaf(self):
        return (not self.leftchild) and (not self.rightchild)


def in_traversal(tree):
    if tree != None:
        in_traversal(tree.leftchild)
        print(tree.seq)
        in_traversal(tree.rightchild)


def lr_path(root):  # leaf -> root
    re = []
    if not root:
        return re
    if root.is_leaf():
        tmp_re = []
        tmp_re.append(root.seq)
        re.append(tmp_re)
        return re

    for path in lr_path(root.leftchild):
        path.append(root.seq)
        re.append(path)
    for path in lr_path(root.rightchild):
        path.append(root.seq)
        re.append(path)
    return re


def reversed_path(node, root_seq):
    lr_p = lr_path(node)  # leaf to root path, lr_p = [[leaf1,..., root], [lesf2,..., root], ..., [leafn, ..., root]]
    for i in range(len(lr_p)):
        lr_p[i].append(root_seq)
        lr_p[i].reverse()
    return lr_p  # root to leaf
