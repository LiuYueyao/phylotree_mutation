#!/usr/bin/env python3
# coding: utf-8

import BinaryTree
import xlwt
import csv
import time

seq_len = 17213
align_file = 'data/1519_alignment.fasta'  # 1519
nwk_file = 'data/Newick_1519.nwk'
out_file= 'mutation_1519.csv'


def ternary_deduce(child_1, child_2, child_3):  # this is NOT recursion
    dict = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0, 'y': 0, 'r': 0, 'k': 0, 'm': 0, 's': 0, 'b': 0, 'd': 0, 'h': 0,
            'v': 0, 'n': 0}
    re = ''  # result for one char
    for c1 in child_1:
        dict[c1] += 1 / len(child_1)
    for c2 in child_2:
        dict[c2] += 1 / len(child_2)
    for c3 in child_3:
        dict[c3] += 1 / len(child_3)

    m = max(dict.values())
    for key, value in dict.items():
        if value == m:
            re += key  # return characters with highest possibility

    return re


def deduceRoot(root, i):  # deduce i-th char
    if root.is_leaf():
        return root.get_seq_i(i)
    else:
        l = deduceRoot(root.leftchild, i)
        r = deduceRoot(root.rightchild, i)

        # deduce
        if l==r:
            root.append_seq(l)
            return l
        else:
            dict = {'a':0, 't':0, 'g':0, 'c':0, '-':0, 'y':0, 'r':0, 'k':0, 'm':0, 's':0, 'w':0, 'b':0, 'd':0, 'h':0,
                    'v':0, 'n':0}
            re = ''  # result
            for cl in l:
                dict[cl] += 1/len(l)
            for cr in r:
                dict[cr] += 1/len(r)
            m = max(dict.values())
            for key, value in dict.items():
                if value == m:
                    re += key  # return characters with highest possibility

            root.append_seq(re)
            return re


def search_align_file(header):  # to get sequence
    l_seq = []
    flag = 0
    with open(align_file) as f_align:
        for line in f_align:
            if line.startswith('>'):
                if line[1:11] == header[0:10]: # match GenBank name
                    flag = 1  # sequence found. add to list from next line
                else:
                    if l_seq:
                        break
            elif flag == 1:
                l_seq.append(line.strip())
    return "".join(l_seq)  # return string


def build_tree(f):
    str_nwk = open(f).read().strip(';')  # Newick_110.nwk
    stack_nwk = []
    header = ''

    for i, c in enumerate(str_nwk):
        if c == ',':
            if header:  # end of header, push a leaf
                tmp_node = BinaryTree.BinaryTree(header)
                tmp_node.set_seq(list(search_align_file(header)))  # get sequence
                stack_nwk.append(tmp_node)
                header = ''
        elif c == ')':  # new blank node
            if i == len(str_nwk) - 1:  # the last ), left 3 nodes
                break

            if header:  # push a leaf
                tmp_node = BinaryTree.BinaryTree(header)
                tmp_node.set_seq(list(search_align_file(header)))  # get sequence
                stack_nwk.append(tmp_node)
                header = ''

            # if c.next == ';': break
            # print("before", stack_nwk)
            tmp_node = BinaryTree.BinaryTree(0)
            tmp_node.rightchild = stack_nwk.pop()
            tmp_node.leftchild = stack_nwk.pop()
            stack_nwk.pop()
            stack_nwk.append(tmp_node)
            # print("after",stack_nwk)

        elif c == '(':
            stack_nwk.append('(')

        else:  # header
            header += c

    # 3 nodes left in stack
    # deal with tree nodes
    node_1 = stack_nwk.pop()
    node_2 = stack_nwk.pop()
    node_3 = stack_nwk.pop()

    return node_1, node_2, node_3


def main():

    node_1, node_2, node_3 = build_tree(nwk_file)

    #for sequence[0] to sequence[len], deduce
    for i in range(seq_len):
        deduceRoot(node_1, i)
        deduceRoot(node_2, i)
        deduceRoot(node_3, i)

    root_seq = []  # final result, root sequence
    for i in range(seq_len):
        root_seq.append(ternary_deduce(node_1.get_seq_i(i), node_2.get_seq_i(i), node_3.get_seq_i(i)))

    path = BinaryTree.reversed_path(node_1,root_seq) + BinaryTree.reversed_path(node_2,root_seq) + BinaryTree.reversed_path(node_3,root_seq)

    # output
    with open(out_file, 'w') as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(['path', 'depth', 'position', 'parent', 'child'])
        for i in range(len(path)):  # i- path No.
            for j in range(len(path[i])-1):  # j- depth of parent node, from root to leaf
                for k in range(len(path[i][j])):  # k-th character of alignment
                    if path[i][j][k] != path[i][j+1][k]:
                        writer.writerow([i, j, k, path[i][j][k], path[i][j+1][k]])


if __name__ == '__main__':
    start_time = time.clock()
    main()
    print("time:", time.clock()-start_time)