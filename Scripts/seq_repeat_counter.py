#!/usr/bin/env python3

import itertools
import sys

def string_parser(str):
    longest_counts = []
    doubles_dict = {"AT":"LL", "AC":"MM", "AG":"NN", "TC":"OO", "TG":"PP", "CG":"QQ", "GC":"RR",
                    "GT":"SS", "CT":"TT", "GA":"UU", "CA":"VV", "TA":"WW"}
    longest_counts.append(max(len(list(y)) for (c,y) in itertools.groupby(str)))
    for i in doubles_dict.keys():
        new_string = str.replace(i, doubles_dict[i])
        longest_counts.append(max(len(list(y)) for (c, y) in itertools.groupby(new_string))) # multiplies the doubles by 2 to get total characters repeated
    final_max = max(longest_counts)
    return  final_max

if __name__ == "__main__":
    input_fasta = sys.path[0] + "/" + sys.argv[1]
    count_dict ={}
    with open(input_fasta, "r") as f:
        while True:
            next_n_lines = list(itertools.islice(f, 2))
            if not next_n_lines:
                break
            name = next_n_lines[0].strip()
            seq = next_n_lines[1].strip()
            count_dict[name] = string_parser(seq)

    with open (sys.argv[2], "a") as f:
        for key, value in count_dict.items():
            f.write("{}\t{}\n".format(key, value))
