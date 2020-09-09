import re

from .scan import *


def yylex(scanner):
    scan_flag = 1
    for key in switcher:
        if re.match('^' + key + '$', item):
            tokval = switcher[key]
            tokstr = item
            scan_flag = 0
            break
    if scan_flag:
        return -1, "Invalid argument"
    return tokval, tokstr
